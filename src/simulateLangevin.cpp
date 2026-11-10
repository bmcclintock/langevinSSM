#include <cmath>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

#define IS_RCPP_BUILD
#include "TMB/include/raster_helpers.hpp"

// Multivariate normal with given standard normal vector z (no RNG use; safe inside threads).
static arma::vec rmvnorm_z(const arma::vec& mean, const arma::mat& sigma, const arma::vec& z) {
  arma::mat chol_sigma;

  // Try Cholesky Decomposition first
  if (arma::chol(chol_sigma, sigma)) {
    return mean + arma::trans(chol_sigma) * z;
  }

  // If Cholesky fails due to floating-point truncation,
  // safely compute using Eigen Decomposition (robust but slower)
  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, sigma);

  // Clamp any microscopic negative eigenvalues to exactly 0
  eigval.elem(arma::find(eigval < 0)).zeros();

  return mean + eigvec * arma::diagmat(arma::sqrt(eigval)) * z;
}

// Helper function for multivariate normal sampling (Hybrid Approach)
arma::vec rmvnorm(const arma::vec& mean, const arma::mat& sigma) {
  int n = mean.n_elem;

  // Generate standard normal random vector using R's RNG
  arma::vec z(n);
  for(int i = 0; i < n; i++) {
    z(i) = R::rnorm(0, 1);
  }
  return rmvnorm_z(mean, sigma, z);
}

// Helper function for multivariate t-distribution sampling
arma::vec rmvt(const arma::vec& mean, const arma::mat& sigma, double df) {
  arma::vec y = rmvnorm(arma::zeros(mean.n_elem), sigma);
  double u = R::rchisq(df);
  return mean + y * std::sqrt(df / u);
}

// Helper function for bivariate pseudo-Huber (hyperbolic) sampling with threshold k:
// density proportional to exp(-k sqrt(k^2 + d^2)), d^2 = r' Sigma^-1 r. The Mahalanobis radius is
// d = k sqrt(s^2 - 1), where s - 1 is a mixture of Exp(k^2) (weight k^2/(k^2+1)) and Gamma(2, k^2);
// the direction is uniform.
arma::vec rmvhyp(const arma::vec& mean, const arma::mat& sigma, double k) {
  double k2 = k * k;
  double t = (R::unif_rand() < k2 / (k2 + 1.0)) ? R::rexp(1.0 / k2) : R::rgamma(2.0, 1.0 / k2);
  double s = 1.0 + t;
  double d = k * std::sqrt(s * s - 1.0);
  double ph = R::runif(0.0, 2.0 * M_PI);
  arma::vec z = {d * std::cos(ph), d * std::sin(ph)};
  arma::mat L = arma::chol(sigma, "lower");
  return mean + L * z;
}

// ---------------------------------------------------------------------------
// Track simulator.
//
// Animals are simulated independently, so with OpenMP (ncores > 1) they are distributed over threads.
// R's RNG is not thread safe, and the R API must not be called from a thread.  Therefore:
//   * all standard normal draws of an animal are made serially from R's RNG before the parallel loop,
//     in exactly the order of the serial algorithm (per animal: 2 initial-velocity draws, then 2
//     (overdamped) or 4 (underdamped) draws per step; animals in order).  The parallel loop only runs
//     the deterministic recursion, so the result is bit-identical for any ncores, and identical to
//     the serial simulator of previous versions for the same seed;
//   * the parallel loop works on plain C++ buffers (no Rcpp objects are created or modified), and an
//     error (position out of bounds, failed decomposition) is stored and re-thrown after the loop, for
//     the first failing animal.
// Animals are processed in batches of ncores, so at most ncores animals' draws are held at a time.
// ---------------------------------------------------------------------------
namespace {

struct SimTrack {
  const double* dt;    // time steps (length nsteps)
  int nsteps;
  const double* pen;   // barrier penalty: length 1 (all steps) or nsteps
  int npen;
  double x0, y0;       // initial position
  size_t start;        // first output row
};

struct SimOut {
  std::vector<double> id, time, dt, mux, muy, velx, vely;
  explicit SimOut(size_t n) : id(n), time(n), dt(n), mux(n), muy(n), velx(n), vely(n) {}
};

struct SimEnv {
  int model;
  double gamma, sigma;
  std::vector<double> beta;
  langevin_raster::SplinePtr cov_spline, barrier_spline;
  const NumericVector* all_z_values;
  const IntegerVector* n_zvals_cov;
  const IntegerVector* cov_offset;
  int n_covs;
  double xmin, xmax, ymin, ymax;
};

inline int draws_per_step(int model) { return model == 1 ? 4 : 2; }

std::string bounds_message(double x, double y, int animal_id, double time, const SimEnv& E) {
  return "Position out of bounds for animal " + std::to_string(animal_id) +
    " at time " + std::to_string(time) +
    "\nPosition: (" + std::to_string(x) + ", " + std::to_string(y) + ")" +
    "\nBounds: x[" + std::to_string(E.xmin) + ", " + std::to_string(E.xmax) + "]" +
    " y[" + std::to_string(E.ymin) + ", " + std::to_string(E.ymax) + "]";
}
inline bool out_of_bounds(double x, double y, const SimEnv& E) {
  return x < E.xmin || x > E.xmax || y < E.ymin || y > E.ymax;
}

// Simulate one animal from its pre-drawn standard normals z.  Thread safe (no R API, no RNG, writes
// only rows [tr.start, tr.start + nsteps] of the output).  Returns "" or an error message.
std::string sim_one_track(int i, const SimTrack& tr, const double* z, const SimEnv& E, SimOut& O) {
  const size_t s0 = tr.start;
  const double sigma = E.sigma, gamma = E.gamma, s2 = sigma * sigma;
  const double sd_v0 = sigma / sqrt(2. * gamma);
  size_t iz = 0;

  // Initialize first observation (R::rnorm(mu, sd) is mu + sd * norm_rand())
  O.id[s0] = i + 1;
  O.time[s0] = 0;
  O.dt[s0] = 0;
  O.mux[s0] = tr.x0;
  O.muy[s0] = tr.y0;
  O.velx[s0] = 0.0 + sd_v0 * z[iz++];
  O.vely[s0] = 0.0 + sd_v0 * z[iz++];

  for (int t = 1; t <= tr.nsteps; t++) {
    const size_t idx = s0 + t;
    const double dt_step = tr.dt[t - 1];

    O.id[idx] = i + 1;
    O.time[idx] = O.time[idx - 1] + dt_step;
    O.dt[idx] = dt_step;

    // Calculate gradients using previous position AND previous time
    arma::mat grad = extract_raster_values(O.mux[idx - 1], O.muy[idx - 1], O.time[idx - 1],
                                           E.cov_spline, *E.all_z_values, *E.n_zvals_cov,
                                           *E.cov_offset, E.n_covs);
    // Calculate force vector
    double h0 = 0.0, h1 = 0.0;
    for (int c = 0; c < E.n_covs; c++) {
      h0 += E.beta[c] * grad(0, c);
      h1 += E.beta[c] * grad(1, c);
    }

    // --- BARRIER PENALTY FORCE ---
    const double lambda_t = tr.npen > 1 ? tr.pen[t - 1] : tr.pen[0];
    apply_barrier_penalty(O.mux[idx - 1], O.muy[idx - 1], E.barrier_spline, lambda_t, h0, h1);

    if (E.model == 0) {  // overdamped Langevin
      const double m0 = O.mux[idx - 1] + s2 * dt_step / 2.0 * h0;
      const double m1 = O.muy[idx - 1] + s2 * dt_step / 2.0 * h1;
      const double sd = sigma * sqrt(dt_step);
      const double nx = m0 + sd * z[iz++];
      const double ny = m1 + sd * z[iz++];
      if (out_of_bounds(nx, ny, E)) return bounds_message(nx, ny, i + 1, O.time[idx], E);
      O.mux[idx] = nx;
      O.muy[idx] = ny;

    } else {  // underdamped Langevin
      const double exp_gdt = exp(-gamma * dt_step);
      const double exp_2gdt = exp(-2 * gamma * dt_step);

      arma::vec mean(4);
      mean(0) = O.mux[idx - 1] + O.velx[idx - 1] / gamma * (1 - exp_gdt) +
        s2 * h0 / gamma * (dt_step - (1 - exp_gdt) / gamma);
      mean(1) = O.velx[idx - 1] * exp_gdt + s2 * h0 / gamma * (1 - exp_gdt);
      mean(2) = O.muy[idx - 1] + O.vely[idx - 1] / gamma * (1 - exp_gdt) +
        s2 * h1 / gamma * (dt_step - (1 - exp_gdt) / gamma);
      mean(3) = O.vely[idx - 1] * exp_gdt + s2 * h1 / gamma * (1 - exp_gdt);

      const double var_x = s2 / (gamma * gamma) * (2 * gamma * dt_step - 3 + 4 * exp_gdt - exp_2gdt);
      const double var_v = s2 * (1 - exp_2gdt);
      const double cov_xv = s2 / gamma * (1 - 2 * exp_gdt + exp_2gdt);

      arma::mat Sigma(4, 4, arma::fill::zeros);
      Sigma(0,0) = Sigma(2,2) = var_x;
      Sigma(1,1) = Sigma(3,3) = var_v;
      Sigma(0,1) = Sigma(1,0) = cov_xv;
      Sigma(2,3) = Sigma(3,2) = cov_xv;

      arma::vec zz(4);
      for (int k = 0; k < 4; k++) zz(k) = z[iz++];
      arma::vec new_state = rmvnorm_z(mean, Sigma, zz);

      if (out_of_bounds(new_state(0), new_state(2), E)) return bounds_message(new_state(0), new_state(2), i + 1, O.time[idx], E);
      O.mux[idx] = new_state(0);
      O.velx[idx] = new_state(1);
      O.muy[idx] = new_state(2);
      O.vely[idx] = new_state(3);
    }
  }
  return std::string();
}

int sim_threads(int ncores) {
#ifdef _OPENMP
  return ncores > 1 ? ncores : 1;
#else
  (void) ncores;
  return 1;
#endif
}

void sim_all_tracks(const std::vector<SimTrack>& tracks, const SimEnv& E, int ncores, SimOut& O) {
  const int n = (int) tracks.size();
  const int nthr = sim_threads(ncores);
  const int k = draws_per_step(E.model);
  for (int b0 = 0; b0 < n; b0 += nthr) {
    const int b1 = std::min(n, b0 + nthr);
    // Serial draws from R's RNG, in animal order
    std::vector<size_t> zoff(b1 - b0 + 1, 0);
    for (int i = b0; i < b1; i++) zoff[i - b0 + 1] = zoff[i - b0] + 2 + (size_t) tracks[i].nsteps * k;
    std::vector<double> z(zoff.back());
    for (size_t j = 0; j < z.size(); j++) z[j] = norm_rand();
    std::vector<std::string> err(b1 - b0);
#ifdef _OPENMP
    const int nt = std::min(nthr, b1 - b0);
#pragma omp parallel for num_threads(nt) schedule(static)
#endif
    for (int i = b0; i < b1; i++) {
      try {
        err[i - b0] = sim_one_track(i, tracks[i], z.data() + zoff[i - b0], E, O);
      } catch (std::exception& e) {
        err[i - b0] = std::string("simulation failed for animal ") + std::to_string(i + 1) + ": " + e.what();
      } catch (...) {
        err[i - b0] = std::string("simulation failed for animal ") + std::to_string(i + 1);
      }
    }
    for (int i = b0; i < b1; i++) if (!err[i - b0].empty()) throw Rcpp::exception(err[i - b0].c_str());
  }
}

// Common set-up of simulate_langevin_cpp / simulate_langevin_tracks_cpp (serial; R API allowed).
void sim_setup(SimEnv& E, int model, double gamma, double sigma, const NumericVector& beta, const List& raster_data,
               const NumericMatrix& barrier_dist, bool any_barrier, int spline_degree, double barrier_onset) {
  NumericVector raster_extent = raster_data["raster_extent"];
  E.model = model;
  E.gamma = gamma;
  E.sigma = sigma;
  E.beta.assign(beta.begin(), beta.end());
  E.xmin = raster_extent[0]; E.xmax = raster_extent[1];
  E.ymin = raster_extent[2]; E.ymax = raster_extent[3];
  E.n_covs = raster_data["n_covs"];
  if ((int) E.beta.size() < E.n_covs) stop("beta must have one value per covariate.");
  NumericVector raster_vals = raster_data["raster_vals"];
  NumericVector raster_resolution = raster_data["raster_resolution"];
  arma::mat barrier_dist_arma = Rcpp::as<arma::mat>(barrier_dist);
  // Spline surfaces (same representation and degree as the TMB likelihood)
  if (spline_degree != 1 && spline_degree != 3 && spline_degree != 5) stop("spline_degree must be 1, 3 or 5.");
  E.cov_spline = build_cov_spline(raster_vals, raster_resolution, raster_extent, spline_degree);
  if (any_barrier) E.barrier_spline = build_barrier_spline(barrier_dist_arma, raster_resolution, raster_extent, spline_degree, barrier_onset);
}

DataFrame sim_output(int model, SimOut& O) {
  NumericVector ID(O.id.begin(), O.id.end()), time(O.time.begin(), O.time.end()), dt(O.dt.begin(), O.dt.end()),
    mu_x(O.mux.begin(), O.mux.end()), mu_y(O.muy.begin(), O.muy.end());
  if (model == 0) {
    return DataFrame::create(Named("id") = ID, Named("date") = time, Named("dt") = dt,
                             Named("mu.x") = mu_x, Named("mu.y") = mu_y);
  }
  NumericVector velx(O.velx.begin(), O.velx.end()), vely(O.vely.begin(), O.vely.end());
  return DataFrame::create(Named("id") = ID, Named("date") = time, Named("dt") = dt,
                           Named("mu.x") = mu_x, Named("mu.y") = mu_y,
                           Named("vel.x") = velx, Named("vel.y") = vely);
}

} // namespace

// [[Rcpp::export]]
DataFrame simulate_langevin_cpp(int model,
                                int nbAnimals,
                                int obsPerAnimal,
                                NumericVector dt_vec,
                                double gamma,
                                double sigma,
                                NumericVector beta,
                                List raster_data,
                                NumericMatrix initialPosition,
                                NumericMatrix barrier_dist,
                                NumericVector barrier_penalty,
                                int spline_degree = 5,
                                double barrier_onset = 0.0,
                                int ncores = 1) {
  // barrier_penalty: length 1 (same lambda on every step) or length(dt_vec) (lambda for each step,
  // e.g. the stability limit at that step's dt for lambda = "tied")
  // ncores: number of OpenMP threads (animals are simulated in parallel); 1 if compiled without OpenMP.
  if (obsPerAnimal != dt_vec.size()) stop("obsPerAnimal must equal length(dt_vec).");
  if (initialPosition.nrow() < nbAnimals) stop("initialPosition must have nbAnimals rows.");
  if (barrier_penalty.size() < 1) stop("barrier_penalty must have length 1 or length(dt_vec).");
  bool any_barrier = false;
  for (int k = 0; k < barrier_penalty.size(); k++) if (barrier_penalty[k] > 0.0) any_barrier = true;
  const bool per_step_penalty = barrier_penalty.size() > 1;
  if (per_step_penalty && barrier_penalty.size() != dt_vec.size()) stop("barrier_penalty must have length 1 or length(dt_vec).");

  SimEnv E;
  NumericVector raster_extent = raster_data["raster_extent"];
  E.xmin = raster_extent[0]; E.xmax = raster_extent[1]; E.ymin = raster_extent[2]; E.ymax = raster_extent[3];
  // Check initial positions
  for (int i = 0; i < nbAnimals; i++)
    if (out_of_bounds(initialPosition(i,0), initialPosition(i,1), E))
      throw Rcpp::exception(bounds_message(initialPosition(i,0), initialPosition(i,1), i + 1, 0.0, E).c_str());

  sim_setup(E, model, gamma, sigma, beta, raster_data, barrier_dist, any_barrier, spline_degree, barrier_onset);
  NumericVector all_z_values = raster_data["all_z_values"];
  IntegerVector n_zvals_cov = raster_data["n_zvals_cov"];
  IntegerVector cov_offset = raster_data["cov_offset"];
  E.all_z_values = &all_z_values; E.n_zvals_cov = &n_zvals_cov; E.cov_offset = &cov_offset;

  const int n_points = obsPerAnimal + 1;
  std::vector<SimTrack> tracks(nbAnimals);
  for (int i = 0; i < nbAnimals; i++) {
    SimTrack& tr = tracks[i];
    tr.dt = dt_vec.begin(); tr.nsteps = obsPerAnimal;
    tr.pen = barrier_penalty.begin(); tr.npen = barrier_penalty.size();
    tr.x0 = initialPosition(i,0); tr.y0 = initialPosition(i,1);
    tr.start = (size_t) i * n_points;
  }
  SimOut O((size_t) nbAnimals * n_points);
  sim_all_tracks(tracks, E, ncores, O);
  return sim_output(model, O);
}

// Several tracks with their own time steps and barrier penalties (used by simLangevin.fitLangevin).
// Equivalent to calling simulate_langevin_cpp(nbAnimals = 1) for each track in turn (same draws in the
// same order), but the tracks can be simulated in parallel.  Track i has id i + 1.
// [[Rcpp::export]]
DataFrame simulate_langevin_tracks_cpp(int model,
                                       List dt_list,
                                       List penalty_list,
                                       double gamma,
                                       double sigma,
                                       NumericVector beta,
                                       List raster_data,
                                       NumericMatrix initialPosition,
                                       NumericMatrix barrier_dist,
                                       int spline_degree = 5,
                                       double barrier_onset = 0.0,
                                       int ncores = 1) {
  const int n = dt_list.size();
  if (penalty_list.size() != n) stop("dt_list and penalty_list must have the same length.");
  if (initialPosition.nrow() < n) stop("initialPosition must have one row per track.");
  std::vector<NumericVector> dts(n), pens(n);
  bool any_barrier = false;
  for (int i = 0; i < n; i++) {
    dts[i] = as<NumericVector>(dt_list[i]);
    pens[i] = as<NumericVector>(penalty_list[i]);
    if (pens[i].size() < 1 || (pens[i].size() > 1 && pens[i].size() != dts[i].size()))
      stop("each barrier penalty must have length 1 or the length of its time steps.");
    for (int k = 0; k < pens[i].size(); k++) if (pens[i][k] > 0.0) any_barrier = true;
  }
  SimEnv E;
  NumericVector raster_extent = raster_data["raster_extent"];
  E.xmin = raster_extent[0]; E.xmax = raster_extent[1]; E.ymin = raster_extent[2]; E.ymax = raster_extent[3];
  for (int i = 0; i < n; i++)
    if (out_of_bounds(initialPosition(i,0), initialPosition(i,1), E))
      throw Rcpp::exception(bounds_message(initialPosition(i,0), initialPosition(i,1), i + 1, 0.0, E).c_str());

  sim_setup(E, model, gamma, sigma, beta, raster_data, barrier_dist, any_barrier, spline_degree, barrier_onset);
  NumericVector all_z_values = raster_data["all_z_values"];
  IntegerVector n_zvals_cov = raster_data["n_zvals_cov"];
  IntegerVector cov_offset = raster_data["cov_offset"];
  E.all_z_values = &all_z_values; E.n_zvals_cov = &n_zvals_cov; E.cov_offset = &cov_offset;

  std::vector<SimTrack> tracks(n);
  size_t start = 0;
  for (int i = 0; i < n; i++) {
    SimTrack& tr = tracks[i];
    tr.dt = dts[i].begin(); tr.nsteps = dts[i].size();
    tr.pen = pens[i].begin(); tr.npen = pens[i].size();
    tr.x0 = initialPosition(i,0); tr.y0 = initialPosition(i,1);
    tr.start = start;
    start += (size_t) tr.nsteps + 1;
  }
  SimOut O(start);
  sim_all_tracks(tracks, E, ncores, O);
  return sim_output(model, O);
}

// Number of processors available to OpenMP (0 if the package was compiled without OpenMP).
// [[Rcpp::export]]
int lssm_omp_num_procs() {
#ifdef _OPENMP
  return omp_get_num_procs();
#else
  return 0;
#endif
}

// [[Rcpp::export]]
List measurementError_rcpp(DataFrame data,
                           double smaj_sd,
                           double smin_sd,
                           NumericVector eor,
                           double psi,
                           IntegerVector obs_dist,
                           NumericVector df,
                           bool knownError = false,
                           bool ell_power = false,
                           double psi_maj = 1.0,
                           double b_maj = 1.0,
                           double b_min = 1.0,
                           double ref_maj = 0.0,
                           double ref_min = 0.0) {
  // obs_dist: distribution of each observation (0 Gaussian, 1 t, 2 Huber; length n, or 1 = shared)
  // df: shape parameter (df or k) of each observation (length n, or 1 = shared)
  // ell_power: effective axes psi_maj * ref_maj * (smaj / ref_maj)^b_maj and psi * ref_min * (smin / ref_min)^b_min;
  // reference lengths <= 0 are replaced by the geometric means of the axes

  int n = data.nrows();
  NumericVector M_rand(n);
  NumericVector m_rand(n);
  NumericVector c_rand(n);

  NumericVector mux = data["mu.x"];
  NumericVector muy = data["mu.y"];

  // Initialize new positions. If knownError=true and an earlier step already
  // generated x/y, preserve them. Otherwise start from true locations.
  NumericVector new_mux;
  NumericVector new_muy;
  if (data.containsElementNamed("x")) {
    new_mux = clone(as<NumericVector>(data["x"]));
    new_muy = clone(as<NumericVector>(data["y"]));
  } else {
    new_mux = clone(mux);
    new_muy = clone(muy);
  }

  if (knownError) {
    NumericVector orig_smaj = data["smaj"];
    NumericVector orig_smin = data["smin"];
    NumericVector orig_eor = data["eor"];

    for(int i = 0; i < n; i++) {
      M_rand[i] = orig_smaj[i];
      m_rand[i] = orig_smin[i];
      if (!NumericVector::is_na(orig_eor[i])) {
        c_rand[i] = orig_eor[i];
      } else {
        c_rand[i] = NA_REAL;
      }
    }
  } else {
    for(int i = 0; i < n; i++) {
      M_rand[i] = std::abs(R::rnorm(0.0, smaj_sd));
      m_rand[i] = std::abs(R::rnorm(0.0, smin_sd));
      if(M_rand[i] < m_rand[i]){
        double tmp = M_rand[i];
        M_rand[i] = m_rand[i];
        m_rand[i] = tmp;
      }
      c_rand[i] = R::runif(eor[0], eor[1]) * M_PI / 180.0;
    }
  }

  double z = std::sqrt(2.0);

  if (ell_power && (ref_maj <= 0.0 || ref_min <= 0.0)) {
    double lM = 0.0, lm = 0.0;
    int nk = 0;
    for(int i = 0; i < n; i++) {
      if (NumericVector::is_na(M_rand[i]) || NumericVector::is_na(m_rand[i])) continue;
      lM += std::log(M_rand[i]); lm += std::log(m_rand[i]); nk++;
    }
    if (nk > 0) {
      if (ref_maj <= 0.0) ref_maj = std::exp(lM / nk);
      if (ref_min <= 0.0) ref_min = std::exp(lm / nk);
    }
  }

  // Apply measurement error
  for(int i = 0; i < n; i++) {
    if (NumericVector::is_na(M_rand[i])) continue; // Skip if no KF error for this row

    double Ma = M_rand[i];
    double ma = m_rand[i] * psi;
    if (ell_power) {
      Ma = psi_maj * ref_maj * std::pow(M_rand[i] / ref_maj, b_maj);
      ma = psi * ref_min * std::pow(m_rand[i] / ref_min, b_min);
    }
    double s2c = std::sin(c_rand[i]) * std::sin(c_rand[i]);
    double c2c = std::cos(c_rand[i]) * std::cos(c_rand[i]);
    double M2 = (Ma / z) * (Ma / z);
    double m2 = (ma / z) * (ma / z);

    arma::mat cov_obs(2, 2, arma::fill::zeros);
    cov_obs(0,0) = M2 * s2c + m2 * c2c;
    cov_obs(1,1) = M2 * c2c + m2 * s2c;
    cov_obs(0,1) = (0.5 * (Ma * Ma - ma * ma)) * std::cos(c_rand[i]) * std::sin(c_rand[i]);
    cov_obs(1,0) = cov_obs(0,1);

    double df_i = df.size() == 1 ? df[0] : df[i];
    arma::vec mean = {mux[i], muy[i]};
    arma::vec new_pos;
    int dist_i = obs_dist.size() == 1 ? obs_dist[0] : obs_dist[i];
    if(dist_i == 1) new_pos = rmvt(mean, cov_obs, df_i);
    else if(dist_i == 2) new_pos = rmvhyp(mean, cov_obs, df_i);
    else new_pos = rmvnorm(mean, cov_obs);

    new_mux[i] = new_pos(0);
    new_muy[i] = new_pos(1);
  }

  if (knownError) {
    return List::create(Named("x") = new_mux, Named("y") = new_muy);
  } else {
    return List::create(Named("x") = new_mux, Named("y") = new_muy,
                        Named("smaj") = M_rand, Named("smin") = m_rand, Named("eor") = c_rand);
  }
}

// [[Rcpp::export]]
List measurementError_LS_rcpp(DataFrame data,
                              double x_sd,
                              double y_sd,
                              NumericVector tau_x,
                              NumericVector tau_y,
                              double rho_o,
                              IntegerVector obs_dist,
                              NumericVector df,
                              bool knownError = false) {
  // tau_x, tau_y, obs_dist, df: values for each observation (length n, or 1 = shared)

  int n = data.nrows();
  NumericVector mux = data["mu.x"];
  NumericVector muy = data["mu.y"];

  NumericVector obs_x;
  NumericVector obs_y;
  if (data.containsElementNamed("x")) {
    obs_x = clone(as<NumericVector>(data["x"]));
    obs_y = clone(as<NumericVector>(data["y"]));
  } else {
    obs_x = clone(mux);
    obs_y = clone(muy);
  }

  NumericVector x_err_vec(n);
  NumericVector y_err_vec(n);

  if (knownError) {
    x_err_vec = data["x.err"];
    y_err_vec = data["y.err"];
  } else {
    for(int i = 0; i < n; i++) {
      x_err_vec[i] = std::abs(R::rnorm(0.0, x_sd));
      y_err_vec[i] = std::abs(R::rnorm(0.0, y_sd));
    }
  }

  for(int i = 0; i < n; i++) {
    if (NumericVector::is_na(x_err_vec[i])) continue; // Skip if no LS error for this row

    double s = (tau_x.size() == 1 ? tau_x[0] : tau_x[i]) * x_err_vec[i];
    double q = (tau_y.size() == 1 ? tau_y[0] : tau_y[i]) * y_err_vec[i];
    double df_i = df.size() == 1 ? df[0] : df[i];

    arma::mat cov_obs(2, 2, arma::fill::zeros);
    cov_obs(0,0) = s * s;
    cov_obs(1,1) = q * q;
    cov_obs(0,1) = s * q * rho_o;
    cov_obs(1,0) = cov_obs(0,1);

    arma::vec mean = {mux[i], muy[i]};
    arma::vec new_pos;
    int dist_i = obs_dist.size() == 1 ? obs_dist[0] : obs_dist[i];
    if(dist_i == 1) new_pos = rmvt(mean, cov_obs, df_i);
    else if(dist_i == 2) new_pos = rmvhyp(mean, cov_obs, df_i);
    else new_pos = rmvnorm(mean, cov_obs);

    obs_x[i] = new_pos(0);
    obs_y[i] = new_pos(1);
  }

  if (knownError) {
    return List::create(Named("x") = obs_x, Named("y") = obs_y);
  } else {
    return List::create(Named("x") = obs_x, Named("y") = obs_y,
                        Named("x.err") = x_err_vec, Named("y.err") = y_err_vec);
  }
}

// [[Rcpp::export]]
NumericVector extract_sdf_rcpp(NumericMatrix pts,
                               NumericMatrix barrier_dist,
                               NumericVector raster_extent,
                               NumericVector raster_resolution) {
  int n = pts.nrow();
  NumericVector out(n);
  arma::mat grid = Rcpp::as<arma::mat>(barrier_dist);
  for(int i=0; i<n; i++) {
    out[i] = get_bilinear_val(pts(i,0), pts(i,1), grid, raster_extent, raster_resolution);
  }
  return out;
}
