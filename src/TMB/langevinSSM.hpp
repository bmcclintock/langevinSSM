// src/TMB/langevinSSM.hpp

#ifndef langevinSSM_hpp
#define langevinSSM_hpp

#include "include/raster_helpers.hpp"

using namespace density;


// Largest barrier strength for which the frozen-force one-step map of a land excursion is stable
// (spectral radius <= 1).  Underdamped, with x = gamma*dt, e = exp(-x), a = 1 - e, in units of
// kappa = sigma^2 lambda / gamma^2:
//   kappa_1 = a / (a - x e)                         (det = 1; binding for x <= 3.72)
//   kappa_2 = 2 (1 + e) / ((1 + e)(x - a) - a^2)     (eigenvalue -1; binding for x > 3.72)
// Overdamped: lambda = 2 / (sigma^2 dt) (one-step removal; unstable above 2x this).
template<class Type>
Type lambda_stability_limit(Type sigma, Type gamma, Type dt, int process_model) {
  if (process_model != 1) return Type(2.0) / (sigma * sigma * dt);
  Type x = gamma * dt;
  Type e = exp(-x);
  Type a = -expm1(-x);
  Type k1 = a / (a - x * e);
  Type d2 = (Type(1.0) + e) * (x - a) - a * a;
  Type k2 = Type(2.0) * (Type(1.0) + e) / CppAD::CondExpGt(d2, Type(1e-12), d2, Type(1e-12));
  Type k = CppAD::CondExpLt(k1, k2, k1, k2);
  return k * gamma * gamma / (sigma * sigma);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// ---------------------------------------------------------------------------
// CDF helpers for TMB oneStepPredict(method = "cdf").
// Both return log CDF values; the tiny offset prevents log(0) in extreme tails.
// ---------------------------------------------------------------------------
// log P(T <= z) and log P(T > z) for a standard Student t with nu df, via the
// regularised incomplete beta function (AD-capable atomic pbeta):
//   P(T > |z|) = 0.5 * I_{nu/(nu+z^2)}(nu/2, 1/2)
template<class Type>
void log_pt_both(Type z, Type nu, Type &logF, Type &logS)
{
  Type z2   = z * z;
  Type tail = Type(0.5) * pbeta(nu / (nu + z2), Type(0.5) * nu, Type(0.5));            // P(T > |z|)
  Type body = Type(0.5) + Type(0.5) * pbeta(z2 / (nu + z2), Type(0.5), Type(0.5) * nu); // P(T < |z|)
  Type F = CppAD::CondExpLt(z, Type(0.0), tail, body);
  Type S = CppAD::CondExpLt(z, Type(0.0), body, tail);
  logF = log(F + Type(1e-300));
  logS = log(S + Type(1e-300));
}
// Same for a standard normal (using symmetry to keep upper-tail precision).
template<class Type>
void log_pnorm_both(Type z, Type &logF, Type &logS)
{
  logF = log(pnorm(z) + Type(1e-300));
  logS = log(pnorm(-z) + Type(1e-300));
}

// log K_nu(z) (modified Bessel function of the second kind) for z > 0, safe for large z.
// For z > 50 the asymptotic expansion is used (relative error < 1e-12 for the orders used here),
// so the result stays finite where besselK itself underflows. Both branches are finite on the tape.
template<class Type>
Type log_besselK(Type z, Type nu)
{
  Type zc = CppAD::CondExpGt(z, Type(50.0), Type(50.0), z);
  Type direct = log(besselK(zc, nu));
  Type mu4 = Type(4.0) * nu * nu;
  Type term = Type(1.0), ser = Type(1.0);
  for (int j = 1; j <= 8; j++) {
    Type a = Type(2.0 * j - 1.0);
    term = term * (mu4 - a * a) / (Type(8.0 * j) * z);
    ser += term;
  }
  Type asym = Type(0.5) * log(Type(M_PI) / (Type(2.0) * z)) - z + log(ser);
  return CppAD::CondExpGt(z, Type(50.0), asym, direct);
}

// Log density of the 1D generalised hyperbolic distribution (symmetric, location 0, unit scale)
//   f(z) = (psi/chi)^(lam/2) / (sqrt(2 pi) K_lam(sqrt(chi psi))) *
//          ((chi + z^2)/psi)^((lam - 1/2)/2) * K_{lam-1/2}(sqrt(psi (chi + z^2)))
// Used for the x-marginal and the y|x-conditional of the bivariate hyperbolic (pseudo-Huber) model.
template<class Type>
Type log_dgh1(Type z, Type lam, Type chi, Type psi)
{
  return Type(0.5) * lam * log(psi / chi) - Type(0.5) * log(Type(2.0 * M_PI)) - log_besselK(sqrt(chi * psi), lam) +
    Type(0.5) * (lam - Type(0.5)) * log((chi + z * z) / psi) + log_besselK(sqrt(psi * (chi + z * z)), lam - Type(0.5));
}

// log Phi(-x) (standard normal upper tail), accurate for large x where pnorm(-x) underflows.
// For x > 8 the asymptotic series is used (relative truncation error < 5e-8); below that pnorm itself is exact.
template<class Type>
Type log_pnorm_upper(Type x)
{
  Type xc = CppAD::CondExpGt(x, Type(8.0), Type(8.0), x);
  Type direct = log(pnorm(-xc));
  Type xa = CppAD::CondExpGt(x, Type(8.0), x, Type(8.0));
  Type r = Type(1.0) / (xa * xa);
  Type ser = Type(1.0) + r * (Type(-1.0) + r * (Type(3.0) + r * (Type(-15.0) + r * (Type(105.0) +
    r * (Type(-945.0) + r * Type(10395.0))))));
  Type asym = -Type(0.5) * xa * xa - Type(0.5) * log(Type(2.0 * M_PI)) - log(xa) + log(ser);
  return CppAD::CondExpGt(x, Type(8.0), asym, direct);
}

// log P(Z > z) for Z ~ GH(lam, chi, psi) (the 1D distribution of log_dgh1).
// Z is a normal variance mixture, Z | W ~ N(0, W) with W ~ GIG(lam, chi, psi), so
//   P(Z > z) = E[ Phi(-z / sqrt(W)) ].
// The expectation is computed by 24-point Gauss-Hermite quadrature in v = log W, centred and scaled at the
// approximate mode of the integrand (for z > 0 the factor Phi(-z/sqrt(w)) ~ exp(-z^2/(2w)) adds z^2 to chi).
// The rule is fixed, so the result is a smooth AD function of (z, lam, chi, psi); the relative error of the
// tail probability is below 2e-5 for k >= 0.3 and below 1e-6 for k >= 2 (checked against numerical integration).
template<class Type>
Type log_pgh1_upper(Type z, Type lam, Type chi, Type psi)
{
  static const double gh_x[24] = {-6.0159255614257443, -5.2593829276680388, -4.6256627564237869, -4.0536644024481499,
                                  -3.5200068130345272, -3.0125461375655633, -2.5238810170114276, -2.0490035736616981, -1.5842500109616937,
                                  -1.1267608176112449, -0.67417110703721195, -0.22441454747251546, 0.2244145474725153, 0.67417110703721228,
                                  1.1267608176112454, 1.5842500109616942, 2.0490035736616989, 2.5238810170114263, 3.0125461375655638,
                                  3.520006813034525, 4.0536644024481516, 4.6256627564237869, 5.259382927668045, 6.0159255614257434};
  // log(weight) + node^2
  static const double gh_lw[24] = {-0.14055535719000289, -0.38776076553999417, -0.51518226609870865, -0.59745971303239642,
                                   -0.65576687610675322, -0.69901095318921058, -0.73173791030169699, -0.75654526290860691, -0.77502142694586862,
                                   -0.7881713515714166, -0.79662923758200099, -0.8007705983701533, -0.80077059837015296, -0.79662923758199922,
                                   -0.7881713515714206, -0.77502142694586817, -0.7565452629085998, -0.73173791030170143, -0.69901095318921058,
                                   -0.65576687610675677, -0.59745971303238932, -0.51518226609870155, -0.38776076553996575, -0.14055535718998868};
                                   // centring: smooth positive part of z
                                   Type zp = Type(0.5) * (z + sqrt(z * z + Type(1.0)));
                                   Type chi2 = chi + zp * zp;
                                   Type l2 = lam - Type(0.5) * zp * zp / (Type(1.0) + zp * zp);
                                   Type wm = (l2 + sqrt(l2 * l2 + chi2 * psi)) / psi;
                                   Type sc = sqrt(Type(2.0)) / sqrt(Type(0.5) * (chi2 / wm + psi * wm));
                                   Type vm = log(wm);
                                   // GIG normalising constant
                                   Type lc = Type(0.5) * lam * log(psi / chi) - log(Type(2.0)) - log_besselK(sqrt(chi * psi), lam);
                                   vector<Type> lv(24);
                                   for (int j = 0; j < 24; j++) {
                                     Type v = vm + sc * Type(gh_x[j]);
                                     Type w = exp(v);
                                     lv(j) = lc + lam * v - Type(0.5) * (chi / w + psi * w) + log_pnorm_upper(z * exp(-Type(0.5) * v)) + Type(gh_lw[j]);
                                   }
                                   Type m = lv(0);
                                   for (int j = 1; j < 24; j++) m = CppAD::CondExpGt(lv(j), m, lv(j), m);
                                   Type sum = Type(0.0);
                                   for (int j = 0; j < 24; j++) sum += exp(lv(j) - m);
                                   return m + log(sum) + log(sc);
}

template<class Type>
Type langevinSSM(objective_function<Type>* obj)
{
  DATA_INTEGER(process_model);      // 0 = overdamped, 1 = underdamped
  DATA_INTEGER(obs_dist);           // any_t + 2 * any_Huber: which tail parameters are used (reporting only)
  DATA_ARRAY(Y);           // Array of observed locations
  DATA_ARRAY_INDICATOR(keep, Y); // used for TMB oneStepPredict bivariate decomposition
  DATA_VECTOR(delta_t);          // Time steps
  DATA_IVECTOR(isd);        // indexes observations vs. interpolation points
  DATA_IVECTOR(obs_mod);    // indicates which obs error model to be used
  DATA_IVECTOR(ID);         // Track IDs
  DATA_IVECTOR(nbObs);      // Number of observations per step
  DATA_SCALAR(scale_factor);
  DATA_IVECTOR(skip_step);  // indicator to skip extremely small (or 0) time steps (assumes no movement between the observations)

  // Raster covariate data
  DATA_ARRAY(raster_vals);         // 3D array of raster values [layer, nrow, ncol]
  DATA_MATRIX(raster_coords);      // Matrix of raster cell coordinates
  DATA_VECTOR(raster_resolution);  // Vector with x,y resolution
  DATA_VECTOR(raster_extent);      // Vector with xmin,xmax,ymin,ymax
  DATA_INTEGER(n_covs);            // Number of unique covariates
  DATA_INTEGER(spline_degree);     // surface interpolation: 1 = bilinear, 3 = bicubic, 5 = quintic B-spline
  DATA_VECTOR(times);              // Absolute time array associated with 'mu' locations
  DATA_VECTOR(all_z_values);           // Vector of times representing the dynamic raster Z-slices
  DATA_IVECTOR(n_zvals_cov);       // Number of layers per covariate (1 for static, >1 for dynamic)
  DATA_IVECTOR(cov_offset);        // Starting layer index (0-based) for each covariate
  DATA_INTEGER(smoothGradient);    // Boolean flag for using smooth gradient calculation
  DATA_VECTOR(weights);            // Vector of 5 (diagonal) or 9 (queen's) weights that sum to 1
  DATA_SCALAR(zetaScale);          // scale factor for smooth gradient neighborhood (>1 increases, <1 decreases)
  DATA_SCALAR(smooth_sigma);       // fixed sigma (working scale) used for the smoothing radius
  DATA_SCALAR(smooth_gamma);       // fixed gamma used for the smoothing radius (underdamped)

  // Barrier constraint data
  DATA_MATRIX(barrier_dist);       // Grid: Signed Distance Function (Positive = allowed, Negative = restricted)
  DATA_SCALAR(barrier_penalty);    // Tuning parameter (lambda) for the severity of the wall
  DATA_INTEGER(lambda_mode);       // 0 = fixed (barrier_penalty); 1 = tied, lambda_t = stability limit at each step's dt; 2 = free, lambda = exp(log_lambda)
  DATA_SCALAR(median_dt);          // median positive step, used to report a representative tied lambda
  DATA_SCALAR(barrier_onset);      // onset length of the barrier potential (working units); 0 = B = 1/2 min(d,0)^2

  // for KF observation model
  DATA_VECTOR(smin);                 //  smin is the semi-minor axis length
  DATA_VECTOR(smaj);                 //  smaj is the semi-major axis length
  DATA_VECTOR(eor);                 //  eor is the orientation of the error ellipse
  // for LS/GPS observation model
  DATA_MATRIX(K);                 // error weighting factors for LS obs model
  // observation groups (0-based): tau group of each LS/GPS observation (column of tau) and
  // tail group of each observation (element of df or k); all 0 when the parameters are shared
  DATA_IVECTOR(tau_group);
  DATA_IVECTOR(dist_obs);   // distribution of each observation: 0 = Gaussian, 1 = t, 2 = pseudo-Huber (bivariate hyperbolic)
  DATA_IVECTOR(df_group);   // df group of each observation (t)
  DATA_IVECTOR(k_group);    // k group of each observation (Huber)
  // reference lengths (working units) of the error ellipse power law
  //   smaj_eff = psi_maj * ref_maj * (smaj / ref_maj)^b_maj,  smin_eff = psi * ref_min * (smin / ref_min)^b_min
  DATA_INTEGER(ell_power);         // 1 = use the power law (psi_maj, b_maj, b_min); 0 = axes as reported (smin scaled by psi)
  DATA_SCALAR(ell_ref_maj);
  DATA_SCALAR(ell_ref_min);

  // Process parameters (simplified)
  PARAMETER_VECTOR(beta);       // Vector of covariate effects
  PARAMETER(log_sigma);  // Log velocity process SD
  PARAMETER(log_gamma);  // Log friction coefficient, only used if model=1
  PARAMETER(log_lambda); // Log barrier strength (working scale), only used if lambda_mode=2

  PARAMETER_MATRIX(mu);        // Locations [2 x timeSteps]
  PARAMETER_MATRIX(vel);        // Velocity innovations [2 x timeSteps], only used if model=1

  // OBSERVATION PARAMETERS
  // for KF OBS MODEL
  PARAMETER(l_psi); 				    // error SD scaling parameter to account for possible uncertainty in Argos error ellipse variables
  PARAMETER(l_psi_maj);         // log scale of the semi-major axis (power law)
  PARAMETER(l_b_maj);           // log exponent of the semi-major axis (power law)
  PARAMETER(l_b_min);           // log exponent of the semi-minor axis (power law)
  // for LS/GPS OBS MODEL
  PARAMETER_VECTOR(l_tau);     	// error dispersion for LS obs model (log scale); (x, y) for each tau group
  PARAMETER(l_rho_o);           // measurement error correlation b/w x,y
  PARAMETER_VECTOR(l_df);       // Log degrees of freedom for t-distribution (one per tail group)
  PARAMETER_VECTOR(l_k);        // Log threshold k of the pseudo-Huber observation model (one per tail group)

  // PRIOR INPUTS
  DATA_INTEGER(has_prior_beta);
  DATA_VECTOR(prior_mean_beta);
  DATA_VECTOR(prior_sd_beta);

  DATA_INTEGER(has_prior_log_sigma);
  DATA_VECTOR(prior_mean_log_sigma);
  DATA_VECTOR(prior_sd_log_sigma);

  DATA_INTEGER(has_prior_log_gamma);
  DATA_VECTOR(prior_mean_log_gamma);
  DATA_VECTOR(prior_sd_log_gamma);

  DATA_INTEGER(has_prior_log_lambda);
  DATA_VECTOR(prior_mean_log_lambda);
  DATA_VECTOR(prior_sd_log_lambda);

  DATA_INTEGER(has_prior_l_psi);
  DATA_VECTOR(prior_mean_l_psi);
  DATA_VECTOR(prior_sd_l_psi);

  DATA_INTEGER(has_prior_l_psi_maj);
  DATA_VECTOR(prior_mean_l_psi_maj);
  DATA_VECTOR(prior_sd_l_psi_maj);

  DATA_INTEGER(has_prior_l_b_maj);
  DATA_VECTOR(prior_mean_l_b_maj);
  DATA_VECTOR(prior_sd_l_b_maj);

  DATA_INTEGER(has_prior_l_b_min);
  DATA_VECTOR(prior_mean_l_b_min);
  DATA_VECTOR(prior_sd_l_b_min);

  DATA_INTEGER(has_prior_l_tau);
  DATA_VECTOR(prior_mean_l_tau);
  DATA_VECTOR(prior_sd_l_tau);

  DATA_INTEGER(has_prior_l_rho_o);
  DATA_VECTOR(prior_mean_l_rho_o);
  DATA_VECTOR(prior_sd_l_rho_o);

  DATA_INTEGER(has_prior_l_df);
  DATA_VECTOR(prior_mean_l_df);
  DATA_VECTOR(prior_sd_l_df);

  DATA_INTEGER(has_prior_l_k);
  DATA_VECTOR(prior_mean_l_k);
  DATA_VECTOR(prior_sd_l_k);

  DATA_INTEGER(has_prior_mu);
  DATA_IVECTOR(prior_idx_mu);
  DATA_VECTOR(prior_mean_mu_val);
  DATA_VECTOR(prior_sd_mu_val);

  DATA_INTEGER(has_prior_vel);
  DATA_IVECTOR(prior_idx_vel);
  DATA_VECTOR(prior_mean_vel_val);
  DATA_VECTOR(prior_sd_vel_val);

  Type psi = exp(l_psi);
  Type psi_maj = exp(l_psi_maj);
  Type b_maj = exp(l_b_maj);
  Type b_min = exp(l_b_min);
  vector<Type> tau = exp(l_tau);
  Type rho_o = Type(2.0) / (Type(1.0) + exp(-l_rho_o)) - Type(1.0);
  vector<Type> df = exp(l_df);
  vector<Type> k = exp(l_k);

  Type nll = 0.0;

  // Create map of unique IDs to track lengths
  int curr_id = ID(0);
  int track_start = 0;
  vector<int> track_starts(1);
  vector<int> track_lengths(1);
  track_starts(0) = 0;

  int timeSteps = delta_t.size();

  // Find track starts and lengths
  for(int i = 1; i < timeSteps; i++) {
    if(ID(i) != curr_id) {
      track_lengths(track_lengths.size()-1) = i - track_start;
      track_starts.conservativeResize(track_starts.size()+1);
      track_lengths.conservativeResize(track_lengths.size()+1);
      track_starts(track_starts.size()-1) = i;
      track_start = i;
      curr_id = ID(i);
    }
  }
  // Add length of last track
  track_lengths(track_lengths.size()-1) = timeSteps - track_start;

  // Transform parameters
  Type sigma = exp(log_sigma + log(scale_factor));   // Original scale
  Type sigma_sca = exp(log_sigma);
  Type gamma = exp(log_gamma);

  // Spline representations of the covariate stack and barrier distance (built from DATA once
  // per template evaluation; evaluation is an atomic operator, so the raster cell is chosen
  // from the current location at every tape evaluation).
  langevin_raster::SplinePtr cov_spline = build_cov_spline(raster_vals, raster_resolution, raster_extent, spline_degree);
  // Barrier strength.  lambda is a parameter of the drift (it enters only through the force
  // -lambda * min(d,0) * grad d), not an additive penalty on the objective.
  // Working scale (coordinates / scale_factor): lambda_w = lambda * scale_factor^2.  log_lambda and
  // barrier_penalty are on the working scale; the reported lambda is on the natural scale.
  Type lambda_w = barrier_penalty;
  if (lambda_mode == 1) {
    // representative (reported) value at the median step; the likelihood uses each step's own limit
    lambda_w = lambda_stability_limit(sigma_sca, gamma, median_dt, process_model);
  } else if (lambda_mode == 2) {
    lambda_w = exp(log_lambda);
  }
  Type lambda = lambda_w / (scale_factor * scale_factor);   // natural scale (reported)

  langevin_raster::SplinePtr barrier_spline;
  if (barrier_penalty > Type(0.0) || lambda_mode > 0) barrier_spline = build_barrier_spline(barrier_dist, raster_resolution, raster_extent, spline_degree, asDouble(barrier_onset));

  std::vector<int> f_lay, f_ci; std::vector<double> f_wl;   // per-step covariate layers / time weights

  // Process model - loop over tracks
  for(int id = 0; id < track_starts.size(); id++) {
    int start_idx = track_starts(id);
    int track_length = track_lengths(id);

    Type s2 = pow(sigma_sca, Type(2.0));

    // Add prior on initial velocities for underdamped model only
    if(process_model == 1) {
      for(int i = 0; i < 2; i++) {
        nll -= dnorm(vel(i,start_idx), Type(0.0), sigma_sca / sqrt(Type(2.0) * gamma), true);
      }
    }

    for(int t = 0; t < (track_length-1); t++) {

      int idx = start_idx + t;
      Type dt_step = delta_t(idx+1);

      if(skip_step(idx+1) == 1) continue;

      // Get current locations and absolute time
      Type x_prev = mu(0,idx);
      Type y_prev = mu(1,idx);
      Type z_prev = times(idx);

      // Force h = grad log pi(x) - lambda_t grad B(x) = grad Phi, with the movement potential
      //   Phi = sum_c beta_c S_c - lambda_t B
      // evaluated by one fused operator per point (PotentialJetOp): its outputs are the partials of
      // Phi itself rather than of every covariate layer, so the force costs the same number of taped
      // operations whatever the number of covariates.
      Type lambda_t = lambda_w;
      if (lambda_mode == 1) lambda_t = lambda_stability_limit(sigma_sca, gamma, dt_step, process_model);
      langevin_raster::SplinePtr bar_t;
      if (barrier_spline && asDouble(lambda_t) > 0.0) bar_t = barrier_spline;
      cov_layer_weights(asDouble(z_prev), all_z_values, n_zvals_cov, cov_offset, n_covs, f_lay, f_wl, f_ci);
      Type h0 = Type(0.0), h1 = Type(0.0);
      if (smoothGradient) {
        // Ring average of the covariate force (see calculate_smooth_gradient); barrier at the centre only.
        Type rad = smooth_gradient_radius(smooth_sigma, smooth_gamma, dt_step, Type(weights(0)), zetaScale, process_model);
        potential_force(x_prev, y_prev, cov_spline, f_lay, f_wl, f_ci, beta, asDouble(weights(0)), bar_t, lambda_t, h0, h1);
        int n_points = weights.size();
        langevin_raster::SplinePtr nobar;
        for (int ip = 0; ip < n_points - 1; ip++) {
          Type angle = (n_points == 9) ? Type(M_PI) * (Type(2.0) * ip + Type(1.0)) / Type(8.0)
            : Type(M_PI) * (Type(2.0) * ip + Type(1.0)) / Type(4.0);
          potential_force(Type(x_prev + rad * cos(angle)), Type(y_prev + rad * sin(angle)), cov_spline, f_lay, f_wl, f_ci,
                          beta, asDouble(weights(ip + 1)), nobar, lambda_t, h0, h1);
        }
      } else {
        potential_force(x_prev, y_prev, cov_spline, f_lay, f_wl, f_ci, beta, 1.0, bar_t, lambda_t, h0, h1);
      }
      vector<Type> h(2);
      h(0) = h0; h(1) = h1;

      if(process_model == 1) {  // Underdamped
        Type gdt = gamma * dt_step;
        Type exp_gdt = exp(-gdt);

        Type em1_gdt  = expm1(-gdt);
        Type em1_2gdt = expm1(-Type(2.0) * gdt);

        // (1 - exp(-gdt))/gamma, reused in both means
        Type term1 = -em1_gdt / gamma;

        // Variance terms
        Type var_x = s2/(gamma*gamma) * (Type(2.0)*gdt + Type(4.0)*em1_gdt - em1_2gdt);
        Type var_v = -s2 * em1_2gdt;
        Type cov_xv = s2/gamma * (em1_gdt * em1_gdt);

        // The two coordinates are independent with the SAME 2x2 covariance, so its inverse and
        // log-determinant are computed once per step (closed form) rather than by two MVNORM calls.
        // Per coordinate this is exactly MVNORM(Sigma)(nu) = 0.5 log|Sigma| + 0.5 nu' Sigma^-1 nu + log(2 pi).
        // Every quantity below except nx and nv depends only on parameters and dt, so the terms that
        // depend on the latent states (which drive the size of the inner Hessian tape) are minimal.
        Type det_S = var_x * var_v - cov_xv * cov_xv;
        Type qa = Type(0.5) * var_v / det_S;    // 0.5 * (Sigma^-1)_11
        Type qb = -cov_xv / det_S;              //       (Sigma^-1)_12   (counted twice in the quadratic form)
        Type qc = Type(0.5) * var_x / det_S;    // 0.5 * (Sigma^-1)_22
        Type cB = s2 / gamma * (dt_step - term1);   // coefficient of h in the position mean
        Type cV = s2 * term1;                       // coefficient of h in the velocity mean
        nll += log(det_S) + Type(2.0 * log(2.0 * M_PI));

        for(int i = 0; i < 2; i++) {
          // Position and velocity means:
          //   mu_x_pred = mu + vel (1 - e^{-gdt})/gamma + s2 h / gamma (dt - (1 - e^{-gdt})/gamma)
          //   mu_v_pred = vel e^{-gdt} + s2 h (1 - e^{-gdt})/gamma
          Type nx = mu(i,idx+1) - (mu(i,idx) + vel(i,idx) * term1 + h(i) * cB);
          Type nv = vel(i,idx+1) - (vel(i,idx) * exp_gdt + h(i) * cV);
          nll += nx * (qa * nx + qb * nv) + qc * nv * nv;
        }
      } else {  // Overdamped

        // Expected locations (simplified case of underdamped as γ → ∞)
        Type Ex = x_prev + s2 * dt_step / Type(2.0) * h(0);
        Type Ey = y_prev + s2 * dt_step / Type(2.0) * h(1);

        Type sd = sigma_sca * sqrt(dt_step);

        nll -= (dnorm(mu(0,idx+1), Ex, sd, true) +
          dnorm(mu(1,idx+1), Ey, sd, true));
      }
    }
  }

  // PRIOR PENALIZATION

  // beta
  if (has_prior_beta == 1) {
    for(int i=0; i<beta.size(); i++){
      if(!R_IsNA(asDouble(prior_mean_beta(i)))) {
        nll -= dnorm(beta(i), prior_mean_beta(i), prior_sd_beta(i), true);
      }
    }
  }

  // log_sigma
  if (has_prior_log_sigma == 1) {
    if(!R_IsNA(asDouble(prior_mean_log_sigma(0)))) {
      Type log_sigma_user = log_sigma + log(scale_factor);
      nll -= dnorm(log_sigma_user, prior_mean_log_sigma(0), prior_sd_log_sigma(0), true);
    }
  }

  // log_gamma
  if (has_prior_log_gamma == 1 && process_model == 1) {
    if(!R_IsNA(asDouble(prior_mean_log_gamma(0)))) {
      nll -= dnorm(log_gamma, prior_mean_log_gamma(0), prior_sd_log_gamma(0), true);
    }
  }

  // l_psi
  // log_lambda (prior on the natural scale; only used when lambda is estimated)
  if (has_prior_log_lambda == 1 && lambda_mode == 2) {
    if(!R_IsNA(asDouble(prior_mean_log_lambda(0)))) {
      Type log_lambda_user = log_lambda - Type(2.0) * log(scale_factor);
      nll -= dnorm(log_lambda_user, prior_mean_log_lambda(0), prior_sd_log_lambda(0), true);
    }
  }

  if (has_prior_l_psi == 1) {
    if(!R_IsNA(asDouble(prior_mean_l_psi(0)))) {
      nll -= dnorm(l_psi, prior_mean_l_psi(0), prior_sd_l_psi(0), true);
    }
  }

  // l_tau
  if (has_prior_l_psi_maj == 1) {
    if(!R_IsNA(asDouble(prior_mean_l_psi_maj(0)))) nll -= dnorm(l_psi_maj, prior_mean_l_psi_maj(0), prior_sd_l_psi_maj(0), true);
  }
  if (has_prior_l_b_maj == 1) {
    if(!R_IsNA(asDouble(prior_mean_l_b_maj(0)))) nll -= dnorm(l_b_maj, prior_mean_l_b_maj(0), prior_sd_l_b_maj(0), true);
  }
  if (has_prior_l_b_min == 1) {
    if(!R_IsNA(asDouble(prior_mean_l_b_min(0)))) nll -= dnorm(l_b_min, prior_mean_l_b_min(0), prior_sd_l_b_min(0), true);
  }

  if (has_prior_l_tau == 1) {
    for(int i=0; i<l_tau.size(); i++){
      if(!R_IsNA(asDouble(prior_mean_l_tau(i)))) {
        nll -= dnorm(l_tau(i), prior_mean_l_tau(i), prior_sd_l_tau(i), true);
      }
    }
  }

  // l_rho_o
  if (has_prior_l_rho_o == 1) {
    if(!R_IsNA(asDouble(prior_mean_l_rho_o(0)))) {
      nll -= dnorm(l_rho_o, prior_mean_l_rho_o(0), prior_sd_l_rho_o(0), true);
    }
  }

  // l_df
  if (has_prior_l_df == 1) {
    for(int i=0; i<l_df.size(); i++){
      if(!R_IsNA(asDouble(prior_mean_l_df(i)))) nll -= dnorm(l_df(i), prior_mean_l_df(i), prior_sd_l_df(i), true);
    }
  }

  // l_k
  if (has_prior_l_k == 1) {
    for(int i=0; i<l_k.size(); i++){
      if(!R_IsNA(asDouble(prior_mean_l_k(i)))) nll -= dnorm(l_k(i), prior_mean_l_k(i), prior_sd_l_k(i), true);
    }
  }

  // mu
  if (has_prior_mu == 1) {
    for(int k=0; k < prior_idx_mu.size(); k++){
      int idx = prior_idx_mu(k);
      int i = idx % 2; // Row: 0 for x, 1 for y
      int j = idx / 2; // Column (time step)
      Type mu_user = mu(i,j) * scale_factor;
      nll -= dnorm(mu_user, prior_mean_mu_val(k), prior_sd_mu_val(k), true);
    }
  }

  // vel
  if (has_prior_vel == 1 && process_model == 1) {
    for(int k=0; k < prior_idx_vel.size(); k++){
      int idx = prior_idx_vel(k);
      int i = idx % 2;
      int j = idx / 2;
      Type vel_user = vel(i,j) * scale_factor;
      nll -= dnorm(vel_user, prior_mean_vel_val(k), prior_sd_vel_val(k), true);
    }
  }

  // OBSERVATION MODEL
  for(int i = 0; i < timeSteps; ++i) {
    if(isd(i) == 1) {
      matrix<Type> cov_obs(2, 2);
      cov_obs.setZero();
      if(obs_mod(i) == 0) {
        // Argos Least Squares and GPS observations (tau of the observation's tau group)
        int g = tau_group(i);
        Type s = tau(2 * g) * K(i,0);
        Type q = tau(2 * g + 1) * K(i,1);
        cov_obs(0,0) = s * s;
        cov_obs(1,1) = q * q;
        cov_obs(0,1) = s * q * rho_o;
        cov_obs(1,0) = cov_obs(0,1);
      } else if(obs_mod(i) == 1) {
        // Argos Kalman Filter observations
        // Effective axes: power law relative to the reference lengths (identity when psi_maj = b_maj = b_min = 1)
        Type Ma = smaj(i);
        Type ma = smin(i) * psi;
        if (ell_power == 1) {
          Ma = psi_maj * ell_ref_maj * exp(b_maj * log(smaj(i) / ell_ref_maj));
          ma = psi * ell_ref_min * exp(b_min * log(smin(i) / ell_ref_min));
        }
        Type z = sqrt(Type(2.0));
        Type s2c = sin(eor(i)) * sin(eor(i));
        Type c2c = cos(eor(i)) * cos(eor(i));
        Type M2  = (Ma / z) * (Ma / z);
        Type m2 = (ma / z) * (ma / z);
        cov_obs(0,0) = (M2 * s2c + m2 * c2c);
        cov_obs(1,1) = (M2 * c2c + m2 * s2c);
        cov_obs(0,1) = (0.5 * (Ma * Ma - ma * ma)) * cos(eor(i)) * sin(eor(i));
        cov_obs(1,0) = cov_obs(0,1);
      }

      // TMB oneStepPredict Bivariate Decomposition
      Type varX = cov_obs(0,0);
      Type varY = cov_obs(1,1);
      Type covXY = cov_obs(0,1);

      Type sdX = sqrt(varX);
      Type sdY = sqrt(varY);

      Type x_obs = Y(0,i);
      Type y_obs = Y(1,i);
      Type mu_x = mu(0,i);
      Type mu_y = mu(1,i);

      Type kX = keep(0,i);
      Type kY = keep(1,i);
      // oneStepPredict(method = "cdf") passes 3*n indicators (keep, lower, upper);
      // the CDF terms are only recorded on the tape in that case.
      // (the lower/upper indicators are AD variables only under method = "cdf").
      bool osa_cdf = keep.osa_active() &&
        (CppAD::Variable(keep.cdf_lower(0,i)) || CppAD::Variable(keep.cdf_lower(1,i)));

      Type df_i = df(df_group(i));   // tail parameters of the observation's groups
      Type k_i = k(k_group(i));
      int dist_i = dist_obs(i);

      if (!keep.osa_active() && dist_i <= 1) {
        // Model fitting (no OSA residuals): all keep indicators are 1, so the marginal x conditional
        // decompositions below reduce to the joint bivariate density, evaluated directly
        // (fewer operations on the latent locations).  With d2 = r' Sigma^-1 r, r = obs - mu:
        //   Gaussian: log f = -log(2 pi) - 1/2 log|Sigma| - d2 / 2
        //   t (df):   log f = -log(2 pi) - 1/2 log|Sigma| - (df + 2)/2 log(1 + d2 / df)
        // (for the t, lgamma(df/2 + 1) - lgamma(df/2) - log(df pi) = -log(2 pi) exactly).
        Type dx = x_obs - mu_x;
        Type ey = (y_obs - mu_y) - (covXY / varX) * dx;
        Type var_y_given_x = varY - covXY * covXY / varX;
        Type d2 = dx * dx / varX + ey * ey / var_y_given_x;
        Type ll_c = -Type(log(2.0 * M_PI)) - Type(0.5) * log(varX * var_y_given_x);
        if (dist_i == 0) nll -= ll_c - Type(0.5) * d2;
        else nll -= ll_c - Type(0.5) * (df_i + Type(2.0)) * log(Type(1.0) + d2 / df_i);
      } else if (dist_i == 0) {
        // Gaussian Model
        nll -= kX * dnorm(x_obs, mu_x, sdX, true);
        Type mu_y_cond = mu_y + kX * (covXY / varX) * (x_obs - mu_x);
        Type sd_y_cond = sqrt(varY - kX * (covXY * covXY / varX));
        nll -= kY * dnorm(y_obs, mu_y_cond, sd_y_cond, true);
        if (osa_cdf) {
          Type lF, lS;
          log_pnorm_both((x_obs - mu_x) / sdX, lF, lS);
          nll -= keep.cdf_lower(0,i) * lF + keep.cdf_upper(0,i) * lS;
          log_pnorm_both((y_obs - mu_y_cond) / sd_y_cond, lF, lS);
          nll -= keep.cdf_lower(1,i) * lF + keep.cdf_upper(1,i) * lS;
        }
      } else if (dist_i == 1) {
        // Bivariate t-distribution model, written as marginal(x) * conditional(y | x)
        // so that it is a smooth function of the OSA indicators kX and kY.
        // (Do NOT branch on kX/kY with 'if': the branch is fixed when the AD tape is
        // recorded, with keep = 1, so oneStepPredict could never switch data terms off.)
        //   x            ~ t_df(mu_x, varX)
        //   y | x        ~ t_{df_i+1}(mu_y + covXY/varX * dx,
        //                           (df_i + dx^2/varX)/(df_i + 1) * (varY - covXY^2/varX))
        // With kX = kY = 1 this equals the joint bivariate t density exactly.
        Type dx = x_obs - mu_x;
        Type ll_x = dt(dx / sdX, df_i, true) - log(sdX);

        Type var_y_given_x = varY - covXY * covXY / varX;
        Type mu_y_cond = mu_y + kX * (covXY / varX) * dx;
        Type scale2_y = kX * ((df_i + dx * dx / varX) / (df_i + Type(1.0))) * var_y_given_x
        + (Type(1.0) - kX) * varY;
        Type df_y = df_i + kX;
        Type ll_y = dt((y_obs - mu_y_cond) / sqrt(scale2_y), df_y, true) - Type(0.5) * log(scale2_y);

        nll -= kX * ll_x + kY * ll_y;
        if (osa_cdf) {
          Type lF, lS;
          log_pt_both(dx / sdX, df_i, lF, lS);
          nll -= keep.cdf_lower(0,i) * lF + keep.cdf_upper(0,i) * lS;
          log_pt_both((y_obs - mu_y_cond) / sqrt(scale2_y), df_y, lF, lS);
          nll -= keep.cdf_lower(1,i) * lF + keep.cdf_upper(1,i) * lS;
        }
      } else if (dist_i == 2) {
        // Pseudo-Huber model (bivariate hyperbolic distribution):
        //   -log f = k sqrt(k^2 + d2) - k^2 + log(2 pi) + 1/2 log|Sigma| + log(1 + 1/k^2),
        // with d2 = r' Sigma^-1 r and r = obs - mu. It is quadratic (Gaussian) for d << k and
        // linear in d for d >> k (exponential tails). sqrt(k^2 + d2) is a norm of (k, Sigma^-1/2 r), so
        // -log f is convex in mu: the observation model can never make the inner (latent state)
        // problem non-convex, unlike the t-distribution. As k -> Inf it tends to the Gaussian.
        // Written as kX * log f(x) + kY * log f(y | x) so that it is smooth in the OSA indicators.
        // Exactly:
        //   x / sdX                       ~ GH(lambda = 3/2, chi = k^2, psi = k^2)
        //   (y - mu_y|x) / sd_y|x | x     ~ GH(lambda = 1,   chi = k^2 + (x/sdX)^2, psi = k^2)
        // log f(y | x) is evaluated as log f(x, y) - log f(x), and log f(x) (Bessel functions) is only
        // recorded on the tape when OSA residuals are being computed.
        Type dx = x_obs - mu_x;
        Type dy = y_obs - mu_y;
        Type var_y_given_x = varY - covXY * covXY / varX;
        Type det = varX * var_y_given_x;
        Type d2 = dx * dx / varX + pow(dy - (covXY / varX) * dx, 2) / var_y_given_x;
        Type k2 = k_i * k_i;
        Type ll_joint = -(k_i * sqrt(k2 + d2) - k2) - log(Type(2.0 * M_PI)) - Type(0.5) * log(det) -
          log(Type(1.0) + Type(1.0) / k2);
        if (!keep.osa_active()) {
          nll -= ll_joint;
        } else {
          Type zx = dx / sdX;
          Type ll_x = log_dgh1(zx, Type(1.5), k2, k2) - log(sdX);
          Type ll_ymarg = log_dgh1(dy / sdY, Type(1.5), k2, k2) - log(sdY);
          Type ll_y = kX * (ll_joint - ll_x) + (Type(1.0) - kX) * ll_ymarg;
          nll -= kX * ll_x + kY * ll_y;
          if (osa_cdf) {
            // x / sdX ~ GH(3/2, k^2, k^2). Given x (kX = 1) the standardised y | x ~ GH(1, k^2 + zx^2, k^2);
            // without x (kX = 0) y / sdY ~ GH(3/2, k^2, k^2). Blended smoothly in kX, exact at kX = 0 and 1.
            Type zy = (dy - kX * (covXY / varX) * dx) / sqrt(kX * var_y_given_x + (Type(1.0) - kX) * varY);
            Type lam_y = Type(1.5) - Type(0.5) * kX;
            Type chi_y = k2 + kX * zx * zx;
            Type lF = log_pgh1_upper(-zx, Type(1.5), k2, k2);
            Type lS = log_pgh1_upper(zx, Type(1.5), k2, k2);
            nll -= keep.cdf_lower(0,i) * lF + keep.cdf_upper(0,i) * lS;
            lF = log_pgh1_upper(-zy, lam_y, chi_y, k2);
            lS = log_pgh1_upper(zy, lam_y, chi_y, k2);
            nll -= keep.cdf_lower(1,i) * lF + keep.cdf_upper(1,i) * lS;
          }
        }
      }
    }
  }

  ADREPORT(beta);
  ADREPORT(sigma);
  if(process_model == 1) ADREPORT(gamma);
  ADREPORT(rho_o);
  ADREPORT(tau);
  ADREPORT(psi);
  if(ell_power == 1) {
    ADREPORT(psi_maj);
    ADREPORT(b_maj);
    ADREPORT(b_min);
  }
  if(obs_dist == 1 || obs_dist == 3) ADREPORT(df);
  if(obs_dist >= 2) ADREPORT(k);
  if(lambda_mode > 0) ADREPORT(lambda);

  // Backup in case sdreport fails
  REPORT(beta);
  REPORT(sigma);
  if(process_model == 1) REPORT(gamma);
  REPORT(rho_o);
  REPORT(tau);
  REPORT(psi);
  if(ell_power == 1) {
    REPORT(psi_maj);
    REPORT(b_maj);
    REPORT(b_min);
  }
  if(obs_dist == 1 || obs_dist == 3) REPORT(df);
  if(obs_dist >= 2) REPORT(k);
  REPORT(lambda);
  REPORT(mu);
  if(process_model == 1) REPORT(vel);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
