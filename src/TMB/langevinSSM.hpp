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

template<class Type>
Type langevinSSM(objective_function<Type>* obj)
{
  DATA_INTEGER(process_model);      // 0 = overdamped, 1 = underdamped
  DATA_INTEGER(obs_dist);           // 0 = Gaussian, 1 = t-distribution
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
  // for LS/GPS OBS MODEL
  PARAMETER_VECTOR(l_tau);     	// error dispersion for LS obs model (log scale)
  PARAMETER(l_rho_o);           // measurement error correlation b/w x,y
  PARAMETER(l_df);       // Log degrees of freedom for t-distribution

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

  DATA_INTEGER(has_prior_l_tau);
  DATA_VECTOR(prior_mean_l_tau);
  DATA_VECTOR(prior_sd_l_tau);

  DATA_INTEGER(has_prior_l_rho_o);
  DATA_VECTOR(prior_mean_l_rho_o);
  DATA_VECTOR(prior_sd_l_rho_o);

  DATA_INTEGER(has_prior_l_df);
  DATA_VECTOR(prior_mean_l_df);
  DATA_VECTOR(prior_sd_l_df);

  DATA_INTEGER(has_prior_mu);
  DATA_IVECTOR(prior_idx_mu);
  DATA_VECTOR(prior_mean_mu_val);
  DATA_VECTOR(prior_sd_mu_val);

  DATA_INTEGER(has_prior_vel);
  DATA_IVECTOR(prior_idx_vel);
  DATA_VECTOR(prior_mean_vel_val);
  DATA_VECTOR(prior_sd_vel_val);

  Type psi = exp(l_psi);
  vector<Type> tau = exp(l_tau);
  Type rho_o = Type(2.0) / (Type(1.0) + exp(-l_rho_o)) - Type(1.0);
  Type df = exp(l_df);

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

      // Calculate gradient of log(π)
      matrix<Type> grad;
      if(smoothGradient) {
        grad = calculate_smooth_gradient(
          x_prev, y_prev, z_prev,
          cov_spline, all_z_values, n_zvals_cov, cov_offset,
          n_covs, weights, smooth_sigma, smooth_gamma, dt_step, zetaScale, process_model
        );
      } else {
        grad = extract_raster_values(
          x_prev, y_prev, z_prev,
          cov_spline, all_z_values, n_zvals_cov, cov_offset, n_covs
        );
      }

      // Force vector h = ∇log[π(x)]
      vector<Type> h(2);
      h.setZero();
      for(int c = 0; c < n_covs; c++) {
        h(0) += beta(c) * grad(0,c);
        h(1) += beta(c) * grad(1,c);
      }

      // --- BARRIER PENALTY FORCE ---
      Type h0 = h(0), h1 = h(1);
      Type lambda_t = lambda_w;
      if (lambda_mode == 1) lambda_t = lambda_stability_limit(sigma_sca, gamma, dt_step, process_model);
      apply_barrier_penalty(x_prev, y_prev, barrier_spline, lambda_t, h0, h1);
      h(0) = h0; h(1) = h1;
      // -----------------------------

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

        for(int i = 0; i < 2; i++) {
          // Position mean
          Type mu_x_pred = mu(i,idx) +
            vel(i,idx) * term1 +
            s2*h(i)/gamma * (dt_step - term1);

          // Velocity mean
          Type mu_v_pred = vel(i,idx) * exp_gdt +
            s2*h(i) * term1;

          // Construct variance-covariance matrix
          matrix<Type> Sigma(2,2);
          Sigma(0,0) = var_x;
          Sigma(1,1) = var_v;
          Sigma(0,1) = cov_xv;
          Sigma(1,0) = cov_xv;

          vector<Type> nu(2);
          nu(0) = mu(i,idx+1) - mu_x_pred;
          nu(1) = vel(i,idx+1) - mu_v_pred;

          nll += MVNORM(Sigma)(nu);
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
    if(!R_IsNA(asDouble(prior_mean_l_df(0)))) nll -= dnorm(l_df, prior_mean_l_df(0), prior_sd_l_df(0), true);
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
        // Argos Least Squares and GPS observations
        Type s = tau(0) * K(i,0);
        Type q = tau(1) * K(i,1);
        cov_obs(0,0) = s * s;
        cov_obs(1,1) = q * q;
        cov_obs(0,1) = s * q * rho_o;
        cov_obs(1,0) = cov_obs(0,1);
      } else if(obs_mod(i) == 1) {
        // Argos Kalman Filter observations
        Type z = sqrt(Type(2.0));
        Type s2c = sin(eor(i)) * sin(eor(i));
        Type c2c = cos(eor(i)) * cos(eor(i));
        Type M2  = (smaj(i) / z) * (smaj(i) / z);
        Type m2 = (smin(i) * psi / z) * (smin(i) * psi / z);
        cov_obs(0,0) = (M2 * s2c + m2 * c2c);
        cov_obs(1,1) = (M2 * c2c + m2 * s2c);
        cov_obs(0,1) = (0.5 * (smaj(i) * smaj(i) - (smin(i) * psi * smin(i) * psi))) * cos(eor(i)) * sin(eor(i));
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

      if (obs_dist == 0) {
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
      } else if (obs_dist == 1) {
        // Bivariate t-distribution model, written as marginal(x) * conditional(y | x)
        // so that it is a smooth function of the OSA indicators kX and kY.
        // (Do NOT branch on kX/kY with 'if': the branch is fixed when the AD tape is
        // recorded, with keep = 1, so oneStepPredict could never switch data terms off.)
        //   x            ~ t_df(mu_x, varX)
        //   y | x        ~ t_{df+1}(mu_y + covXY/varX * dx,
        //                           (df + dx^2/varX)/(df + 1) * (varY - covXY^2/varX))
        // With kX = kY = 1 this equals the joint bivariate t density exactly.
        Type dx = x_obs - mu_x;
        Type ll_x = dt(dx / sdX, df, true) - log(sdX);

        Type var_y_given_x = varY - covXY * covXY / varX;
        Type mu_y_cond = mu_y + kX * (covXY / varX) * dx;
        Type scale2_y = kX * ((df + dx * dx / varX) / (df + Type(1.0))) * var_y_given_x
        + (Type(1.0) - kX) * varY;
        Type df_y = df + kX;
        Type ll_y = dt((y_obs - mu_y_cond) / sqrt(scale2_y), df_y, true) - Type(0.5) * log(scale2_y);

        nll -= kX * ll_x + kY * ll_y;
        if (osa_cdf) {
          Type lF, lS;
          log_pt_both(dx / sdX, df, lF, lS);
          nll -= keep.cdf_lower(0,i) * lF + keep.cdf_upper(0,i) * lS;
          log_pt_both((y_obs - mu_y_cond) / sqrt(scale2_y), df_y, lF, lS);
          nll -= keep.cdf_lower(1,i) * lF + keep.cdf_upper(1,i) * lS;
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
  if(obs_dist == 1) ADREPORT(df);
  if(lambda_mode > 0) ADREPORT(lambda);

  // Backup in case sdreport fails
  REPORT(beta);
  REPORT(sigma);
  if(process_model == 1) REPORT(gamma);
  REPORT(rho_o);
  REPORT(tau);
  REPORT(psi);
  if(obs_dist == 1) REPORT(df);
  REPORT(lambda);
  REPORT(mu);
  if(process_model == 1) REPORT(vel);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
