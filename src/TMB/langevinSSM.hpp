// src/TMB/langevinSSM.hpp

#ifndef langevinSSM_hpp
#define langevinSSM_hpp

#include "include/raster_helpers.hpp"

using namespace density;

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type langevinSSM(objective_function<Type>* obj)
{
  DATA_INTEGER(process_model);      // 0 = overdamped, 1 = underdamped
  DATA_ARRAY(Y);           // Array of observed locations
  DATA_ARRAY_INDICATOR(keep, Y); // used for TMB oneStepPredict bivariate decomposition
  DATA_VECTOR(dt);          // Time steps
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
  DATA_VECTOR(times);              // Absolute time array associated with 'mu' locations
  DATA_VECTOR(all_z_values);           // Vector of times representing the dynamic raster Z-slices
  DATA_IVECTOR(n_zvals_cov);       // Number of layers per covariate (1 for static, >1 for dynamic)
  DATA_IVECTOR(cov_offset);        // Starting layer index (0-based) for each covariate
  DATA_INTEGER(smoothGradient);    // Boolean flag for using smooth gradient calculation
  DATA_VECTOR(weights);            // Vector of 5 (diagonal) or 9 (queen's) weights that sum to 1
  DATA_SCALAR(zetaScale);          // scale factor for smooth gradient neighborhood (>1 increases, <1 decreases)

  // Barrier constraint data
  DATA_MATRIX(barrier_dist);       // Grid: Signed Distance Function (Positive = allowed, Negative = restricted)
  DATA_SCALAR(barrier_penalty);    // Tuning parameter (lambda) for the severity of the wall

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

  PARAMETER_MATRIX(mu);        // Locations [2 x timeSteps]
  PARAMETER_MATRIX(vel);        // Velocity innovations [2 x timeSteps], only used if model=1

  // OBSERVATION PARAMETERS
  // for KF OBS MODEL
  PARAMETER(l_psi); 				    // error SD scaling parameter to account for possible uncertainty in Argos error ellipse variables
  // for LS/GPS OBS MODEL
  PARAMETER_VECTOR(l_tau);     	// error dispersion for LS obs model (log scale)
  PARAMETER(l_rho_o);           // measurement error correlation b/w x,y

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

  DATA_INTEGER(has_prior_l_psi);
  DATA_VECTOR(prior_mean_l_psi);
  DATA_VECTOR(prior_sd_l_psi);

  DATA_INTEGER(has_prior_l_tau);
  DATA_VECTOR(prior_mean_l_tau);
  DATA_VECTOR(prior_sd_l_tau);

  DATA_INTEGER(has_prior_l_rho_o);
  DATA_VECTOR(prior_mean_l_rho_o);
  DATA_VECTOR(prior_sd_l_rho_o);

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

  Type nll = 0.0;

  // Create map of unique IDs to track lengths
  int curr_id = ID(0);
  int track_start = 0;
  vector<int> track_starts(1);
  vector<int> track_lengths(1);
  track_starts(0) = 0;

  int timeSteps = dt.size();

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
      Type dt_step = dt(idx+1);

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
          raster_vals, all_z_values, n_zvals_cov, cov_offset,
          raster_coords, raster_resolution, raster_extent,
          n_covs, weights, sigma_sca, gamma, dt_step, zetaScale, process_model
        );
      } else {
        grad = extract_raster_values(
          x_prev, y_prev, z_prev,
          raster_vals, all_z_values, n_zvals_cov, cov_offset,
          raster_coords, raster_resolution, raster_extent, n_covs
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
      apply_barrier_penalty(x_prev, y_prev, barrier_dist, raster_extent, raster_resolution, barrier_penalty, h0, h1);
      h(0) = h0; h(1) = h1;
      // -----------------------------

      if(process_model == 1) {  // Underdamped
        Type exp_gdt = exp(-gamma * dt_step);
        Type exp_2gdt = exp(-Type(2.0) * gamma * dt_step);

        // Variance terms
        Type var_x = s2/(gamma*gamma) * (Type(2.0)*gamma*dt_step - Type(3.0) +
          Type(4.0)*exp_gdt - exp_2gdt);
        Type var_v = s2 * (Type(1.0) - exp_2gdt);
        Type cov_xv = s2/gamma * (Type(1.0) - Type(2.0)*exp_gdt + exp_2gdt);

        for(int i = 0; i < 2; i++) {
          // Position mean
          Type mu_x_pred = mu(i,idx) +
            vel(i,idx)/gamma * (Type(1.0) - exp_gdt) +
            s2*h(i)/gamma * (dt_step - (Type(1.0) - exp_gdt)/gamma);

          // Velocity mean
          Type mu_v_pred = vel(i,idx) * exp_gdt +
            s2*h(i)/gamma * (Type(1.0) - exp_gdt);

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
      }

      // TMB oneStepPredict Bivariate Decomposition
      Type varX = cov_obs(0,0);
      Type varY = cov_obs(1,1);
      Type covXY = cov_obs(0,1);

      Type sdX = sqrt(varX);

      Type x_obs = Y(0,i);
      Type y_obs = Y(1,i);
      Type mu_x = mu(0,i);
      Type mu_y = mu(1,i);

      Type kX = keep(0,i);
      Type kY = keep(1,i);

      // 1. Marginal Likelihood of X
      nll -= kX * dnorm(x_obs, mu_x, sdX, true);

      // 2. Conditional Likelihood of Y
      Type mu_y_cond = mu_y + kX * (covXY / varX) * (x_obs - mu_x);
      Type sd_y_cond = sqrt(varY - kX * (covXY * covXY / varX));

      nll -= kY * dnorm(y_obs, mu_y_cond, sd_y_cond, true);
    }
  }

  ADREPORT(beta);
  ADREPORT(sigma);
  if(process_model == 1) ADREPORT(gamma);
  ADREPORT(rho_o);
  ADREPORT(tau);
  ADREPORT(psi);

  // Backup in case sdreport fails
  REPORT(beta);
  REPORT(sigma);
  if(process_model == 1) REPORT(gamma);
  REPORT(rho_o);
  REPORT(tau);
  REPORT(psi);
  REPORT(mu);
  if(process_model == 1) REPORT(vel);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
