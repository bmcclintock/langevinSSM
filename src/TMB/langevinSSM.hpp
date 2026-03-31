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
  DATA_MATRIX(Y);           // Matrix of observed locations
  DATA_VECTOR(dt);          // Time steps
  DATA_IVECTOR(isd);        // indexes observations vs. interpolation points
  DATA_IVECTOR(obs_mod);    // indicates which obs error model to be used
  DATA_IVECTOR(ID);         // Track IDs
  DATA_IVECTOR(nbObs);      // Number of observations per step
  DATA_SCALAR(scale_factor);

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

  // for KF observation model
  DATA_VECTOR(m);                 //  m is the semi-minor axis length
  DATA_VECTOR(M);                 //  M is the semi-major axis length
  DATA_VECTOR(c);                 //  c is the orientation of the error ellipse
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
      //h(0) = beta(0);
      //h(1) = beta(0);
      for(int c = 0; c < n_covs; c++) {
        h(0) += beta(c) * grad(0,c);
        h(1) += beta(c) * grad(1,c);
        //h(0) += beta(c+1) * grad(0,c);
        //h(1) += beta(c+1) * grad(1,c);
      }

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


  // OBSERVATION MODEL
  matrix<Type> cov_obs(2, 2);
  cov_obs.setZero();
  for(int i = 0; i < timeSteps; ++i) {
    if(isd(i) == 1) {
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
        Type s2c = sin(c(i)) * sin(c(i));
        Type c2c = cos(c(i)) * cos(c(i));
        Type M2  = (M(i) / z) * (M(i) / z);
        Type m2 = (m(i) * psi / z) * (m(i) * psi / z);
        cov_obs(0,0) = (M2 * s2c + m2 * c2c);
        cov_obs(1,1) = (M2 * c2c + m2 * s2c);
        cov_obs(0,1) = (0.5 * (M(i) * M(i) - (m(i) * psi * m(i) * psi))) * cos(c(i)) * sin(c(i));
        cov_obs(1,0) = cov_obs(0,1);
      }
      nll += MVNORM(cov_obs)(Y.col(i) - mu.col(i));
    }
  }

  ADREPORT(beta);
  ADREPORT(sigma);
  if(process_model == 1) ADREPORT(gamma);
  ADREPORT(rho_o);
  ADREPORT(tau);
  ADREPORT(psi);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
