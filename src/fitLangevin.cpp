// Use 'TMBad' regardless of what's selected by 'TMB::compile'
#undef TMBAD_FRAMEWORK
#undef CPPAD_FRAMEWORK
#define TMBAD_FRAMEWORK

// Just in case we want to use RTMB to inspect tapes...
#undef TMBAD_INDEX_TYPE
#define TMBAD_INDEX_TYPE uint64_t

#include <TMB.hpp>

using namespace density;

template<class Type>
matrix<Type> extract_raster_values(Type x, Type y, Type z,
                                   const array<Type> &raster_vals,
                                   const vector<Type> &all_z_values, // NEW: Flattened vector of ALL times
                                   const vector<int> &n_zvals_cov, 
                                   const vector<int> &cov_offset,  
                                   const matrix<Type> &raster_coords,
                                   const vector<Type> &raster_resolution,
                                   const vector<Type> &raster_extent,
                                   int n_covs) {
  
  int n_cols = CppAD::Integer((raster_extent[1] - raster_extent[0]) / raster_resolution[0]);
  int n_rows = CppAD::Integer((raster_extent[3] - raster_extent[2]) / raster_resolution[1]);
  
  Type x_prop = (x - (raster_extent[0] + raster_resolution[0]/2)) / raster_resolution[0];
  Type y_prop = (y - (raster_extent[2] + raster_resolution[1]/2)) / raster_resolution[1];
  
  int col = CppAD::Integer(x_prop);
  int row = CppAD::Integer(y_prop);
  
  matrix<Type> grad_values(2, n_covs);
  grad_values.setZero();
  
  if(col < 0 || col >= (n_cols-1) || row < 0 || row >= (n_rows-1)) return grad_values;
  
  Type x1 = raster_extent[0] + (col + Type(0.5)) * raster_resolution[0];
  Type x2 = x1 + raster_resolution[0];
  Type y1 = raster_extent[2] + (row + Type(0.5)) * raster_resolution[1];
  Type y2 = y1 + raster_resolution[1];
  
  int rev_row = n_rows - 1 - row;
  
  for(int i = 0; i < n_covs; i++) {
    int layers = n_zvals_cov[i];
    int offset = cov_offset[i];
    
    int z_idx1 = 0, z_idx2 = 0;
    Type z_weight = 0.0;
    
    // --- Z-axis (Time) Interpolation Logic PER COVARIATE ---
    if (layers > 1) {
      // Time values for this specific covariate start at all_z_values[offset]
      if (z <= all_z_values[offset]) {
        z_idx1 = 0; z_idx2 = 0; z_weight = 0.0;
      } else if (z >= all_z_values[offset + layers - 1]) {
        z_idx1 = layers - 1; z_idx2 = layers - 1; z_weight = 0.0;
      } else {
        for (int k = 0; k < layers - 1; k++) {
          Type t1 = all_z_values[offset + k];
          Type t2 = all_z_values[offset + k + 1];
          if (z >= t1 && z <= t2) {
            z_idx1 = k;
            z_idx2 = k + 1;
            z_weight = (z - t1) / (t2 - t1);
            break;
          }
        }
      }
    }
    
    // Evaluate first time slice (or static layer)
    int layer_idx1 = offset + z_idx1;
    int idx1 = layer_idx1 * (n_rows * n_cols) + rev_row * n_cols + col;
    
    Type f1_00 = raster_vals.data()[idx1];
    Type f1_10 = raster_vals.data()[idx1 + 1];
    Type f1_01 = raster_vals.data()[idx1 - n_cols];
    Type f1_11 = raster_vals.data()[idx1 - n_cols + 1];
    
    Type grad_x_1 = ((y2 - y) * (f1_10 - f1_00) + (y - y1) * (f1_11 - f1_01)) / ((y2 - y1) * (x2 - x1));
    Type grad_y_1 = ((x2 - x) * (f1_01 - f1_00) + (x - x1) * (f1_11 - f1_10)) / ((y2 - y1) * (x2 - x1));
    
    // If dynamic and between two distinct time slices, interpolate
    if (layers > 1 && z_idx1 != z_idx2) {
      int layer_idx2 = offset + z_idx2;
      int idx2 = layer_idx2 * (n_rows * n_cols) + rev_row * n_cols + col;
      
      Type f2_00 = raster_vals.data()[idx2];
      Type f2_10 = raster_vals.data()[idx2 + 1];
      Type f2_01 = raster_vals.data()[idx2 - n_cols];
      Type f2_11 = raster_vals.data()[idx2 - n_cols + 1];
      
      Type grad_x_2 = ((y2 - y) * (f2_10 - f2_00) + (y - y1) * (f2_11 - f2_01)) / ((y2 - y1) * (x2 - x1));
      Type grad_y_2 = ((x2 - x) * (f2_01 - f2_00) + (x - x1) * (f2_11 - f2_10)) / ((y2 - y1) * (x2 - x1));
      
      grad_values(0,i) = grad_x_1 * (Type(1.0) - z_weight) + grad_x_2 * z_weight;
      grad_values(1,i) = grad_y_1 * (Type(1.0) - z_weight) + grad_y_2 * z_weight;
    } else {
      // Static covariate OR dynamic evaluated outside bounds
      grad_values(0,i) = grad_x_1;
      grad_values(1,i) = grad_y_1;
    }
  }
  return grad_values;
}

template<class Type>
matrix<Type> calculate_smooth_gradient(Type x, Type y, Type z, 
                                       const array<Type> &raster_vals,
                                       const vector<Type> &all_z_values,
                                       const vector<int> &n_zvals_cov, 
                                       const vector<int> &cov_offset, 
                                       const matrix<Type> &raster_coords,
                                       const vector<Type> &raster_resolution,
                                       const vector<Type> &raster_extent,
                                       int n_covs,
                                       const vector<Type> &weights,
                                       Type sigma,
                                       Type gamma,
                                       Type dt_step,
                                       Type zetaScale,
                                       int model) {
  
  // Calculate scale based on sigma
  Type zeta;
  if(model == 1) {  // Underdamped
    zeta = zetaScale * sqrt(Type(2.0) * Type(M_PI)) * sigma; 
  } else {  // Overdamped
    zeta = zetaScale * sigma / sqrt(Type(2.0));
  }
  
  Type neighborhood_size = zeta * sqrt(dt_step); 
  
  // Calculate gradient at current location
  matrix<Type> smooth_grads = weights(0) * extract_raster_values(x, y, z, raster_vals, 
                                      all_z_values, n_zvals_cov, cov_offset, 
                                      raster_coords, raster_resolution, 
                                      raster_extent, n_covs);
  
  int n_points = weights.size();
  
  // Add gradients for points around center
  for(int i = 0; i < (n_points-1); i++) {
    Type angle;
    if(n_points == 9) {
      // Queen's case: NNE, NEE, SEE, SSE, SSW, SWW, NWW, NNW
      angle = Type(M_PI) * (Type(2.0) * i + Type(1.0)) / Type(8.0);
    } else {
      // Diagonal case: NE, NW, SW, SE
      angle = Type(M_PI) * (Type(2.0) * i + Type(1.0)) / Type(4.0);
    }
    
    smooth_grads += weights(i+1) * extract_raster_values(
      x + neighborhood_size * cos(angle),
      y + neighborhood_size * sin(angle),
      z,
      raster_vals,
      all_z_values,
      n_zvals_cov,
      cov_offset,
      raster_coords,
      raster_resolution,
      raster_extent,
      n_covs);
  }
  
  return smooth_grads;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(model);      // 0 = overdamped, 1 = underdamped
  DATA_MATRIX(Y);           // Matrix of observed locations
  DATA_VECTOR(dt);          // Time steps
  DATA_IVECTOR(isd);        // indexes observations vs. interpolation points
  DATA_IVECTOR(obs_mod);    // indicates which obs error model to be used
  DATA_IVECTOR(ID);         // Track IDs
  DATA_INTEGER(nbStates);   // Number of states
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
  PARAMETER_VECTOR(log_sigma);  // Log velocity process SD
  PARAMETER_MATRIX(beta);       // Matrix of covariate effects [nstates, ncovs]
  
  PARAMETER_MATRIX(mu);        // Locations [2 x timeSteps]
  PARAMETER_MATRIX(v_mu);        // Velocity innovations [2 x timeSteps], only used if model=1
  PARAMETER_VECTOR(log_gamma);  // Log friction coefficient, only used if model=1
  
  // HMM parameters
  PARAMETER_VECTOR(l_delta);    // Initial state distribution
  PARAMETER_MATRIX(l_gamma);    // Transition probabilities
  
  // OBSERVATION PARAMETERS
  // for KF OBS MODEL
  PARAMETER(l_psi); 				    // error SD scaling parameter to account for possible uncertainty in Argos error ellipse variables
  // for LS/GPS OBS MODEL
  PARAMETER_VECTOR(l_tau);     	// error dispersion for LS obs model (log scale)
  PARAMETER(l_rho_o);           // measurement error correlation b/w x,y
  
  Type psi = exp(l_psi);
  vector<Type> tau = exp(l_tau);
  Type rho_o = Type(2.0) / (Type(1.0) + exp(-l_rho_o)) - Type(1.0);
  
  // Transform HMM parameters 
  vector<Type> delta(nbStates);
  matrix<Type> trMat(nbStates,nbStates);
  delta(0) = Type(1.0);
  for(int i=1; i < nbStates; i++){
    delta(i) = exp(l_delta(i-1));
  }
  delta /= delta.sum();
  
  // Pre-calculate log transition probabilities
  matrix<Type> log_trMat(nbStates,nbStates);
  int cpt = 0;
  for(int i=0; i<nbStates; i++){
    for(int j=0; j<nbStates; j++){
      if(i==j) {
        trMat(i,j) = Type(1.0);
        cpt++;
      } else trMat(i,j) = exp(l_gamma(0,i*nbStates+j-cpt));
    }
    trMat.row(i) /= trMat.row(i).sum();
    for(int j=0; j<nbStates; j++){
      log_trMat(i,j) = log(trMat(i,j));
    }
  }
  
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
  vector<Type> sigma(nbStates);
  vector<Type> sigma_sca(nbStates);
  vector<Type> gamma(nbStates);
  for(int i=0; i < nbStates; i++){
    sigma(i) = exp(log_sigma(i) + log(scale_factor));   // Original scale
    sigma_sca(i) = exp(log_sigma(i));
    gamma(i) = exp(log_gamma(i));
  }
  
  // Process model - loop over tracks
  for(int id = 0; id < track_starts.size(); id++) {
    int start_idx = track_starts(id);
    int track_length = track_lengths(id);
    
    for(int state = 0; state < nbStates; state++) {
      
      Type s2 = pow(sigma_sca(state), Type(2.0));
      
      // Add prior on initial velocities for underdamped model only
      if(model == 1) {
        for(int i = 0; i < 2; i++) {
          nll -= dnorm(v_mu(i,start_idx), Type(0.0), sigma_sca(state) / sqrt(Type(2.0) * gamma(state)), true);
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
            n_covs, weights, sigma_sca(state), gamma(state), dt_step, zetaScale, model
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
        h(0) = beta(state,0);
        h(1) = beta(state,0);
        for(int c = 0; c < n_covs; c++) {
          h(0) += beta(state,c+1) * grad(0,c);
          h(1) += beta(state,c+1) * grad(1,c);
        }
        
        if(model == 1) {  // Underdamped
          Type exp_gdt = exp(-gamma(state) * dt_step);
          Type exp_2gdt = exp(-Type(2.0) * gamma(state) * dt_step);
          
          // Variance terms
          Type var_x = s2/(gamma(state)*gamma(state)) * (Type(2.0)*gamma(state)*dt_step - Type(3.0) + 
            Type(4.0)*exp_gdt - exp_2gdt);
          Type var_v = s2 * (Type(1.0) - exp_2gdt);
          Type cov_xv = s2/gamma(state) * (Type(1.0) - Type(2.0)*exp_gdt + exp_2gdt);
          
          for(int i = 0; i < 2; i++) {
            // Position mean
            Type mu_x_pred = mu(i,idx) + 
              v_mu(i,idx)/gamma(state) * (Type(1.0) - exp_gdt) +
              s2*h(i)/gamma(state) * (dt_step - (Type(1.0) - exp_gdt)/gamma(state));
            
            // Velocity mean
            Type mu_v_pred = v_mu(i,idx) * exp_gdt +
              s2*h(i)/gamma(state) * (Type(1.0) - exp_gdt);
            
            // Construct variance-covariance matrix
            matrix<Type> Sigma(2,2);
            Sigma(0,0) = var_x;
            Sigma(1,1) = var_v;
            Sigma(0,1) = cov_xv;
            Sigma(1,0) = cov_xv;
            
            vector<Type> nu(2);
            nu(0) = mu(i,idx+1) - mu_x_pred;
            nu(1) = v_mu(i,idx+1) - mu_v_pred;
            
            nll += MVNORM(Sigma)(nu);
          }
        } else {  // Overdamped
          
          // Expected locations (simplified case of underdamped as γ → ∞)
          Type Ex = x_prev + s2 * dt_step / Type(2.0) * h(0);
          Type Ey = y_prev + s2 * dt_step / Type(2.0) * h(1);
          
          Type sd = sigma_sca(state) * sqrt(dt_step);
          
          nll -= (dnorm(mu(0,idx+1), Ex, sd, true) + 
            dnorm(mu(1,idx+1), Ey, sd, true));
        }
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
  if(model == 1) ADREPORT(gamma);
  ADREPORT(rho_o);
  ADREPORT(tau);
  ADREPORT(psi);
  
  REPORT(beta);
  REPORT(sigma);
  if(model == 1) REPORT(gamma);
  REPORT(mu);
  if(model == 1) REPORT(v_mu);
  
  if(nbStates>1){
    ADREPORT(delta);
    ADREPORT(trMat);
    REPORT(delta);
    REPORT(trMat);
  }
  
  return nll;
}
