// Use 'TMBad' regardless of what's selected by 'TMB::compile'
#undef TMBAD_FRAMEWORK
#undef CPPAD_FRAMEWORK
#define TMBAD_FRAMEWORK

// Just in case we want to use RTMB to inspect tapes...
#undef TMBAD_INDEX_TYPE
#define TMBAD_INDEX_TYPE uint64_t

#include <TMB.hpp>

// Make locally linearized log(sum(exp(.)))
// - sparseDeriv marks expression to have sparse derivatives
// - Turn on/off using 'lseflag_' to compare with exact version
int lseflag_;
extern "C" double Rf_logspace_add(double, double);
TMB_ATOMIC_STATIC_FUNCTION(lse,
                           // Input dim
                           2,
                           // Atomic double
                           ty[0] = Rf_logspace_add(tx[0], tx[1]);,
                           // Atomic reverse
                           // Exact
                           px[0] = exp(tx[0] - ty[0]) * py[0];
                           px[1] = exp(tx[1] - ty[0]) * py[0];
                           if (lseflag_ == 0) {
                             // Force sparse derivatives
                             px[0] = TMBad::sparseDeriv(px[0]);
                             px[1] = TMBad::sparseDeriv(px[1]);
                           }
)
  template<class Type>
  Type lse(Type x, Type y) {
    CppAD::vector<Type> arg(2);
    arg[0] = x;
    arg[1] = y;
    return lse(arg)[0];
  }
  // Override logspace_add
#define logspace_add(x,y) lse(x,y)
  
using namespace density;

/* Numerically stable forward algorithm based on http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf */
template<class Type>
Type forward_alg(vector<Type> delta, matrix<Type> log_trMat, matrix<Type> lnProbs, int nbSteps) {
  int nbStates = log_trMat.cols();
  Type logalpha;
  vector<Type> ldeltaG(nbStates);
  vector<Type> lalpha(nbStates);
  vector<Type> lnewalpha(nbStates);
  
  Type sumalpha  = -INFINITY;
  for(int j=0; j < nbStates; j++){
    //ldeltaG(j) = -INFINITY;
    //for(int i=0; i < nbStates; i++){
    //  ldeltaG(j) = logspace_add(ldeltaG(j),log(delta(i))); 
    //  Rprintf("t 0 state %d ldeltaG %f delta %f log(delta) %f \n",i,asDouble(ldeltaG(j)),asDouble(delta(i)),asDouble(log(delta(i))));
    //}
    lalpha(j) = log(delta(j)) + lnProbs(0,j);
    sumalpha  = logspace_add(sumalpha,lalpha(j));
    //Rprintf("t 0 state %d lalpha %f delta %f l_delta %f lnProbs %f sumalpha %f \n",j,asDouble(lalpha(j)),asDouble(delta(j)),asDouble(log(delta(j))),asDouble(lnProbs(0,j)),asDouble(sumalpha));
  }
  Type jnll = -sumalpha;
  //Rprintf("t 0 delta_1 %f delta_2 %f sumalpha %f jnll %f \n",asDouble(delta(0)),asDouble(delta(1)),asDouble(sumalpha),asDouble(jnll));
  lalpha -= sumalpha;
  for(int t=1; t < nbSteps; t++){
    sumalpha = -INFINITY;
    for(int j=0; j < nbStates; j++){
      logalpha = -INFINITY;
      for(int i=0; i < nbStates; i++){
        logalpha = logspace_add(logalpha,lalpha(i)+log_trMat(i,j));     // does not recognize sparsity pattern even when nbStates=1
        //logalpha = lalpha(i)+log_trMat(i,j);                          // recognizes sparsity pattern (jnll still correct for nbStates=1 but not for nbStates>1)
        //logalpha = log(exp(logalpha)+exp(lalpha(i)+log_trMat(i,j)));  // doesn't fail but does not recognize sparsity pattern even when nbStates=1
        //Rprintf("t %d j %d i %d logalpha %f log_trMat %f \n",t,j,i,asDouble(logalpha),asDouble(log_trMat(i,j)));
      }
      lnewalpha(j) = logalpha + lnProbs(t,j);
      sumalpha  = logspace_add(sumalpha,lnewalpha(j));               // does not recognize sparsity pattern even when nbStates=1
      //Rprintf("t %d j %d lnewalpha %f logalpha %f lnProbs %f sumalpha %f \n",t,j,asDouble(lnewalpha(j)),asDouble(logalpha),asDouble(lnProbs(t,j)),asDouble(sumalpha));
      //sumalpha  = lnewalpha(j);                                    // recognizes sparsity pattern (jnll still correct for nbStates=1 but not for nbStates>1)
      //sumalpha = log(exp(sumalpha)+exp(lnewalpha(j)));             // fails due to NaN
    }
    jnll -= sumalpha;
    lalpha = lnewalpha - sumalpha;
    //Rprintf("t %d sumalpha %f jnll %f \n",t,asDouble(sumalpha),asDouble(jnll));
  }
  return jnll;
}

template<class Type>
matrix<Type> extract_raster_values(Type x, Type y, 
                                   const array<Type> &raster_vals,
                                   const matrix<Type> &raster_coords,
                                   const vector<Type> &raster_resolution,
                                   const vector<Type> &raster_extent,
                                   int n_covs) {
  // Calculate raster dimensions
  int n_cols = CppAD::Integer((raster_extent[1] - raster_extent[0]) / raster_resolution[0]);
  int n_rows = CppAD::Integer((raster_extent[3] - raster_extent[2]) / raster_resolution[1]);
  
  // Calculate cell indices
  Type x_prop = (x - (raster_extent[0] + raster_resolution[0]/2)) / raster_resolution[0];
  Type y_prop = (y - (raster_extent[2] + raster_resolution[1]/2)) / raster_resolution[1];
  
  int col = CppAD::Integer(x_prop);
  int row = CppAD::Integer(y_prop);
  
  matrix<Type> grad_values(2, n_covs);
  grad_values.setZero();
  
  // Bounds checking
  if(col < 0 || col >= (n_cols-1) || row < 0 || row >= (n_rows-1)) {
    return grad_values;
  }
  
  // Get cell coordinates (matching collapseRaster's cell centers)
  Type x1 = raster_extent[0] + (col + Type(0.5)) * raster_resolution[0];
  Type x2 = x1 + raster_resolution[0];
  Type y1 = raster_extent[2] + (row + Type(0.5)) * raster_resolution[1];
  Type y2 = y1 + raster_resolution[1];
  
  for(int i = 0; i < n_covs; i++) {
    matrix<Type> f(2,2);
    
    // Need to transform indices to match t(apply(as.matrix(rast), 2, rev))
    int rev_row = n_rows - 1 - row;
    
    // Calculate linear indices for 3D array [layer, row, col]
    int idx = i * (n_rows * n_cols) + rev_row * n_cols + col;
    
    // Access values in the transformed space
    f(0,0) = raster_vals.data()[idx];                          // bottom-left
    f(1,0) = raster_vals.data()[idx + 1];                      // bottom-right
    f(0,1) = raster_vals.data()[idx - n_cols];                 // top-left
    f(1,1) = raster_vals.data()[idx - n_cols + 1];             // top-right
    
    // Calculate gradients exactly as in biGrad
    grad_values(0,i) = ((y2 - y) * (f(1,0) - f(0,0)) + 
      (y - y1) * (f(1,1) - f(0,1))) / 
      ((y2 - y1) * (x2 - x1));
    
    grad_values(1,i) = ((x2 - x) * (f(0,1) - f(0,0)) + 
      (x - x1) * (f(1,1) - f(1,0))) / 
      ((y2 - y1) * (x2 - x1));
  }
  return grad_values;
}

template<class Type>
matrix<Type> calculate_smooth_gradient(Type x, Type y, 
                                       const array<Type> &raster_vals,
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
    zeta = zetaScale * sqrt(Type(2.0) * Type(M_PI)) * sigma; // Type zeta = sigma / sqrt(gamma);
  } else {  // Overdamped
    zeta = zetaScale * sigma / sqrt(Type(2.0));
  }
  
  Type neighborhood_size = zeta * sqrt(dt_step); // would use: neighborhood_size = zeta * sqrt(Type(2.0) * dt_step) to match Blackwell & Matthiopoulos (2024) when zeta is calculated from mu
  
  // Calculate gradient at current location
  matrix<Type> smooth_grads = weights(0) * extract_raster_values(x, y, raster_vals, 
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
      raster_vals,
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
  DATA_ARRAY(raster_vals);         // 3D array of raster values [ncov, nrow, ncol]
  DATA_MATRIX(raster_coords);      // Matrix of raster cell coordinates
  DATA_VECTOR(raster_resolution);  // Vector with x,y resolution
  DATA_VECTOR(raster_extent);      // Vector with xmin,xmax,ymin,ymax
  DATA_INTEGER(n_covs);           // Number of covariates
  DATA_INTEGER(smoothGradient);  // Boolean flag for using smooth gradient calculation
  DATA_VECTOR(weights);          // Vector of 5 (diagonal) or 9 (queen's) weights that sum to 1
  DATA_SCALAR(zetaScale);        // scale factor for smooth gradient neighborhood (>1 increases, <1 decreases)
  
  // for KF observation model
  DATA_VECTOR(m);             // m is the semi-minor axis length
  DATA_VECTOR(M);             // M is the semi-major axis length
  DATA_VECTOR(c);             // c is the orientation of the error ellipse
  // for LS observation model
  DATA_MATRIX(K);            // error weighting factors for LS obs model
  
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
  PARAMETER(l_psi);            // error SD scaling parameter
  PARAMETER_VECTOR(l_tau);      // error dispersion for LS obs model (log scale)
  PARAMETER(l_rho_o);          // error correlation
  
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
    
    //matrix<Type> allProbs(track_length-1, nbStates);
    //allProbs.setZero();
    
    for(int state = 0; state < nbStates; state++) {
      
      Type s2 = pow(sigma_sca(state), Type(2.0));
      
      // Add prior on initial velocities for underdamped model only
      if(model == 1) {
        for(int i = 0; i < 2; i++) {
          nll -= dnorm(v_mu(i,start_idx), Type(0.0), sigma_sca(state) / sqrt(Type(2.0) * gamma(state)), true); // sigma_sca(state), true);
        }
      }
      
      for(int t = 0; t < (track_length-1); t++) {
        
        int idx = start_idx + t;
        Type dt_step = dt(idx+1); 
        
        // Get current locations
        Type x_prev = mu(0,idx);
        Type y_prev = mu(1,idx);
        
        // Calculate gradient of log(π)
        matrix<Type> grad;
        if(smoothGradient) {
          grad = calculate_smooth_gradient(
            x_prev, y_prev,
            raster_vals,
            raster_coords,
            raster_resolution,
            raster_extent,
            n_covs,
            weights,
            sigma_sca(state),
            gamma(state),
            dt_step,
            zetaScale,
            model
          );
        } else {
          grad = extract_raster_values(
            x_prev, y_prev,
            raster_vals,
            raster_coords,
            raster_resolution,
            raster_extent,
            n_covs
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
          Type var_x = s2/(gamma(state)*gamma(state)) * 
            (Type(2.0)*gamma(state)*dt_step - Type(3.0) + 
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
            
            //allProbs(t,state) -= MVNORM(Sigma)(nu);
            
            nll += MVNORM(Sigma)(nu);
          }
        } else {  // Overdamped
          
          // Expected locations (simplified case of underdamped as γ → ∞)
          Type Ex = x_prev + s2 * dt_step / Type(2.0) * h(0);
          Type Ey = y_prev + s2 * dt_step / Type(2.0) * h(1);
          
          Type sd = sigma_sca(state) * sqrt(dt_step);
          
          // Add to state probability
          //allProbs(t,state) += dnorm(x_curr, Ex, sd, true) + 
          //  dnorm(y_curr, Ey, sd, true);
          
          nll -= (dnorm(mu(0,idx+1), Ex, sd, true) + 
            dnorm(mu(1,idx+1), Ey, sd, true));
        }
      }
    }
    //nll += forward_alg(delta, log_trMat, allProbs, track_length-1);
  }
  
  
  // OBSERVATION MODEL
  matrix<Type> cov_obs(2, 2);
  cov_obs.setZero();
  for(int i = 0; i < timeSteps; ++i) {
    if(isd(i) == 1) {
      if(obs_mod(i) == 0) {
        // Argos Least Squares observations
        Type s = tau(0) * K(i,0);
        Type q = tau(1) * K(i,1);
        cov_obs(0,0) = s * s;
        cov_obs(1,1) = q * q;
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
  //ADREPORT(mu);
  //ADREPORT(v_mu);
  
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
  