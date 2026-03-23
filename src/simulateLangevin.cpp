#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

NumericMatrix extract_raster_values(double x, double y, double z,
                                    const NumericVector& raster_vals,
                                    const NumericVector& all_z_values,
                                    const IntegerVector& n_zvals_cov,
                                    const IntegerVector& cov_offset,
                                    const NumericMatrix& raster_coords,
                                    const NumericVector& raster_resolution,
                                    const NumericVector& raster_extent,
                                    int n_covs) {
  // Calculate raster dimensions
  int n_cols = (raster_extent[1] - raster_extent[0]) / raster_resolution[0];
  int n_rows = (raster_extent[3] - raster_extent[2]) / raster_resolution[1];
  
  // Calculate cell indices
  double x_prop = (x - (raster_extent[0] + raster_resolution[0]/2)) / raster_resolution[0];
  double y_prop = (y - (raster_extent[2] + raster_resolution[1]/2)) / raster_resolution[1];
  
  int col = floor(x_prop);
  int row = floor(y_prop);
  
  NumericMatrix grad_values(2, n_covs);
  std::fill(grad_values.begin(), grad_values.end(), 0.0);
  
  // Bounds checking
  if(col < 0 || col >= (n_cols-1) || row < 0 || row >= (n_rows-1)) {
    return grad_values;
  }
  
  // Get cell coordinates
  double x1 = raster_extent[0] + (col + 0.5) * raster_resolution[0];
  double x2 = x1 + raster_resolution[0];
  double y1 = raster_extent[2] + (row + 0.5) * raster_resolution[1];
  double y2 = y1 + raster_resolution[1];
  
  // Transform row index
  int rev_row = n_rows - 1 - row;
  
  for(int i = 0; i < n_covs; i++) {
    int layers = n_zvals_cov[i];
    int offset = cov_offset[i];
    
    int z_idx1 = 0, z_idx2 = 0;
    double z_weight = 0.0;
    
    // --- Z-axis (Time) Interpolation Logic PER COVARIATE ---
    if (layers > 1) {
      // Time values for this specific covariate start at all_z_values[offset]
      if (z <= all_z_values[offset]) {
        z_idx1 = 0; z_idx2 = 0; z_weight = 0.0;
      } else if (z >= all_z_values[offset + layers - 1]) {
        z_idx1 = layers - 1; z_idx2 = layers - 1; z_weight = 0.0;
      } else {
        for (int k = 0; k < layers - 1; k++) {
          double t1 = all_z_values[offset + k];
          double t2 = all_z_values[offset + k + 1];
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
    
    double f1_00 = raster_vals[idx1];
    double f1_10 = raster_vals[idx1 + 1];
    double f1_01 = raster_vals[idx1 - n_cols];
    double f1_11 = raster_vals[idx1 - n_cols + 1];
    
    double grad_x_1 = ((y2 - y) * (f1_10 - f1_00) + (y - y1) * (f1_11 - f1_01)) / ((y2 - y1) * (x2 - x1));
    double grad_y_1 = ((x2 - x) * (f1_01 - f1_00) + (x - x1) * (f1_11 - f1_10)) / ((y2 - y1) * (x2 - x1));
    
    // If dynamic and between two distinct time slices, interpolate
    if (layers > 1 && z_idx1 != z_idx2) {
      int layer_idx2 = offset + z_idx2;
      int idx2 = layer_idx2 * (n_rows * n_cols) + rev_row * n_cols + col;
      
      double f2_00 = raster_vals[idx2];
      double f2_10 = raster_vals[idx2 + 1];
      double f2_01 = raster_vals[idx2 - n_cols];
      double f2_11 = raster_vals[idx2 - n_cols + 1];
      
      double grad_x_2 = ((y2 - y) * (f2_10 - f2_00) + (y - y1) * (f2_11 - f2_01)) / ((y2 - y1) * (x2 - x1));
      double grad_y_2 = ((x2 - x) * (f2_01 - f2_00) + (x - x1) * (f2_11 - f2_10)) / ((y2 - y1) * (x2 - x1));
      
      grad_values(0,i) = grad_x_1 * (1.0 - z_weight) + grad_x_2 * z_weight;
      grad_values(1,i) = grad_y_1 * (1.0 - z_weight) + grad_y_2 * z_weight;
    } else {
      // Static covariate OR dynamic evaluated outside bounds
      grad_values(0,i) = grad_x_1;
      grad_values(1,i) = grad_y_1;
    }
  }
  return grad_values;
}

// Helper function for multivariate normal sampling 
arma::vec rmvnorm(const arma::vec& mean, const arma::mat& sigma) {
  int n = mean.n_elem;
  arma::vec z(n);
  for(int i = 0; i < n; i++) {
    z(i) = R::rnorm(0, 1);
  }
  return mean + arma::trans(arma::chol(sigma)) * z;
}

// [[Rcpp::export]]
DataFrame simulate_langevin_cpp(int model,
                                int nbAnimals, 
                                int obsPerAnimal,
                                double timeStep,
                                double gamma,
                                double sigma,
                                NumericVector beta,
                                List raster_data,
                                NumericMatrix initialPosition) {
  
  // Extract raster extent for boundary checking
  NumericVector raster_extent = raster_data["raster_extent"];
  double xmin = raster_extent[0];
  double xmax = raster_extent[1];
  double ymin = raster_extent[2];
  double ymax = raster_extent[3];
  
  // Bounds checking function
  auto check_bounds = [xmin, xmax, ymin, ymax](double x, double y, int animal_id, double time) {
    if(x < xmin || x > xmax || y < ymin || y > ymax) {
      std::string err_msg = "Position out of bounds for animal " + 
        std::to_string(animal_id) + 
        " at time " + std::to_string(time) + 
        "\nPosition: (" + std::to_string(x) + ", " + std::to_string(y) + ")" +
        "\nBounds: x[" + std::to_string(xmin) + ", " + std::to_string(xmax) + "]" +
        " y[" + std::to_string(ymin) + ", " + std::to_string(ymax) + "]";
      throw Rcpp::exception(err_msg.c_str());
    }
  };
  
  // Check initial positions
  for(int i = 0; i < nbAnimals; i++) {
    check_bounds(initialPosition(i,0), initialPosition(i,1), i+1, 0.0);
  }
  
  // Extract raster information and time metadata
  NumericVector raster_vals = raster_data["raster_vals"];
  NumericMatrix raster_coords = raster_data["raster_coords"];
  NumericVector raster_resolution = raster_data["raster_resolution"];
  NumericVector all_z_values = raster_data["all_z_values"];
  IntegerVector n_zvals_cov = raster_data["n_zvals_cov"];
  IntegerVector cov_offset = raster_data["cov_offset"];
  int n_covs = raster_data["n_covs"];
  
  // Create output dataframe
  int total_obs = nbAnimals * obsPerAnimal;
  NumericVector ID(total_obs);
  NumericVector time(total_obs);
  NumericVector dt(total_obs);
  NumericVector mu_x(total_obs);
  NumericVector mu_y(total_obs);
  NumericVector v_mux(total_obs);
  NumericVector v_muy(total_obs);
  
  for(int i = 0; i < nbAnimals; i++) {
    int start_idx = i * obsPerAnimal;
    
    // Initialize first observation
    ID[start_idx] = i + 1;
    time[start_idx] = 0;
    dt[start_idx] = 0;
    mu_x[start_idx] = initialPosition(i,0);
    mu_y[start_idx] = initialPosition(i,1);
    v_mux[start_idx] = R::rnorm(0, sigma / sqrt(2. * gamma));
    v_muy[start_idx] = R::rnorm(0, sigma / sqrt(2. * gamma));
    
    double s2 = sigma * sigma;
    
    // Simulate remaining observations
    for(int t = 1; t < obsPerAnimal; t++) {
      int idx = start_idx + t;
      ID[idx] = i + 1;
      time[idx] = time[idx-1] + timeStep; 
      dt[idx] = timeStep; 
      
      double dt_step = dt[idx];
      
      // Calculate gradients using previous position AND previous time
      NumericMatrix grad = extract_raster_values(mu_x[idx-1], mu_y[idx-1], time[idx-1],
                                                 raster_vals,
                                                 all_z_values,
                                                 n_zvals_cov,
                                                 cov_offset,
                                                 raster_coords,
                                                 raster_resolution,
                                                 raster_extent,
                                                 n_covs);
      
      // Calculate force vector
      NumericVector h(2, 0.0);
      for(int c = 0; c < n_covs; c++) {
        h[0] += beta[c] * grad(0,c);
        h[1] += beta[c] * grad(1,c);
      }
      
      if(model == 0) {  // overdamped Langevin
        // Calculate means with scaled parameters
        arma::vec mean(2);
        mean(0) = mu_x[idx-1] + s2*dt_step/2.0 * h[0];
        mean(1) = mu_y[idx-1] + s2*dt_step/2.0 * h[1];
        
        // Calculate scaled standard deviation
        double sd = sigma * sqrt(dt_step);
        
        // Generate new positions
        arma::vec new_state(2);
        new_state(0) = R::rnorm(mean(0), sd);
        new_state(1) = R::rnorm(mean(1), sd);
        
        check_bounds(new_state(0), new_state(1), i+1, time[idx]);
        
        mu_x[idx] = new_state(0);
        mu_y[idx] = new_state(1);
        
      } else if(model == 1) {  // underdamped Langevin
        double exp_gdt = exp(-gamma * dt_step);
        double exp_2gdt = exp(-2 * gamma * dt_step);
        
        // Calculate means with scaled parameters
        arma::vec mean(4);
        mean(0) = mu_x[idx-1] + v_mux[idx-1]/gamma * (1 - exp_gdt) +
          s2*h[0]/gamma * (dt_step - (1 - exp_gdt)/gamma);
        mean(1) = v_mux[idx-1] * exp_gdt + s2*h[0]/gamma * (1 - exp_gdt);
        mean(2) = mu_y[idx-1] + v_muy[idx-1]/gamma * (1 - exp_gdt) +
          s2*h[1]/gamma * (dt_step - (1 - exp_gdt)/gamma);
        mean(3) = v_muy[idx-1] * exp_gdt + s2*h[1]/gamma * (1 - exp_gdt);
        
        // Calculate covariance matrix with scaled sigma
        double var_x = s2/(gamma*gamma) * (2*gamma*dt_step - 3 + 4*exp_gdt - exp_2gdt);
        double var_v = s2 * (1 - exp_2gdt);
        double cov_xv = s2/gamma * (1 - 2*exp_gdt + exp_2gdt);
        
        arma::mat Sigma(4,4, arma::fill::zeros);
        Sigma(0,0) = Sigma(2,2) = var_x;
        Sigma(1,1) = Sigma(3,3) = var_v;
        Sigma(0,1) = Sigma(1,0) = cov_xv;
        Sigma(2,3) = Sigma(3,2) = cov_xv;
        
        // Generate new state
        arma::vec new_state = rmvnorm(mean, Sigma);
        
        check_bounds(new_state(0), new_state(2), i+1, time[idx]);
        
        mu_x[idx] = new_state(0);
        v_mux[idx] = new_state(1);
        mu_y[idx] = new_state(2);
        v_muy[idx] = new_state(3);
      }
    }
  }
  
  if(model == 0) {
    return DataFrame::create(
      Named("ID") = ID,
      Named("time") = time,
      Named("dt") = dt,
      Named("mu.x") = mu_x,
      Named("mu.y") = mu_y
    );
  } else {
    return DataFrame::create(
      Named("ID") = ID,
      Named("time") = time,
      Named("dt") = dt,
      Named("mu.x") = mu_x,
      Named("mu.y") = mu_y,
      Named("v_mux") = v_mux,
      Named("v_muy") = v_muy
    );
  }
}

// [[Rcpp::export]]
DataFrame measurementError_rcpp(DataFrame data,
                                double M,
                                double m,
                                NumericVector c,
                                double psi,
                                int model) {
  
  // Get dimensions
  int n = data.nrows();
  
  // Generate random values for error parameters
  NumericVector M_rand(n);
  NumericVector m_rand(n);
  NumericVector c_rand(n);
  
  for(int i = 0; i < n; i++) {
    M_rand[i] = abs(R::rnorm(0.0, M));
    m_rand[i] = abs(R::rnorm(0.0, m));
    if(M_rand[i] < m_rand[i]){
      double tmpM = M_rand[i];
      double tmpm = m_rand[i];
      M_rand[i] = tmpm;
      m_rand[i] = tmpM;
    }
    // Convert degrees to radians (like momentuHMM:::radian)
    c_rand[i] = R::runif(c[0], c[1]) * M_PI / 180.0;
  }
  
  // Constants
  double z = sqrt(2.0);
  
  // Get position vectors
  NumericVector mux = data["mu.x"];
  NumericVector muy = data["mu.y"];
  
  // Create new vectors for results
  NumericVector new_mux = clone(mux);
  NumericVector new_muy = clone(muy);
  
  // Calculate error parameters
  NumericVector s2c(n);
  NumericVector c2c(n);
  NumericVector M2(n);
  NumericVector m2(n);
  
  for(int i = 0; i < n; i++) {
    s2c[i] = sin(c_rand[i]) * sin(c_rand[i]);
    c2c[i] = cos(c_rand[i]) * cos(c_rand[i]);
    M2[i] = (M_rand[i] / z) * (M_rand[i] / z);
    m2[i] = (m_rand[i] * psi / z) * (m_rand[i] * psi / z);
  }
  
  // Apply measurement error
  for(int i = 0; i < n; i++) {
    arma::mat cov_obs(2, 2, arma::fill::zeros);
    
    cov_obs(0,0) = M2[i] * s2c[i] + m2[i] * c2c[i];
    cov_obs(1,1) = M2[i] * c2c[i] + m2[i] * s2c[i];
    cov_obs(0,1) = (0.5 * (M_rand[i] * M_rand[i] - 
      (m_rand[i] * psi * m_rand[i] * psi))) * cos(c_rand[i]) * sin(c_rand[i]);
    cov_obs(1,0) = cov_obs(0,1);
    
    arma::vec mean = {mux[i], muy[i]};
    arma::vec new_pos = rmvnorm(mean, cov_obs);
    
    new_mux[i] = new_pos(0);
    new_muy[i] = new_pos(1);
  }
  
  // Create output DataFrame
  if(model==0){
    return DataFrame::create(
      Named("ID") = data["ID"],
                        Named("time") = data["time"],
                                            Named("dt") = data["dt"],
                                                              Named("mu.x") = new_mux,
                                                              Named("mu.y") = new_muy,
                                                              Named("error_semimajor_axis") = M_rand,
                                                              Named("error_semiminor_axis") = m_rand,
                                                              Named("error_ellipse_orientation") = c_rand,
                                                              Named("mux") = mux, // true location
                                                              Named("muy") = muy); // true location   
  } else {
    return DataFrame::create(
      Named("ID") = data["ID"],
                        Named("time") = data["time"],
                                            Named("dt") = data["dt"],
                                                              Named("mu.x") = new_mux,
                                                              Named("mu.y") = new_muy,
                                                              Named("error_semimajor_axis") = M_rand,
                                                              Named("error_semiminor_axis") = m_rand,
                                                              Named("error_ellipse_orientation") = c_rand,
                                                              Named("mux") = mux, // true location
                                                              Named("muy") = muy, // true location   
                                                              Named("v_mux") = data["v_mux"],  // true velocity
                                                                                   Named("v_muy") = data["v_muy"]); // true velocity
  }
}
