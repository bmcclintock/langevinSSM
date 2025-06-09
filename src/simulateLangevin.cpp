#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

NumericMatrix extract_raster_values(double x, double y, 
                                    const NumericVector& raster_vals,
                                    const NumericMatrix& raster_coords,
                                    const NumericVector& raster_resolution,
                                    const NumericVector& raster_extent,
                                    int n_covs) {
  // Calculate raster dimensions
  int n_cols = (raster_extent[1] - raster_extent[0]) / raster_resolution[0];
  int n_rows = (raster_extent[3] - raster_extent[2]) / raster_resolution[1];
  
  // Calculate cell indices exactly as in TMB code
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
  
  // Get cell coordinates exactly as in TMB code
  double x1 = raster_extent[0] + (col + 0.5) * raster_resolution[0];
  double x2 = x1 + raster_resolution[0];
  double y1 = raster_extent[2] + (row + 0.5) * raster_resolution[1];
  double y2 = y1 + raster_resolution[1];
  
  // Transform row index to match TMB code
  int rev_row = n_rows - 1 - row;
  
  for(int i = 0; i < n_covs; i++) {
    // Calculate linear indices for 3D array [layer, row, col]
    int idx = i * (n_rows * n_cols) + rev_row * n_cols + col;
    
    // Access values in the transformed space exactly as in TMB code
    double f00 = raster_vals[idx];                          // bottom-left
    double f10 = raster_vals[idx + 1];                      // bottom-right
    double f01 = raster_vals[idx - n_cols];                 // top-left
    double f11 = raster_vals[idx - n_cols + 1];             // top-right
    
    // Calculate gradients exactly as in TMB code
    grad_values(0,i) = ((y2 - y) * (f10 - f00) + 
      (y - y1) * (f11 - f01)) / 
      ((y2 - y1) * (x2 - x1));
    
    grad_values(1,i) = ((x2 - x) * (f01 - f00) + 
      (x - x1) * (f11 - f10)) / 
      ((y2 - y1) * (x2 - x1));
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
                                  double lambda,
                                  double gamma,
                                  double sigma,
                                  NumericVector beta,
                                  List raster_data,
                                  NumericMatrix initialPosition,
                                  double min_dt = 4e-5) {
  
  // Extract raster extent early for boundary checking
  NumericVector raster_extent = raster_data["raster_extent"];
  double xmin = raster_extent[0];
  double xmax = raster_extent[1];
  double ymin = raster_extent[2];
  double ymax = raster_extent[3];
  
  // Function to check if position is within bounds
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
  
  // Extract raster information
  NumericVector raster_vals = raster_data["raster_vals"];
  NumericMatrix raster_coords = raster_data["raster_coords"];
  NumericVector raster_resolution = raster_data["raster_resolution"];
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
  
  double s2 = sigma * sigma;
  
  for(int i = 0; i < nbAnimals; i++) {
    int start_idx = i * obsPerAnimal;
    
    // Generate wait times
    NumericVector waitTimes(obsPerAnimal-1);
    for(int j = 0; j < obsPerAnimal-1; j++) {
      double wait;
      do {
        wait = R::rexp(lambda);
      } while(wait < min_dt);
      waitTimes[j] = wait;
    }
    
    // Initialize first observation
    ID[start_idx] = i + 1;
    time[start_idx] = 0;
    dt[start_idx] = 0;
    mu_x[start_idx] = initialPosition(i,0);
    mu_y[start_idx] = initialPosition(i,1);
    v_mux[start_idx] = R::rnorm(0, sigma);
    v_muy[start_idx] = R::rnorm(0, sigma);
    
    // Simulate remaining observations
    for(int t = 1; t < obsPerAnimal; t++) {
      int idx = start_idx + t;
      ID[idx] = i + 1;
      time[idx] = time[idx-1] + waitTimes[t-1];
      dt[idx] = waitTimes[t-1];
      
      double dt_step = dt[idx];
      
      // Calculate gradients
      NumericMatrix grad = extract_raster_values(mu_x[idx-1], mu_y[idx-1],
                                                 raster_vals,
                                                 raster_coords,
                                                 raster_resolution,
                                                 raster_extent,
                                                 n_covs);
      
      // Calculate force vector
      NumericVector h(2, 0.0);  // Initialize to zero
      for(int c = 0; c < n_covs; c++) {
        h[0] += beta[c] * grad(0,c);
        h[1] += beta[c] * grad(1,c);
      }
      
      if(model == 0){  // original (overdamped) Langevin diffusion
        
        // Calculate means
        arma::vec mean(2);
        mean(0) = mu_x[idx-1] + s2*dt_step/2.0 * h[0];
        mean(1) = mu_y[idx-1] + s2*dt_step/2.0 * h[1];
        
        // Calculate covariance matrix
        double sd = sigma * sqrt(dt_step);
        
        // Generate new positions and velocities jointly
        arma::vec new_state(2);
        new_state(0) = R::rnorm(mean(0), sd);
        new_state(1) = R::rnorm(mean(1), sd);
        
        // Check if new position is within bounds before assigning
        check_bounds(new_state(0), new_state(1), i+1, time[idx]);
        
        mu_x[idx] = new_state(0);
        mu_y[idx] = new_state(1);
        
      } else if(model == 1){    // underdamped Langevin diffusion
        
        double exp_gdt = exp(-gamma * dt_step);
        double exp_2gdt = exp(-2 * gamma * dt_step);
        
        // Calculate means
        arma::vec mean(4);
        mean(0) = mu_x[idx-1] + v_mux[idx-1]/gamma * (1 - exp_gdt) +
          s2*h[0]/gamma * (dt_step - (1 - exp_gdt)/gamma);
        mean(1) = v_mux[idx-1] * exp_gdt + s2*h[0]/gamma * (1 - exp_gdt);
        mean(2) = mu_y[idx-1] + v_muy[idx-1]/gamma * (1 - exp_gdt) +
          s2*h[1]/gamma * (dt_step - (1 - exp_gdt)/gamma);
        mean(3) = v_muy[idx-1] * exp_gdt + s2*h[1]/gamma * (1 - exp_gdt);
        
        // Calculate covariance matrix
        double var_x = s2/(gamma*gamma) * (2*gamma*dt_step - 3 + 4*exp_gdt - exp_2gdt);
        double var_v = s2 * (1 - exp_2gdt);
        double cov_xv = s2/gamma * (1 - 2*exp_gdt + exp_2gdt);
        
        arma::mat Sigma(4,4, arma::fill::zeros);
        Sigma(0,0) = Sigma(2,2) = var_x;
        Sigma(1,1) = Sigma(3,3) = var_v;
        Sigma(0,1) = Sigma(1,0) = cov_xv;
        Sigma(2,3) = Sigma(3,2) = cov_xv;
        
        // Generate new positions and velocities jointly
        arma::vec new_state = rmvnorm(mean, Sigma);
        
        // Check if new position is within bounds before assigning
        check_bounds(new_state(0), new_state(2), i+1, time[idx]);
        
        mu_x[idx] = new_state(0);
        v_mux[idx] = new_state(1);
        mu_y[idx] = new_state(2);
        v_muy[idx] = new_state(3);
      }
    }
  }
  
  if(model == 0){
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
