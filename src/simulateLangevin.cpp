#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

#define IS_RCPP_BUILD
#include "TMB/include/raster_helpers.hpp"

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

  // Cast raster_coords to an Armadillo matrix to match the header signature
  NumericMatrix raster_coords_rcpp = raster_data["raster_coords"];
  arma::mat raster_coords = Rcpp::as<arma::mat>(raster_coords_rcpp);

  NumericVector raster_resolution = raster_data["raster_resolution"];
  NumericVector all_z_values = raster_data["all_z_values"];
  IntegerVector n_zvals_cov = raster_data["n_zvals_cov"];
  IntegerVector cov_offset = raster_data["cov_offset"];
  int n_covs = raster_data["n_covs"];

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
      arma::mat grad = extract_raster_values(mu_x[idx-1], mu_y[idx-1], time[idx-1],
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
      Named("id") = ID,
      Named("date") = time,
      Named("dt") = dt,
      Named("mu.x") = mu_x,
      Named("mu.y") = mu_y
    );
  } else {
    return DataFrame::create(
      Named("id") = ID,
      Named("date") = time,
      Named("dt") = dt,
      Named("mu.x") = mu_x,
      Named("mu.y") = mu_y,
      Named("v.x") = v_mux,
      Named("v.y") = v_muy
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
    // Convert degrees to radians
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
      Named("id") = data["id"],
                        Named("date") = data["date"],
                                            Named("dt") = data["dt"],
                                                              Named("x") = new_mux,
                                                              Named("y") = new_muy,
                                                              Named("smaj") = M_rand,
                                                              Named("smin") = m_rand,
                                                              Named("eor") = c_rand,
                                                              Named("mu.x") = mux, // true location
                                                              Named("mu.y") = muy); // true location
  } else {
    return DataFrame::create(
      Named("id") = data["id"],
                        Named("date") = data["date"],
                                            Named("dt") = data["dt"],
                                                              Named("x") = new_mux,
                                                              Named("y") = new_muy,
                                                              Named("smaj") = M_rand,
                                                              Named("smin") = m_rand,
                                                              Named("eor") = c_rand,
                                                              Named("mu.x") = mux, // true location
                                                              Named("mu.y") = muy, // true location
                                                              Named("v.x") = data["v.x"],  // true velocity
                                                              Named("v.y") = data["v.y"]); // true velocity
  }
}

// [[Rcpp::export]]
DataFrame measurementError_LS_rcpp(DataFrame data,
                                   double x_sd,
                                   double y_sd,
                                   double tau_x,
                                   double tau_y,
                                   double rho_o,
                                   int model) {

  int n = data.nrows();

  // Get true position vectors
  NumericVector mux = data["mu.x"];
  NumericVector muy = data["mu.y"];

  // Vectors for generated SDs and observed positions
  NumericVector x_sd_vec(n);
  NumericVector y_sd_vec(n);
  NumericVector obs_x(n);
  NumericVector obs_y(n);

  // Generate random SDs and apply measurement error
  for(int i = 0; i < n; i++) {
    // Generate observed standard deviations based on the global x_sd/y_sd
    x_sd_vec[i] = std::abs(R::rnorm(0.0, x_sd));
    y_sd_vec[i] = std::abs(R::rnorm(0.0, y_sd));

    // Calculate covariance matrix components
    double s = tau_x * x_sd_vec[i];
    double q = tau_y * y_sd_vec[i];

    arma::mat cov_obs(2, 2, arma::fill::zeros);
    cov_obs(0,0) = s * s;
    cov_obs(1,1) = q * q;
    cov_obs(0,1) = s * q * rho_o;
    cov_obs(1,0) = cov_obs(0,1);

    // Sample error
    arma::vec mean = {mux[i], muy[i]};
    arma::vec new_pos = rmvnorm(mean, cov_obs);

    obs_x[i] = new_pos(0);
    obs_y[i] = new_pos(1);
  }

  // Create output DataFrame
  if(model == 0){ // Overdamped
    return DataFrame::create(
      Named("id") = data["id"],
                        Named("date") = data["date"],
                                            Named("dt") = data["dt"],
                                                              Named("x") = obs_x,
                                                              Named("y") = obs_y,
                                                              Named("x.sd") = x_sd_vec,
                                                              Named("y.sd") = y_sd_vec,
                                                              Named("mu.x") = mux,
                                                              Named("mu.y") = muy
    );
  } else { // Underdamped
    return DataFrame::create(
      Named("id") = data["id"],
                        Named("date") = data["date"],
                                            Named("dt") = data["dt"],
                                                              Named("x") = obs_x,
                                                              Named("y") = obs_y,
                                                              Named("x.sd") = x_sd_vec,
                                                              Named("y.sd") = y_sd_vec,
                                                              Named("mu.x") = mux,
                                                              Named("mu.y") = muy,
                                                              Named("v.x") = data["v.x"],
                                                              Named("v.y") = data["v.y"]
    );
  }
}
