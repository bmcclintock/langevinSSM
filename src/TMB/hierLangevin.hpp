// src/TMB/hierLangevin.hpp
//
// Stage II model for the multistage (Conditionally Independent Hierarchical
// Model, CIHM) fitting approach of Johnson, Brost & Hooten (2022, JABES
// 27:382-400), "Greater Than the Sum of its Parts: Computationally Flexible
// Bayesian Hierarchical Modeling".

#ifndef hierLangevin_hpp
#define hierLangevin_hpp

using namespace density;

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type hierLangevin(objective_function<Type>* obj)
{
  DATA_MATRIX(theta_hat);      // n x p matrix of Stage I point estimates (working scale)
  DATA_ARRAY(Shat);            // p x p x n array of Stage I covariance matrices (fixed/known)
  DATA_INTEGER(use_mpl);       // 0/1: apply MPL log-gamma penalty to log_sd
  DATA_VECTOR(mpl_alpha);      // length-p MPL shape parameters (only used if use_mpl==1)
  // Natural-scale reporting: col_type(j) = 0 means parameter j is already on
  // its natural/additive scale (e.g. a habitat-selection coefficient beta);
  // col_type(j) = 1 means parameter j is a log-scale movement parameter
  // (e.g. log_sigma, log_gamma) whose natural-scale population mean and
  // individual-level predictions require exponentiation (with an optional
  // additive shift col_shift(j), used for log_sigma's log(scaleFactor) term).
  DATA_IVECTOR(col_type);      // length p, 0 = identity, 1 = exponential
  DATA_VECTOR(col_shift);      // length p, additive shift applied before exp() for type==1

  PARAMETER_VECTOR(mu);        // p population-level means (working scale)
  PARAMETER_VECTOR(log_sd);    // p population-level log SDs (diagonal between-individual D)
  PARAMETER_MATRIX(u);         // n x p individual-level deviations (random effects)

  int n = theta_hat.rows();
  int p = theta_hat.cols();

  Type nll = 0.0;

  vector<Type> sd = exp(log_sd);

  // Maximum Penalized Likelihood (MPL) penalty on the population-level SDs
  // (Chung et al. 2013): a log-gamma(alpha, 1) penalty on log_sd guarantees
  // a nondegenerate (non-boundary) posterior mode even when a parameter's
  // between-individual signal is weak.
  if (use_mpl == 1) {
    for (int j = 0; j < p; j++) {
      nll -= (mpl_alpha(j) - Type(1.0)) * log_sd(j);
    }
  }

  // Random-effects prior: u_ij ~ N(0, sd_j^2), iid across individuals.
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      nll -= dnorm(u(i, j), Type(0), sd(j), true);
    }
  }

  // Measurement-error model: theta_hat_i ~ MVN(mu + u_i, Shat_i), with
  // Shat_i the (fixed, known) Stage I Hessian-based covariance matrix.
  for (int i = 0; i < n; i++) {
    matrix<Type> Sigma_i(p, p);
    for (int a = 0; a < p; a++) {
      for (int b = 0; b < p; b++) {
        // Manually calculate 1D index to prevent Eigen from vectorizing TMB's 3D array lookup.
        // This bypasses the GCC 14.3 array-bounds false positive without disabling SIMD.
        int idx = a + (b * p) + (i * p * p);
        Sigma_i(a, b) = Shat[idx];
      }
    }
    vector<Type> resid(p);
    for (int j = 0; j < p; j++) {
      resid(j) = theta_hat(i, j) - (mu(j) + u(i, j));
    }
    MVNORM_t<Type> mvn(Sigma_i);
    nll += mvn(resid);
  }

  // Natural-scale population means and per-individual predictions
  // (delta-method SEs obtained automatically via ADREPORT).
  vector<Type> mu_nat(p);
  matrix<Type> ind_nat(n, p);
  for (int j = 0; j < p; j++) {
    if (col_type(j) == 1) {
      mu_nat(j) = exp(mu(j) + col_shift(j));
    } else {
      mu_nat(j) = mu(j);
    }
    for (int i = 0; i < n; i++) {
      Type ind_working = mu(j) + u(i, j);
      if (col_type(j) == 1) {
        ind_nat(i, j) = exp(ind_working + col_shift(j));
      } else {
        ind_nat(i, j) = ind_working;
      }
    }
  }

  REPORT(mu);
  REPORT(sd);
  REPORT(u);
  REPORT(mu_nat);
  REPORT(ind_nat);
  ADREPORT(mu);
  ADREPORT(sd);
  ADREPORT(mu_nat);
  ADREPORT(ind_nat);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
