#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List simulate_ud_cpp(int nSims, int n_cells, int n_ud_layers, int n_covs,
                           const arma::mat& beta_draws, const Rcpp::List& cov_mats_list,
                           bool show_progress) {

  // Extract covariates into an std::vector of arma::mat
  std::vector<arma::mat> cov_mats(n_covs);
  for(int j = 0; j < n_covs; j++) {
    SEXP obj = cov_mats_list[j];
    if(Rf_isMatrix(obj)) {
      cov_mats[j] = Rcpp::as<arma::mat>(obj);
    } else {
      arma::vec v = Rcpp::as<arma::vec>(obj);
      cov_mats[j] = arma::mat(v);
    }
  }

  arma::mat mean_pi(n_cells, n_ud_layers, arma::fill::zeros);
  arma::mat M2_pi(n_cells, n_ud_layers, arma::fill::zeros);
  arma::mat W(n_cells, n_ud_layers);
  arma::mat pi_sim(n_cells, n_ud_layers);

  int last_pct = -1;

  for(int i = 0; i < nSims; i++) {

    if (i % 100 == 0 || i == nSims - 1) {
      Rcpp::checkUserInterrupt();
    }

    W.zeros();
    for(int j = 0; j < n_covs; j++) {
      double b = beta_draws(i, j);
      if(cov_mats[j].n_cols == 1 && n_ud_layers > 1) {
        for(int k = 0; k < n_ud_layers; k++) {
          W.col(k) += cov_mats[j].col(0) * b;
        }
      } else {
        W += cov_mats[j] * b;
      }
    }

    for(int k = 0; k < n_ud_layers; k++) {
      arma::vec W_k = W.col(k);
      arma::uvec valid_idx = arma::find_finite(W_k);

      arma::vec pi_k(n_cells);
      pi_k.fill(arma::datum::nan);

      if(valid_idx.n_elem > 0) {
        double max_W = arma::max(W_k.elem(valid_idx));
        arma::vec exp_W = arma::exp(W_k.elem(valid_idx) - max_W);
        double sum_pi = arma::sum(exp_W);
        pi_k.elem(valid_idx) = exp_W / sum_pi;
      }
      pi_sim.col(k) = pi_k;
    }

    // Welford's algorithm
    arma::mat delta = pi_sim - mean_pi;
    mean_pi += delta / (i + 1.0);
    arma::mat delta2 = pi_sim - mean_pi;
    M2_pi += delta % delta2;

    // Mimic utils::txtProgressBar(style = 3)
    if (show_progress) {
      int pct = (i + 1) * 100 / nSims;

      if (pct > last_pct) {
        std::string bar(pct / 2, '=');
        std::string space(50 - pct / 2, ' ');
        Rprintf("\r  |%s%s| %3d%%", bar.c_str(), space.c_str(), pct);
        last_pct = pct;
      }
    }
  }

  if (show_progress) Rprintf("\n");

  return Rcpp::List::create(
    Rcpp::Named("mean_pi") = mean_pi,
    Rcpp::Named("M2_pi") = M2_pi
  );
}
