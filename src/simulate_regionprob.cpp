#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat simulate_regionprob_cpp(int nSims, int n_cells, int n_layers, int n_covs,
                                  const arma::mat& beta_draws, const Rcpp::List& cov_mats_list,
                                  const arma::mat& mask_mat, bool show_progress) {

  // Safely extract covariates into a vector of armadillo matrices
  std::vector<arma::mat> cov_mats(n_covs);
  for(int j = 0; j < n_covs; j++) {
    SEXP obj = cov_mats_list[j];
    if(Rf_isMatrix(obj)) {
      cov_mats[j] = Rcpp::as<arma::mat>(obj);
    } else {
      arma::vec v = Rcpp::as<arma::vec>(obj);
      cov_mats[j] = arma::mat(v); // Deep copy to prevent memory decay
    }
  }

  // Matrix to store the regional probability for each draw and layer
  arma::mat P_sims(nSims, n_layers, arma::fill::zeros);
  arma::mat W(n_cells, n_layers);

  int last_pct = -1; // Tracker to limit console I/O

  for(int i = 0; i < nSims; i++) {

    // Throttle interrupt check
    if (i % 100 == 0 || i == nSims - 1) {
      Rcpp::checkUserInterrupt();
    }

    W.zeros();
    for(int j = 0; j < n_covs; j++) {
      double b = beta_draws(i, j);
      // Handle dynamic layers broadcasting
      if(cov_mats[j].n_cols == 1 && n_layers > 1) {
        for(int k = 0; k < n_layers; k++) {
          W.col(k) += cov_mats[j].col(0) * b;
        }
      } else {
        W += cov_mats[j] * b;
      }
    }

    for(int k = 0; k < n_layers; k++) {
      arma::vec W_k = W.col(k);
      arma::uvec valid_idx = arma::find_finite(W_k);

      if(valid_idx.n_elem > 0) {
        // Log-sum-exp normalization
        double max_W = arma::max(W_k.elem(valid_idx));
        arma::vec exp_W = arma::exp(W_k.elem(valid_idx) - max_W);
        double sum_pi = arma::sum(exp_W);
        arma::vec pi_k_valid = exp_W / sum_pi;

        // Extract the correct mask layer
        arma::vec m_k;
        if(mask_mat.n_cols == 1) {
          m_k = mask_mat.col(0);
        } else {
          m_k = mask_mat.col(k);
        }

        // The regional probability is the dot product of the distribution and the mask
        P_sims(i, k) = arma::dot(pi_k_valid, m_k.elem(valid_idx));
      } else {
        P_sims(i, k) = arma::datum::nan;
      }
    }

    // Progress Bar
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

  return P_sims;
}
