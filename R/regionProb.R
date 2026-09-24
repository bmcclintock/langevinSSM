#' Calculate regional probability from fitted Langevin model
#'
#' This function calculates the probability of an animal being in a specified region (defined by a mask) based on the fitted Langevin model. It provides both a point estimate and uncertainty quantification using the Delta method and Monte Carlo simulations. The Delta method is a quick approximation, while the Monte Carlo approach captures the full uncertainty from the covariance matrix. If dynamic (time-varying) covariates are provided, it calculates the regional probability for each temporal layer.
#'
#' @details
#' \strong{Calculations for \code{fitLangevin} objects:}
#' \itemize{
#'   \item \strong{Point Estimate:} Calculated directly from the point estimates of the habitat selection coefficients.
#'   \item \strong{Standard Errors (Delta Method):} The analytical gradient of the masked spatial softmax is calculated and multiplied by the Hessian-based covariance matrix.
#'   \item \strong{Standard Errors (Monte Carlo):} If \code{nSims > 0}, draws are generated from the multivariate normal distribution of the coefficients \eqn{MVN(\hat{\beta}, \hat{\Sigma}_\beta)}, and the regional probability is computed for each draw to form a credible interval.
#' }
#'
#' \strong{Calculations for \code{hierLangevin} objects:}
#' \itemize{
#'   \item \strong{Population Probability:} Because spatial selection is non-linear, the regional probability is integrated over the estimated between-individual variance \eqn{\mathrm{N}(\mu_\beta, \sigma^2_\beta)} using Monte Carlo integration (with \code{nSims} draws).
#'   \item \strong{Population SE:} Estimated using a nested Monte Carlo approach. Outer draws capture uncertainty in the hyperparameters \eqn{(\mu_\beta, \log \sigma_\beta)} from the Stage II working covariance matrix. For each outer draw, an inner Monte Carlo integration computes the marginal regional probability. The SE is the standard deviation across the outer draws.
#'   \item \strong{Individual Probability (BLUP):} Computed using the specific individual's Best Linear Unbiased Predictor (BLUP) for the selection coefficients: \eqn{\beta_i = \mu_\beta + u_i}.
#'   \item \strong{Individual SE:} Draws are generated directly from the sparse joint precision matrix of the fixed effects and random effects to properly propagate the conditional uncertainty in the BLUPs.
#' }
#'
#' @param fit A \code{fitLangevin} or \code{hierLangevin} object to supply the covariance matrix.
#' @param spatialCovs List of \code{\link[terra]{SpatRaster}} spatial covariates used to fit the model.
#' @param mask A \code{\link[terra]{SpatRaster}} with 1s in the region of interest and NAs/0s elsewhere.
#' @param nSims Integer. Number of simulations to generate credible intervals. For \code{hierLangevin} objects, this dictates the number of draws from the between-individual MVN distribution used to integrate the marginal population-level probability, and \strong{must} be > 0. Default: \code{0}.
#' @param level Numeric. The confidence level for the intervals. Default: \code{0.95}.
#' @param individual Optional vector of individual IDs (character or numeric) to calculate individual-level probabilities for hierarchical models. If \code{NULL} (default), only the population-level probability is calculated. Default: \code{NULL}.
#' @param show_progress Logical. If \code{TRUE}, displays a progress bar and messages. Default: \code{TRUE}.
#' @return A \code{regLangevin} object (which is a list) containing the point estimate(s), Delta method SE, Monte Carlo SE/CIs, and the underlying spatial rasters (\code{prob_raster} and \code{mask}) used for plotting. If multiple individuals are requested, a named list of \code{regLangevin} objects is returned.
#'
#' @examples
#' \donttest{
#' # fit the underdamped Langevin model
#' fit <- fitLangevin(data = exampleDat,
#'                    spatialCovs = exampleCovs)
#'
#' # create a spatial mask for the region of interest
#' d2c <- exampleCovs$d2c < 2.5
#'
#' # calculate the probability of the animal being in the region
#' reg_prob <- regionProb(fit = fit,
#'                        spatialCovs = exampleCovs,
#'                        mask = d2c,
#'                        nSims = 1000,
#'                        level = 0.95)
#'
#' # point estimate and 95% Monte Carlo credible interval
#' reg_prob
#'
#' # plot the regional probability
#' plot(reg_prob, log = FALSE)
#' }
#' @seealso \code{\link{getUD}} for calculating the utilization distribution.
#' @importFrom stats quantile sd qnorm
#' @importFrom terra as.matrix ncell nlyr time
#' @export
regionProb <- function(fit, spatialCovs, mask, nSims = 0, level = 0.95, individual = NULL, show_progress = TRUE) {

  if(!inherits(fit, "fitLangevin") && !inherits(fit, "hierLangevin")) {
    stop("'fit' must be a fitLangevin object or a hierLangevin object")
  }
  is_hierLangevin <- inherits(fit, "hierLangevin")

  if(!is.list(spatialCovs) || !all(sapply(spatialCovs, inherits, "SpatRaster"))) stop("'spatialCovs' must be a list of SpatRaster objects")
  if(!inherits(mask, "SpatRaster")) stop("'mask' must be a SpatRaster object")

  if (is_hierLangevin && nSims <= 0) {
    stop("For 'hierLangevin' models, population-level regional probabilities must be integrated over the between-individual variance using Monte Carlo methods. Please specify 'nSims > 0' (e.g., nSims = 1000).")
  }

  verify_signatures(fit, spatialCovs = spatialCovs)
  if(!is_hierLangevin) boundsWarning(fit)

  if (!terra::compareGeom(spatialCovs[[1]], mask, stopOnError = FALSE)) {
    stop("The 'mask' raster must share the same projection (CRS), extent, and resolution as the rasters in 'spatialCovs'.")
  }
  if(length(nSims) != 1 || !is.numeric(nSims) || nSims < 0) stop("'nSims' must be a single non-negative integer")
  if(length(level) != 1 || !is.numeric(level) || level <= 0 || level >= 1) stop("'level' must be a single numeric value between 0 and 1")

  barrier <- NULL
  lambda <- NULL
  scaleFactor <- 1

  if(!is.null(fit$conditions$barrier)) barrier <- fit$conditions$barrier
  if(!is.null(fit$conditions$lambda)) lambda <- fit$conditions$lambda
  if(!is.null(fit$conditions$scaleFactor)) scaleFactor <- fit$conditions$scaleFactor

  bar_info <- .prep_barrier_raster(spatialCovs, barrier, lambda, scaleFactor)
  spatialCovs <- bar_info$spatialCovs
  mod_spatialCovs <- bar_info$mod_spatialCovs

  beta_info <- .extract_langevin_beta_list(fit = fit, individual = individual)
  beta_list <- beta_info$beta_list
  hierarchical_logical <- beta_info$hierarchical_logical

  cov_info <- .extract_langevin_cov(fit, is_hierLangevin, hierarchical_logical)
  beta_cov <- cov_info$beta_cov
  beta_idx_cov <- cov_info$beta_idx_cov
  can_calc_se <- cov_info$can_calc_se

  if (can_calc_se && is_hierLangevin) {
    message("   Skipping Delta Method (mathematically invalid for hierLangevin marginal probabilities). Relying on Monte Carlo...")
  }

  if (!can_calc_se && nSims > 0) stop("Cannot estimate uncertainty (nSims > 0) without a valid covariance matrix.")
  run_sim_se <- can_calc_se && (nSims > 0)

  n_cells <- terra::ncell(spatialCovs[[1]])
  n_covs <- length(mod_spatialCovs)
  n_layers <- max(sapply(spatialCovs, terra::nlyr))

  mask_mat <- terra::as.matrix(mask, wide = FALSE)
  mask_mat[is.na(mask_mat)] <- 0

  cov_mats <- lapply(mod_spatialCovs, function(x) {
    m <- terra::as.matrix(x, wide = FALSE)
    if(ncol(m) > 1) return(m)
    return(as.vector(m))
  })

  a <- (1 - level) / 2
  z_val <- stats::qnorm(1 - a)

  draw_data <- if (run_sim_se && hierarchical_logical && length(beta_list) > 1) {
    .sample_joint_precision_re(fit, is_hierLangevin, beta_idx_cov, nSims)
  } else {
    list(has_Q = FALSE)
  }

  out_list <- list()

  for (i in seq_along(beta_list)) {
    nm <- names(beta_list)[i]
    item <- beta_list[[i]]

    P_est <- numeric(n_layers)
    SE_delta <- NULL
    CI_delta <- NULL
    SE_sim <- NULL
    CI_sim <- NULL
    simulated_draws <- NULL

    pi_mat <- matrix(NA, nrow = n_cells, ncol = n_layers)

    if (is.list(item) && !is.null(item$type) && item$type == "hierLangevin_population") {
      n_mc <- nSims
      message("   Using Monte Carlo integration (", n_mc, " draws) for Population region probability...")

      pop_cov <- diag(item$sd^2, nrow = length(item$sd))
      beta_draws_marg <- as.matrix(MASS::mvrnorm(n_mc, item$mu, pop_cov))
      if (!is.null(barrier)) beta_draws_marg <- cbind(beta_draws_marg, 1)

      cpp_res_marg <- simulate_ud_cpp(
        nSims = n_mc, n_cells = n_cells, n_ud_layers = n_layers,
        n_covs = n_covs, beta_draws = beta_draws_marg,
        cov_mats_list = cov_mats, show_progress = FALSE
      )
      pi_mat <- cpp_res_marg$mean_pi

      for (k in 1:n_layers) {
        mk <- if (ncol(mask_mat) == 1) mask_mat[, 1] else mask_mat[, k]
        P_est[k] <- sum(pi_mat[, k] * mk, na.rm = TRUE)
      }

      if (run_sim_se) {
        message("     Running nested Monte Carlo for Population regional uncertainty...")
        M_inner <- 30
        pop_params_draws <- MASS::mvrnorm(nSims, c(item$mu, log(item$sd)), beta_cov)

        P_sims_mat <- matrix(NA, nrow = nSims, ncol = n_layers)

        if (show_progress) pb <- txtProgressBar(min = 0, max = nSims, style = 3, width=50)

        for (s in 1:nSims) {
          mu_s <- pop_params_draws[s, 1:length(item$mu)]
          sd_s <- exp(pop_params_draws[s, (length(item$mu)+1):ncol(pop_params_draws)])

          inner_cov <- diag(sd_s^2, nrow = length(sd_s))
          inner_beta_draws <- as.matrix(MASS::mvrnorm(M_inner, mu_s, inner_cov))
          if (!is.null(barrier)) inner_beta_draws <- cbind(inner_beta_draws, 1)

          inner_res <- simulate_regionprob_cpp(
            nSims = M_inner, n_cells = n_cells, n_layers = n_layers,
            n_covs = n_covs, beta_draws = inner_beta_draws,
            cov_mats_list = cov_mats, mask_mat = mask_mat,
            show_progress = FALSE
          )

          P_sims_mat[s, ] <- colMeans(inner_res, na.rm = TRUE)
          if (show_progress) setTxtProgressBar(pb, s)
        }
        if (show_progress) close(pb)

        SE_sim <- apply(P_sims_mat, 2, stats::sd, na.rm = TRUE)
        if (n_layers == 1) {
          CI_sim <- matrix(stats::quantile(P_sims_mat[, 1], probs = c(a, 1 - a), na.rm = TRUE), nrow = 1)
        } else {
          CI_sim <- t(apply(P_sims_mat, 2, stats::quantile, probs = c(a, 1 - a), na.rm = TRUE))
        }
        simulated_draws <- P_sims_mat
      }

    } else {
      b_vec <- if (is.list(item)) item$vec else item
      if (!is.null(barrier)) b_vec <- c(b_vec, 1)

      if (can_calc_se && !hierarchical_logical) {
        SE_delta <- numeric(n_layers)
        CI_delta <- matrix(NA, nrow = n_layers, ncol = 2)
      }

      for (k in 1:n_layers) {
        mk <- if (ncol(mask_mat) == 1) mask_mat[, 1] else mask_mat[, k]
        Ck <- matrix(NA, nrow = n_cells, ncol = n_covs)
        for (j in 1:n_covs) {
          Ck[, j] <- if (is.matrix(cov_mats[[j]]) && ncol(cov_mats[[j]]) > 1) cov_mats[[j]][, k] else cov_mats[[j]]
        }

        W <- as.numeric(Ck %*% b_vec)
        pi_vec <- exp(W - max(W, na.rm = TRUE))
        pi_vec <- pi_vec / sum(pi_vec, na.rm = TRUE)
        pi_mat[, k] <- pi_vec

        P_est[k] <- sum(pi_vec * mk, na.rm = TRUE)

        # Delta Method (only if not hierLangevin Population marginals)
        if (can_calc_se && !hierarchical_logical) {
          mu_C <- colSums(Ck * pi_vec, na.rm = TRUE)
          mu_C_mask <- colSums(Ck * (pi_vec * mk), na.rm = TRUE)

          if (!is.null(barrier)) {
            g_P <- mu_C_mask[1:(n_covs-1)] - (P_est[k] * mu_C[1:(n_covs-1)])
          } else {
            g_P <- mu_C_mask - (P_est[k] * mu_C)
          }
          var_P <- as.numeric(t(g_P) %*% beta_cov %*% g_P)
          SE_delta[k] <- sqrt(pmax(var_P, 0))
          CI_delta[k, ] <- c(max(0, P_est[k] - z_val * SE_delta[k]), min(1, P_est[k] + z_val * SE_delta[k]))
        }
      }

      if (run_sim_se) {
        if (hierarchical_logical) {
          if (nm != "Population" && draw_data$has_Q) {
            message("   Simulating region probability uncertainty for individual: ", sub("^ID_", "", nm), "...")
            beta_draws <- .get_individual_beta_draws(nm, fit, is_hierLangevin, draw_data, beta_idx_cov, nSims)
          } else {
            beta_draws <- NULL
          }
        } else {
          message("   Simulating", ifelse(show_progress," ",paste0(" ",nSims," ")), "draws for regional probability...")
          beta_draws <- MASS::mvrnorm(nSims, b_vec[1:nrow(beta_cov)], beta_cov)
        }

        if (!is.null(beta_draws)) {
          if (!is.null(barrier)) beta_draws <- cbind(beta_draws, 1)

          P_sims_mat <- simulate_regionprob_cpp(
            nSims = nSims, n_cells = n_cells, n_layers = n_layers,
            n_covs = n_covs, beta_draws = beta_draws,
            cov_mats_list = cov_mats, mask_mat = mask_mat,
            show_progress = show_progress
          )

          SE_sim <- apply(P_sims_mat, 2, stats::sd, na.rm = TRUE)
          if (n_layers == 1) {
            CI_sim <- matrix(stats::quantile(P_sims_mat[, 1], probs = c(a, 1 - a), na.rm = TRUE), nrow = 1)
          } else {
            CI_sim <- t(apply(P_sims_mat, 2, stats::quantile, probs = c(a, 1 - a), na.rm = TRUE))
          }
          simulated_draws <- P_sims_mat
        }
      }
    }

    dyn_idx <- which(sapply(spatialCovs, terra::nlyr) == n_layers)[1]
    prob_rast <- if(!is.na(dyn_idx)) spatialCovs[[dyn_idx]] else spatialCovs[[1]]
    terra::values(prob_rast) <- pi_mat
    names(prob_rast) <- rep(ifelse(nm == "UD", "Probability", paste0("Probability_", nm)), n_layers)

    res_obj <- list(Point_Estimate = P_est, prob_raster = prob_rast, mask = mask, level = level)
    if (!is.null(SE_delta)) {
      res_obj$SE_delta <- SE_delta
      res_obj$CI_delta <- CI_delta
    }
    if (nSims > 0 && !is.null(SE_sim)) {
      res_obj$SE_sim <- SE_sim
      res_obj$CI_sim <- CI_sim
      res_obj$simulated_draws <- simulated_draws
    }

    class(res_obj) <- c("regLangevin", "list")
    out_list[[nm]] <- res_obj
  }

  if (length(out_list) == 1) return(out_list[[1]])
  class(out_list) <- c("regLangevin", "list")
  return(out_list)
}
