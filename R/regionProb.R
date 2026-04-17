#' Calculate regional probability from fitted Langevin model
#'
#' This function calculates the probability of an animal being in a specified region (defined by a mask) based on the fitted Langevin model. It provides both a point estimate and uncertainty quantification using the Delta method and Monte Carlo simulations. The Delta method is a quick approximation, while the Monte Carlo approach captures the full uncertainty from the covariance matrix. If dynamic (time-varying) covariates are provided, it calculates the regional probability for each temporal layer.
#'
#' @param fit A \code{fitLangevin} object to supply the covariance matrix.
#' @param spatialCovs List of \code{\link[terra]{SpatRaster}} spatial covariates used to fit the model.
#' @param mask A \code{\link[terra]{SpatRaster}} with 1s in the region of interest and NAs/0s elsewhere.
#' @param nSims Integer. Number of simulations to generate credible intervals. Default: \code{0}.
#' @param level Numeric. The confidence level for the intervals. Default: \code{0.95}.
#' @param show_progress Logical. If \code{TRUE}, displays a progress bar and messages. Default: \code{TRUE}.
#' @return A \code{regLangevin} object (which is a list) containing the point estimate(s), Delta method SE, Monte Carlo SE/CIs, and the underlying spatial rasters (\code{prob_raster} and \code{mask}) used for plotting.
#'
#' @examples
#' \donttest{
#' # fit the underdamped Langevin model
#' fit <- fitLangevin(data = exampleDat,
#'                    spatialCovs = exampleCovs,
#'                    silent = TRUE)
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
# #' @importFrom MASS mvrnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats quantile sd qnorm
#' @importFrom terra as.matrix ncell nlyr time
#' @export
regionProb <- function(fit, spatialCovs, mask, nSims = 0, level = 0.95, show_progress = TRUE) {

  if(!inherits(fit, "fitLangevin")) stop("'fit' must be a fitLangevin object")
  if(!is.list(spatialCovs) || !all(sapply(spatialCovs, inherits, "SpatRaster"))) stop("'spatialCovs' must be a list of SpatRaster objects")
  if(!inherits(mask, "SpatRaster")) stop("'mask' must be a SpatRaster object")

  verify_signatures(fit, spatialCovs = spatialCovs)

  if (!terra::compareGeom(spatialCovs[[1]], mask, stopOnError = FALSE)) {
    stop("The 'mask' raster must share the same projection (CRS), extent, and resolution as the rasters in 'spatialCovs'.")
  }
  if(length(nSims) != 1 || !is.numeric(nSims) || nSims < 0) stop("'nSims' must be a single non-negative integer")
  if(length(level) != 1 || !is.numeric(level) || level <= 0 || level >= 1) stop("'level' must be a single numeric value between 0 and 1")

  rn <- rownames(fit$estimates$natural)
  beta_idx <- which(grepl("^beta", rn))
  beta <- fit$estimates$natural[beta_idx, "Estimate"]
  beta_cov <- fit$covariance$natural[beta_idx, beta_idx, drop = FALSE]

  n_cells <- terra::ncell(spatialCovs[[1]])
  n_covs <- length(spatialCovs)
  n_layers <- max(sapply(spatialCovs, terra::nlyr))

  mask_mat <- terra::as.matrix(mask, wide = FALSE)
  cov_mats <- lapply(spatialCovs, function(x) terra::as.matrix(x, wide = FALSE))

  P_est <- numeric(n_layers)
  SE_delta <- numeric(n_layers)
  CI_delta <- matrix(NA, nrow = n_layers, ncol = 2)
  a <- (1 - level) / 2
  z_val <- stats::qnorm(1 - a)

  SE_sim <- if(nSims > 0) numeric(n_layers) else NULL
  CI_sim <- if(nSims > 0) matrix(NA, nrow = n_layers, ncol = 2) else NULL
  simulated_draws <- if(nSims > 0) matrix(NA, nrow = nSims, ncol = n_layers) else NULL

  pi_mat <- matrix(NA, nrow = n_cells, ncol = n_layers)

  if (nSims > 0) {
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' needed for simulation. Please install it.", call. = FALSE)
    message("Simulating ", nSims, " draws for regional probability...")
    beta_draws <- MASS::mvrnorm(nSims, beta, beta_cov)
  }

  total_steps <- if(nSims > 0) n_layers * nSims else n_layers
  if (show_progress) pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)
  counter <- 0

  for (k in 1:n_layers) {
    mk <- if (ncol(mask_mat) == 1) mask_mat[, 1] else mask_mat[, k]
    mk[is.na(mk)] <- 0

    Ck <- matrix(NA, nrow = n_cells, ncol = n_covs)
    for (j in 1:n_covs) {
      Ck[, j] <- if (ncol(cov_mats[[j]]) == 1) cov_mats[[j]][, 1] else cov_mats[[j]][, k]
    }

    W <- as.numeric(Ck %*% beta)
    pi_vec <- exp(W - max(W, na.rm = TRUE))
    pi_vec <- pi_vec / sum(pi_vec, na.rm = TRUE)
    pi_mat[, k] <- pi_vec

    P_est[k] <- sum(pi_vec * mk, na.rm = TRUE)

    # Delta Method
    mu_C <- colSums(Ck * pi_vec, na.rm = TRUE)
    mu_C_mask <- colSums(Ck * (pi_vec * mk), na.rm = TRUE)
    g_P <- mu_C_mask - (P_est[k] * mu_C)
    var_P <- as.numeric(t(g_P) %*% beta_cov %*% g_P)

    SE_delta[k] <- sqrt(pmax(var_P, 0))
    CI_delta[k, ] <- c(max(0, P_est[k] - z_val * SE_delta[k]), min(1, P_est[k] + z_val * SE_delta[k]))

    # Monte Carlo Simulations
    if (nSims > 0) {
      P_sims <- numeric(nSims)
      for (i in 1:nSims) {
        W_sim <- as.numeric(Ck %*% beta_draws[i, ])
        pi_sim <- exp(W_sim - max(W_sim, na.rm = TRUE))
        pi_sim <- pi_sim / sum(pi_sim, na.rm = TRUE)
        P_sims[i] <- sum(pi_sim * mk, na.rm = TRUE)
        counter <- counter + 1
        if (show_progress) utils::setTxtProgressBar(pb, counter)
      }
      SE_sim[k] <- stats::sd(P_sims)
      CI_sim[k, ] <- as.numeric(stats::quantile(P_sims, probs = c(a, 1 - a)))
      simulated_draws[, k] <- P_sims
    } else {
      counter <- counter + 1
      if (show_progress) utils::setTxtProgressBar(pb, counter)
    }
  }
  if (show_progress) close(pb)

  # Prepare multi-layer raster output
  dyn_idx <- which(sapply(spatialCovs, terra::nlyr) == n_layers)[1]
  prob_rast <- if(!is.na(dyn_idx)) spatialCovs[[dyn_idx]] else spatialCovs[[1]]
  terra::values(prob_rast) <- pi_mat
  names(prob_rast) <- rep("Probability", n_layers)

  out <- list(Point_Estimate = P_est, SE_delta = SE_delta, CI_delta = CI_delta,
              prob_raster = prob_rast, mask = mask, level = level)

  if (nSims > 0) {
    out$SE_sim <- SE_sim
    out$CI_sim <- CI_sim
    out$simulated_draws <- simulated_draws
  }

  class(out) <- c("regLangevin", "list")
  return(out)
}
