#' Calculate regional probability from fitted Langevin model
#'
#' This function calculates the probability of an animal being in a specified region (defined by a mask) based on the fitted Langevin model. It provides both a point estimate and uncertainty quantification using the Delta method and Monte Carlo simulations. The Delta method is a quick approximation, while the Monte Carlo approach captures the full uncertainty from the covariance matrix.
#'
#' @param fit A \code{fitLangevin} object to supply the covariance matrix.
#' @param spatialCovs List of \code{\link[terra]{SpatRaster}} spatial covariates used to fit the model.
#' @param mask A \code{\link[terra]{SpatRaster}} with 1s in the region of interest and NAs/0s elsewhere.
#' @param nSims Integer. Number of simulations to generate credible intervals.
#' @param show_progress Logical. If \code{TRUE}, displays a progress bar and messages. Default: \code{TRUE}.
#' @return A list containing the point estimate, Delta method SE, and Monte Carlo SE/CIs.
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
#'                        nSims = 1000)
#'
#' # point estimate and 95% Monte Carlo credible interval
#' reg_prob$Point_Estimate
#' reg_prob$CI_sim_95
#' }
#' @seealso \code{\link{getUD}} for calculating the utilization distribution.
# #' @importFrom MASS mvrnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats quantile sd
#' @export
regionProb <- function(fit, spatialCovs, mask, nSims = 1000, show_progress = TRUE) {

  if(!inherits(fit, "fitLangevin")) stop("'fit' must be a fitLangevin object")
  if(!is.list(spatialCovs) || !all(sapply(spatialCovs, inherits, "SpatRaster"))) stop("'spatialCovs' must be a list of SpatRaster objects")
  if(!inherits(mask, "SpatRaster")) stop("'mask' must be a SpatRaster object")

  verify_signatures(fit, spatialCovs = spatialCovs)

  if (!terra::compareGeom(spatialCovs[[1]], mask, stopOnError = FALSE)) {
    stop("The 'mask' raster must share the same projection (CRS), extent, and resolution as the rasters in 'spatialCovs'.")
  }
  if(length(nSims) != 1 || !is.numeric(nSims) || nSims < 0) stop("'nSims' must be a single non-negative integer")

  # 1. Extract beta estimates and covariance
  rn <- rownames(fit$estimates$natural)
  beta_idx <- which(grepl("^beta", rn))
  beta <- fit$estimates$natural[beta_idx, "Estimate"]
  beta_cov <- fit$covariance$natural[beta_idx, beta_idx, drop = FALSE]

  # 2. Prepare data matrices (using the first layer if dynamic)
  mask_vec <- terra::as.matrix(mask, wide = FALSE)[, 1]
  mask_vec[is.na(mask_vec)] <- 0

  C_mat <- sapply(spatialCovs, function(x) terra::as.matrix(x, wide = FALSE)[, 1])

  # 3. Calculate Base UD - Flatten to numeric to avoid conformable array errors
  W <- as.numeric(C_mat %*% beta)
  pi_vec <- exp(W - max(W, na.rm = TRUE))
  pi_vec <- pi_vec / sum(pi_vec, na.rm = TRUE)

  # Point Estimate (sum of probabilities in the mask)
  P_est <- sum(pi_vec * mask_vec, na.rm = TRUE)
  out <- list(Point_Estimate = P_est)

  # ==========================================
  # 4. Delta Method Approximation for the Sum
  # ==========================================
  # The gradient of the sum is the sum of the gradients!
  mu_C <- colSums(C_mat * pi_vec, na.rm = TRUE)
  mu_C_mask <- colSums(C_mat * (pi_vec * mask_vec), na.rm = TRUE)

  g_P <- mu_C_mask - (P_est * mu_C)
  var_P <- as.numeric(t(g_P) %*% beta_cov %*% g_P)

  out$SE_delta <- sqrt(pmax(var_P, 0))
  # Standard Wald 95% Confidence Interval (capped at 0 and 1)
  out$CI_delta_95 <- c(max(0, P_est - 1.96 * out$SE_delta),
                       min(1, P_est + 1.96 * out$SE_delta))

  # ==========================================
  # 5. Monte Carlo Simulations (Highly Recommended)
  # ==========================================
  if (nSims > 0) {
    message("Simulating ", nSims, " draws for regional probability...")

    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package \"MASS\" needed for simulation Please install it.", call. = FALSE)

    beta_draws <- MASS::mvrnorm(nSims, beta, beta_cov)
    P_sims <- numeric(nSims)

    if (show_progress) pb <- utils::txtProgressBar(min = 0, max = nSims, style = 3)

    for (i in 1:nSims) {
      # Flatten to numeric
      W_sim <- as.numeric(C_mat %*% beta_draws[i, ])
      pi_sim <- exp(W_sim - max(W_sim, na.rm = TRUE))
      pi_sim <- pi_sim / sum(pi_sim, na.rm = TRUE)

      P_sims[i] <- sum(pi_sim * mask_vec, na.rm = TRUE)
      if (show_progress) utils::setTxtProgressBar(pb, i)
    }
    if (show_progress) close(pb)

    out$SE_sim <- stats::sd(P_sims)
    out$CI_sim_95 <- as.numeric(stats::quantile(P_sims, probs = c(0.025, 0.975)))
    out$simulated_draws <- P_sims # You can plot a histogram of this!
  }

  return(out)
}
