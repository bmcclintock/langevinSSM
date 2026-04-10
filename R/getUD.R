#' Compute Utilization Distribution
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates. The covariates must be on the same spatial grid and have the same spatial extent.
#' @param fit A \code{fitLangevin} object return by \code{\link{fitLangevin}}.
#' @param beta Numeric vector of habitat selection coefficients for the spatial covariates. The order of the coefficients must match the order of the covariates in \code{spatialCovs}.
#' @param log Logical indicating whether or not to return the log of the utilization distribution. Default: \code{TRUE}.
#' @param nsims Integer. Number of draws from the covariance matrix to use for estimating uncertainty in the UD. If \code{nsims > 0}, the returned list of raster will include additional elements for the (probability-scale) standard error (\code{SE}) and the coefficient of variation (\code{CV}) for each cell. Default: \code{0}.
#' @param show_progress Logical. If \code{TRUE}, displays a progress bar and messages during simulation. Default: \code{TRUE}.
#' @return An object of class \code{udLangevin}, which is a list containing the (log) utilization distribution (\code{UD}) as a \code{\link[terra]{SpatRaster-class}} object. If \code{nsims > 0}, the list also includes the estimated standard error (\code{SE}) and coefficient of variation (\code{CV}).
#' @seealso \code{\link{plotRaster}} for plotting the utilization distribution and covariates.
#' @examples
#' # exampleCovs included in package; see ?exampleCovs for details
#' UD <- getUD(exampleCovs, beta = c(-4, 6, 5, -0.1) )
#' @importFrom terra global nlyr varnames app
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
getUD <- function(spatialCovs, fit, beta, log = TRUE, nsims = 0, show_progress = TRUE) {
  if((missing(fit) & missing(beta)) | (!missing(fit) & !missing(beta))) stop("Either 'fit' or 'beta' must be provided")
  if(!missing(fit)) verify_signatures(fit, spatialCovs = spatialCovs)

  if(missing(beta)) {
    rn <- rownames(fit$estimates$natural)
    beta_idx <- which(grepl("^beta", rn))
    beta <- fit$estimates$natural[beta_idx, "Estimate"]
  }

  if(length(spatialCovs) != length(beta)) stop("length(spatialCovs) must equal length(beta)")

  # --- Internal Helper: Compute a normalized probability UD ---
  calc_prob_ud_base <- function(b_vec) {
    ud_rast <- spatialCovs[[1]] * b_vec[1]
    if(length(spatialCovs) > 1) {
      for (j in 2:length(spatialCovs)) ud_rast <- ud_rast + (spatialCovs[[j]] * b_vec[j])
    }
    # Log-Sum-Exp for stability
    max_log <- terra::global(ud_rast, "max", na.rm = TRUE)$max
    for(k in 1:terra::nlyr(ud_rast)) ud_rast[[k]] <- exp(ud_rast[[k]] - max_log[k])
    # Normalize to sum to 1
    layer_sums <- terra::global(ud_rast, "sum", na.rm = TRUE)$sum
    for(k in 1:terra::nlyr(ud_rast)) ud_rast[[k]] <- ud_rast[[k]] / layer_sums[k]
    return(ud_rast)
  }

  if (log) {
    ud_base <- spatialCovs[[1]] * beta[1]
    if(length(spatialCovs) > 1) {
      for (j in 2:length(spatialCovs)) ud_base <- ud_base + (spatialCovs[[j]] * beta[j])
    }
    base_name <- "log_UD"
  } else {
    ud_base <- calc_prob_ud_base(beta)
    base_name <- "UD"
  }
  names(ud_base) <- rep(base_name, terra::nlyr(ud_base))

  if (nsims == 0) return(class_udLangevin(list(UD = ud_base)))

  if (is.null(fit$estimates$cov_natural)) stop("fit$estimates$cov_natural not found. Refit model to get cov_natural.")
  beta_cov <- fit$estimates$cov_natural[beta_idx, beta_idx, drop = FALSE]
  beta_draws <- MASS::mvrnorm(nsims, beta, beta_cov)

  if (show_progress) message("   Simulating ", nsims, " draws to estimate UD uncertainty...")

  n_cells <- terra::ncell(spatialCovs[[1]])
  n_covs <- length(spatialCovs)
  n_ud_layers <- max(sapply(spatialCovs, terra::nlyr))

  cov_mats <- lapply(spatialCovs, function(x) {
    m <- terra::as.matrix(x, wide = FALSE)
    if(ncol(m) > 1) return(m)
    return(as.vector(m)) # Return as vector if single layer
  })

  sum_pi <- matrix(0, nrow = n_cells, ncol = n_ud_layers)
  sum_sq_pi <- matrix(0, nrow = n_cells, ncol = n_ud_layers)

  if (show_progress) pb <- utils::txtProgressBar(min = 0, max = nsims, style = 3)

  for (i in 1:nsims) {
    W <- matrix(0, nrow = n_cells, ncol = n_ud_layers)

    for (j in 1:n_covs) {
      b <- beta_draws[i, j]
      m <- cov_mats[[j]]
      W <- W + (m * b)
    }

    # Probability Normalization for each layer (Log-Sum-Exp)
    pi_sim <- W
    for (k in 1:n_ud_layers) {
      max_W <- max(W[, k], na.rm = TRUE)
      pi_sim[, k] <- exp(W[, k] - max_W)
      pi_sim[, k] <- pi_sim[, k] / sum(pi_sim[, k], na.rm = TRUE)
    }

    sum_pi <- sum_pi + pi_sim
    sum_sq_pi <- sum_sq_pi + (pi_sim^2)

    if (show_progress) utils::setTxtProgressBar(pb, i)
  }

  if (show_progress) close(pb) # Close connection to finalize the bar in console

  # Variance, SE, and CV calculation (standard one-pass formula)
  mean_pi_mat <- sum_pi / nsims
  var_pi_mat <- (sum_sq_pi - (sum_pi^2) / nsims) / (nsims - 1)
  se_pi_mat <- sqrt(pmax(var_pi_mat, 0))
  cv_pi_mat <- se_pi_mat / mean_pi_mat

  # Reconstruct Rasters
  template <- ud_base
  ud_se <- terra::setValues(template, se_pi_mat)
  ud_cv <- terra::setValues(template, cv_pi_mat)

  names(ud_se) <- rep("UD_SE", n_ud_layers)
  names(ud_cv) <- rep("UD_CV", n_ud_layers)

  # Ensure time attributes are preserved
  dyn_idx <- which(sapply(spatialCovs, terra::nlyr) == n_ud_layers)[1]
  if (!is.na(dyn_idx) && !is.null(terra::time(spatialCovs[[dyn_idx]]))) {
    terra::time(ud_se) <- terra::time(spatialCovs[[dyn_idx]])
    terra::time(ud_cv) <- terra::time(spatialCovs[[dyn_idx]])
  }

  return(class_udLangevin(list(UD = ud_base, SE = ud_se, CV = ud_cv)))
}
