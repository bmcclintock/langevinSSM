#' Compute Utilization Distribution
#'
#' This function computes the utilization distribution (UD) for a given set of spatial covariates and habitat selection coefficients. It also estimates uncertainty in the UD using the Delta method and Monte Carlo simulations (if \code{nSims>0}). The resulting UD and uncertainty estimates are returned as \code{\link[terra]{SpatRaster-class}} objects, which can be plotted and analyzed further.
#'
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates. The covariates must be on the same spatial grid and have the same spatial extent.
#' @param fit A \code{fitLangevin} object return by \code{\link{fitLangevin}}.
#' @param beta Numeric vector of habitat selection coefficients for the spatial covariates. The order of the coefficients must match the order of the covariates in \code{spatialCovs}.
#' @param barrier Character string. The name of the barrier in \code{spatialCovs} that is represented as a signed distance field (see \code{\link{prepBarrier}}). If provided, this raster is exclusively used for the barrier penalty and is not included in the habitat selection covariates. Default: \code{NULL} (no barrier). If \code{fit} is provided, this is extracted automatically.
#' @param lambda Numeric. The penalty weight for the barrier constraint. Default: \code{NULL}. If \code{fit} is provided, this is extracted automatically.
#' @param scaleFactor Internal scaling factor for the coordinates and parameters (see \code{\link{fitLangevin}}). Automatically extracted if \code{fit} is provided. Default: 1 (no scaling).
#' @param log Logical indicating whether or not to return the log of the utilization distribution. Default: \code{TRUE}.
#' @param nSims Integer. Number of draws from the covariance matrix to use for estimating Monte Carlo uncertainty in the UD. If \code{nSims > 0}, the returned raster stack will include additional layers for simulated SE and CV. Default: \code{0}.
#' @param show_progress Logical. If \code{TRUE}, displays a progress bar for simulations. Default: \code{TRUE}.
#' @param plot Logical. Plot the resulting UD using \code{\link{plotUD}}? Default: \code{TRUE}.
#' @param maskRast \code{\link[terra]{SpatRaster-class}} object for areas to be masked out (set to \code{NA}) before plotting the UD. Default: \code{NULL} (no mask).
#' @return A \code{\link[terra]{SpatRaster}} object. It contains the (log) utilization distribution. It will also contain layers for Delta method standard errors (\code{UD_SE_delta}) and CVs (\code{UD_CV_delta}) if the covariance matrix is available. If \code{nSims > 0}, it adds simulated layers (\code{UD_SE_sim}, \code{UD_CV_sim}).
#' @seealso \code{\link{plotUD}}, \code{\link{regionProb}}.
#' @examples
#' # exampleCovs included in package; see ?exampleCovs for details
#' UD <- getUD(exampleCovs, beta = c(-4, 6, 5, -0.1) )
#' @importFrom terra global nlyr varnames app setValues compareGeom mask
#' @importFrom stats setNames
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
getUD <- function(spatialCovs, fit, beta, barrier = NULL, lambda = NULL, scaleFactor = 1, log = TRUE, nSims = 0, show_progress = TRUE, plot = TRUE, maskRast = NULL) {
  if((missing(fit) & missing(beta)) | (!missing(fit) & !missing(beta))) stop("Either 'fit' or 'beta' must be provided, but not both.")
  if(!missing(fit)) verify_signatures(fit, spatialCovs = spatialCovs)

  if(!missing(fit)) {
    barrier <- fit$conditions$barrier
    lambda <- fit$conditions$lambda
    if (!is.null(fit$conditions$scaleFactor)) scaleFactor <- fit$conditions$scaleFactor
  }

  if (!is.null(barrier)) {
    if (!(barrier %in% names(spatialCovs))) stop(sprintf("Barrier raster '%s' not found in spatialCovs.", barrier))

    barrier_sdf <- spatialCovs[[barrier]]
    if(!isTRUE(attr(barrier_sdf,"barLangevin"))) stop("barrier is not a 'barLangevin' object created by prepBarrier.")
    spatialCovs[[barrier]] <- NULL

    if (is.null(lambda)) stop("To plot a barrier without a fitted model ('fit'), you must manually specify 'lambda'.")
    .validate_lambda(lambda)
  }

  if(missing(beta)) {
    rn <- rownames(fit$estimates$natural)
    beta_idx <- which(grepl("^beta", rn))
    beta <- fit$estimates$natural[beta_idx, "Estimate"]
  }

  if(length(spatialCovs) != length(beta)) stop("length(spatialCovs) must equal length(beta). Note that if a barrier is specified, the barrier does not receive a beta coefficient.")

  if(!missing(fit)){
    boundsWarning(fit)
  }

  can_calc_se <- !missing(fit) && !is.null(fit$covariance$natural)
  if (can_calc_se) {
    beta_cov <- fit$covariance$natural[beta_idx, beta_idx, drop = FALSE]
  } else if (nSims > 0) {
    if (missing(fit)) stop("Cannot estimate uncertainty (nSims > 0) without a fitted model object.")
    else stop("The provided model ('fit') does not contain a covariance matrix.")
  }

  mod_spatialCovs <- spatialCovs
  mod_beta <- beta

  if (!is.null(barrier)) {
    # U = 0.5 * lambda * d^2. Log(pi) is proportional to -U.
    penalty_rast <- terra::app(barrier_sdf, fun = function(x) {
      x_work <- x / scaleFactor
      ifelse(x_work <= 0, -0.5 * lambda * (x_work^2), 0)
    })
    names(penalty_rast) <- "barrier_penalty"

    mod_spatialCovs <- c(mod_spatialCovs, penalty_rast)
    mod_beta <- c(mod_beta, 1) # The penalty acts as a fixed covariate with coefficient = 1
  }

  # --- Internal Helper: Compute a normalized probability UD ---
  calc_prob_ud_base <- function(b_vec, cov_list, maskRast) {
    ud_rast <- cov_list[[1]] * b_vec[1]
    if(length(cov_list) > 1) {
      for (j in 2:length(cov_list)) ud_rast <- ud_rast + (cov_list[[j]] * b_vec[j])
    }
    # Log-Sum-Exp for stability
    max_log <- terra::global(ud_rast, "max", na.rm = TRUE)$max
    for(k in 1:terra::nlyr(ud_rast)) ud_rast[[k]] <- exp(ud_rast[[k]] - max_log[k])

    if(!is.null(maskRast)){
      if(!inherits(maskRast, "SpatRaster")) stop("'maskRast' must be a SpatRaster")
      if (!terra::compareGeom(ud_rast, maskRast, stopOnError = FALSE)) stop("The 'maskRast' raster must share the same projection (CRS), extent, and resolution as the rasters in 'spatialCovs'.")
      m_na <- terra::ifel(maskRast <= 0, NA, 1)
      ud_rast <- terra::mask(ud_rast, m_na, maskvalues = NA)
    }

    # Normalize to sum to 1
    layer_sums <- terra::global(ud_rast, "sum", na.rm = TRUE)$sum
    for(k in 1:terra::nlyr(ud_rast)) ud_rast[[k]] <- ud_rast[[k]] / layer_sums[k]
    return(ud_rast)
  }

  ud_prob_rast <- calc_prob_ud_base(mod_beta, mod_spatialCovs, maskRast)

  if (log) {
    ud_base <- log(ud_prob_rast)
    base_name <- "log_UD"
  } else {
    ud_base <- ud_prob_rast
    base_name <- "UD"
  }

  names(ud_base) <- rep(base_name, terra::nlyr(ud_base))

  if (plot) {
    print(plotUD(ud_base, log = log))
  }

  if (!can_calc_se && nSims == 0) return(ud_base)

  # ==========================================
  # Extract components for uncertainty
  # ==========================================
  n_cells <- terra::ncell(spatialCovs[[1]])
  n_ud_layers <- max(sapply(spatialCovs, terra::nlyr))
  out_rast <- ud_base
  template <- ud_base

  # ==========================================
  # 1. Delta Method Approximation
  # ==========================================
  if (can_calc_se) {
    message("   Calculating Delta Method UD uncertainty...")
    n_covs_orig <- length(spatialCovs)

    # because the penalty term has 0 variance (fixed parameter). Its derivative zeroes out in the chain rule.
    cov_mats_orig <- lapply(spatialCovs, function(x) {
      m <- terra::as.matrix(x, wide = FALSE)
      if(ncol(m) > 1) return(m)
      return(as.vector(m))
    })

    delta_se_mat <- matrix(0, nrow = n_cells, ncol = n_ud_layers)
    ud_prob_mat <- terra::as.matrix(ud_prob_rast, wide = FALSE)

    for (k in 1:n_ud_layers) {
      pi_k <- ud_prob_mat[, k]
      C_k <- matrix(0, nrow = n_cells, ncol = n_covs_orig)

      for(j in 1:n_covs_orig) {
        if(is.matrix(cov_mats_orig[[j]])) C_k[, j] <- cov_mats_orig[[j]][, k]
        else C_k[, j] <- cov_mats_orig[[j]]
      }

      mu_C <- colSums(C_k * pi_k, na.rm = TRUE)
      C_centered <- sweep(C_k, 2, mu_C, FUN = "-")
      g_mat <- sweep(C_centered, 1, pi_k, FUN = "*")

      var_delta <- rowSums((g_mat %*% beta_cov) * g_mat)
      delta_se_mat[, k] <- sqrt(pmax(var_delta, 0))
    }

    delta_cv_mat <- delta_se_mat / ud_prob_mat
    se_rast <- terra::setValues(template, delta_se_mat)
    names(se_rast) <- rep("UD_SE_delta", n_ud_layers)
    cv_rast <- terra::setValues(template, delta_cv_mat)
    names(cv_rast) <- rep("UD_CV_delta", n_ud_layers)
    out_rast <- c(out_rast, se_rast, cv_rast)
  }

  # ==========================================
  # 2. Welford's Online Algorithm (Simulations)
  # ==========================================
  if (nSims > 0) {
    message("   Simulating ", nSims, " draws to estimate UD uncertainty...")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package \"MASS\" needed for simulation.")

    # We use mod_spatialCovs here so the C++ loop applies the barrier penalty to every draw
    n_covs_mod <- length(mod_spatialCovs)
    cov_mats_mod <- lapply(mod_spatialCovs, function(x) {
      m <- terra::as.matrix(x, wide = FALSE)
      if(ncol(m) > 1) return(m)
      return(as.vector(m))
    })

    # Draw the original betas, then append a column of 1s for the fixed penalty "beta"
    beta_draws <- MASS::mvrnorm(nSims, beta, beta_cov)
    if (!is.null(barrier)) beta_draws <- cbind(beta_draws, 1)

    cpp_res <- simulate_ud_cpp(
      nSims = nSims,
      n_cells = n_cells,
      n_ud_layers = n_ud_layers,
      n_covs = n_covs_mod,
      beta_draws = beta_draws,
      cov_mats_list = cov_mats_mod,
      show_progress = show_progress
    )

    var_pi_mat <- cpp_res$M2_pi / (nSims - 1)
    sim_se_mat <- sqrt(pmax(var_pi_mat, 0))
    sim_cv_mat <- sim_se_mat / cpp_res$mean_pi

    sim_se_rast <- terra::setValues(template, sim_se_mat)
    names(sim_se_rast) <- rep("UD_SE_sim", n_ud_layers)
    sim_cv_rast <- terra::setValues(template, sim_cv_mat)
    names(sim_cv_rast) <- rep("UD_CV_sim", n_ud_layers)
    out_rast <- c(out_rast, sim_se_rast, sim_cv_rast)
  }

  dyn_idx <- which(sapply(spatialCovs, terra::nlyr) == n_ud_layers)[1]
  if (!is.na(dyn_idx) && !is.null(terra::time(spatialCovs[[dyn_idx]]))) {
    time_vals <- terra::time(spatialCovs[[dyn_idx]])
    num_blocks <- terra::nlyr(out_rast) / n_ud_layers
    terra::time(out_rast) <- rep(time_vals, num_blocks)
  }

  return(out_rast)
}
