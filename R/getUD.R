#' Compute Utilization Distribution
#'
#' This function computes the utilization distribution (UD) for a given set of spatial covariates and habitat selection coefficients. It also estimates uncertainty in the UD using the Delta method and Monte Carlo simulations (if \code{nSims>0}). The resulting UD and uncertainty estimates are returned as \code{\link[terra]{SpatRaster-class}} objects, which can be plotted and analyzed further.
#'
#' @details
#' \strong{Calculations for \code{fitLangevin} objects:}
#' \itemize{
#'   \item \strong{UD:} Calculated as the linear predictor: \eqn{\pi(x) = \exp(c(x)\beta) / \sum \exp(c(y)\beta)}.
#'   \item \strong{Standard Errors (Delta Method):} The analytical gradient is calculated and multiplied by the Hessian-based covariance matrix of the habitat selection coefficients.
#'   \item \strong{Standard Errors (Monte Carlo):} If \code{nSims > 0}, draws are generated from the multivariate normal distribution of the coefficients \eqn{MVN(\hat{\beta}, \hat{\Sigma}_\beta)}, and the empirical standard deviation of the resulting UDs is computed on the fly using Welford's online algorithm.
#'   \item \strong{Coefficient of Variation (CV):} Computed as the standard error divided by the mean of the Monte Carlo UD draws (or the point estimate for the Delta method). Using the MC mean prevents inflation in low-probability habitats caused by the heavy upper tails of the log-normal transformation.
#' }
#'
#' \strong{Calculations for \code{hierLangevin} objects:}
#' \itemize{
#'   \item \strong{Population UD:} Because the UD is a non-linear transformation, the UD of the "average" individual is not equal to the average UD of the population (\eqn{UD(E[\beta]) \neq E[UD(\beta)]}). Therefore, the marginal population UD is integrated over the estimated between-individual variance \eqn{\mathrm{N}(\mu_\beta, \sigma^2_\beta)} using Monte Carlo integration (with \code{nSims} draws).
#'   \item \strong{Population SE:} Estimated using a nested Monte Carlo approach. Outer draws capture uncertainty in the hyperparameters \eqn{(\mu_\beta, \log \sigma_\beta)} from the Stage II working covariance matrix. For each outer draw, an inner Monte Carlo integration computes the marginal UD. The SE is the standard deviation of these marginal UDs across the outer draws.
#'   \item \strong{Individual UD (BLUP):} Computed using the specific individual's Best Linear Unbiased Predictor (BLUP) for the selection coefficients: \eqn{\beta_i = \mu_\beta + u_i}.
#'   \item \strong{Individual SE:} Draws are generated directly from the sparse joint precision matrix of the fixed effects and random effects to propagate the conditional uncertainty in the BLUPs.
#'   \item \strong{Coefficient of Variation (CV):} Computed as the standard error divided by the mean of the corresponding Monte Carlo UD draws (for both population and individual levels). As with single-track models, dividing by the MC mean stabilizes the metric in low-probability habitats where the UD point estimate approaches zero.
#' }
#'
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates. The covariates must be on the same spatial grid and have the same spatial extent.
#' @param fit A \code{fitLangevin} or \code{hierLangevin} object.
#' @param beta Numeric vector of habitat selection coefficients for the spatial covariates. The order of the coefficients must match the order of the covariates in \code{spatialCovs}.
#' @param barrier Character string. The name of the barrier in \code{spatialCovs} that is represented as a signed distance field (see \code{\link{prepBarrier}}). If provided, this raster is exclusively used for the barrier penalty and is not included in the habitat selection covariates. Default: \code{NULL} (no barrier). If \code{fit} is provided, this is extracted automatically.
#' @param lambda Numeric. The penalty weight for the barrier constraint. Default: \code{NULL}. If \code{fit} is provided, this is extracted automatically.
#' @param scaleFactor Internal scaling factor for the coordinates and parameters (see \code{\link{fitLangevin}}). Automatically extracted if \code{fit} is provided. Default: 1 (no scaling).
#' @param log Logical indicating whether or not to return the log of the utilization distribution. Default: \code{TRUE}.
#' @param nSims Integer. Number of draws from the covariance matrix to use for estimating Monte Carlo uncertainty in the UD. For \code{hierLangevin} objects, this dictates the number of draws from the between-individual MVN distribution used to integrate the marginal population-level UD, and \strong{must} be > 0. If \code{nSims > 0}, the returned raster stack will include additional layers for simulated SE and CV. Default: \code{0}.
#' @param individual Optional vector of individual IDs (character or numeric) to calculate individual-level UDs for hierarchical models. If \code{NULL} (default), only the population-level UD is calculated. For example, \code{individual = c(1,3)}, \code{individual = c("ID1","ID3")}, or \code{individual = "all"}. Default: \code{NULL}.
#' @param show_progress Logical. If \code{TRUE}, displays a progress bar for simulations. Default: \code{TRUE}.
#' @param plot Logical. Plot the resulting UD using \code{\link{plotUD}}? Default: \code{TRUE}.
#' @param maskRast \code{\link[terra]{SpatRaster-class}} object for areas to be masked out (set to \code{NA}) before plotting the UD. Default: \code{NULL} (no mask).
#' @param extent Optional. A numeric vector of length 4 \code{c(xmin, xmax, ymin, ymax)} or a \code{\link[terra]{SpatExtent}} object defining the bounding box. If provided, the returned UD is cropped to this extent. Default: \code{NULL}.
#' @param normalize Logical. If \code{TRUE} and \code{extent} is provided, the UD is normalized specifically across the \code{extent} (ignoring cells outside the extent). If \code{FALSE}, the UD is normalized globally over the extent of \code{spatialCovs} before any cropping. Default: \code{FALSE}.
#' @return A \code{\link[terra]{SpatRaster}} object. It contains the (log) utilization distribution. For hierarchical models, the stack contains the population-level UD and any specified individual-level UDs. It will also contain layers for Delta method standard errors (\code{UD_SE_delta}) and CVs (\code{UD_CV_delta}) if the covariance matrix is available. If \code{nSims > 0}, it adds simulated layers (\code{UD_SE_sim}, \code{UD_CV_sim}).
#' @seealso \code{\link{plotUD}}, \code{\link{regionProb}}.
#' @examples
#' # exampleCovs included in package; see ?exampleCovs for details
#' UD <- getUD(exampleCovs, beta = c(-4, 6, 5, -0.1) )
#' @importFrom terra global nlyr varnames app setValues compareGeom mask crop ext
#' @importFrom stats setNames
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
getUD <- function(spatialCovs, fit, beta, barrier = NULL, lambda = NULL, scaleFactor = 1, log = TRUE, nSims = 0, individual = NULL, show_progress = TRUE, plot = TRUE, maskRast = NULL, extent = NULL, normalize = FALSE) {

  if((missing(fit) & missing(beta)) | (!missing(fit) & !missing(beta))) stop("Either 'fit' or 'beta' must be provided, but not both.")

  is_hierLangevin <- !missing(fit) && inherits(fit, "hierLangevin")

  if (is_hierLangevin && nSims <= 0) {
    stop("For 'hierLangevin' models, population-level UDs must be integrated over the between-individual variance using Monte Carlo methods. Please specify 'nSims > 0' (e.g., nSims = 1000).")
  }

  if(!missing(fit)) {
    verify_signatures(fit, spatialCovs = spatialCovs)
    if(!is_hierLangevin) boundsWarning(fit)
  }

  if(!missing(fit)) {
    if (!is.null(fit$conditions$barrier)) barrier <- fit$conditions$barrier
    if (!is.null(fit$conditions$lambda)) lambda <- fit$conditions$lambda
    if (!is.null(fit$conditions$scaleFactor)) scaleFactor <- fit$conditions$scaleFactor
  }

  # Strip the barrier out to allow the length check to pass, while preserving the modified stack for calculations
  bar_info <- .prep_barrier_raster(spatialCovs, barrier, lambda, scaleFactor)
  spatialCovs <- bar_info$spatialCovs
  mod_spatialCovs <- bar_info$mod_spatialCovs

  if (!missing(beta) && length(spatialCovs) != length(beta)) {
    stop("length(spatialCovs) must equal length(beta). Note that if a barrier is specified, the barrier does not receive a beta coefficient.")
  }

  beta_info <- .extract_langevin_beta_list(fit = fit, beta = beta, individual = individual)
  beta_list <- beta_info$beta_list
  hierarchical_logical <- beta_info$hierarchical_logical
  sd_beta_vec <- beta_info$sd_beta_vec

  cov_info <- .extract_langevin_cov(fit, is_hierLangevin, hierarchical_logical)
  beta_cov <- cov_info$beta_cov
  beta_idx_cov <- cov_info$beta_idx_cov
  can_calc_se <- cov_info$can_calc_se

  if (!can_calc_se && nSims > 0) {
    if (missing(fit)) stop("Cannot estimate uncertainty (nSims > 0) without a fitted model object.")
    else stop("The provided model ('fit') does not contain a valid covariance matrix.")
  }

  run_sim_se <- can_calc_se && (nSims > 0)

  if (!is.null(extent) && normalize) {
    crop_ext <- tryCatch(terra::ext(extent), error = function(e) NULL)
    if (!is.null(crop_ext)) {
      spatialCovs <- lapply(spatialCovs, function(x) terra::crop(x, crop_ext))
      mod_spatialCovs <- lapply(mod_spatialCovs, function(x) terra::crop(x, crop_ext))
      if (!is.null(maskRast)) {
        maskRast <- terra::crop(maskRast, crop_ext)
      }
    }
  }

  calc_prob_ud_base <- function(b_vec, cov_list, maskRast) {
    ud_rast <- cov_list[[1]] * b_vec[1]
    if(length(cov_list) > 1) {
      for (j in 2:length(cov_list)) ud_rast <- ud_rast + (cov_list[[j]] * b_vec[j])
    }
    max_log <- terra::global(ud_rast, "max", na.rm = TRUE)$max
    for(k in 1:terra::nlyr(ud_rast)) ud_rast[[k]] <- exp(ud_rast[[k]] - max_log[k])

    if(!is.null(maskRast)){
      if(!inherits(maskRast, "SpatRaster")) stop("'maskRast' must be a SpatRaster")
      if (!terra::compareGeom(ud_rast, maskRast, stopOnError = FALSE)) stop("The 'maskRast' raster must share the same projection (CRS), extent, and resolution as the rasters in 'spatialCovs'.")
      m_na <- terra::ifel(maskRast <= 0, NA, 1)
      ud_rast <- terra::mask(ud_rast, m_na, maskvalues = NA)
    }

    layer_sums <- terra::global(ud_rast, "sum", na.rm = TRUE)$sum
    for(k in 1:terra::nlyr(ud_rast)) ud_rast[[k]] <- ud_rast[[k]] / layer_sums[k]
    return(ud_rast)
  }

  n_cells <- terra::ncell(spatialCovs[[1]])
  n_ud_layers <- max(sapply(spatialCovs, terra::nlyr))
  n_covs_mod <- length(mod_spatialCovs)

  cov_mats_mod <- lapply(mod_spatialCovs, function(x) {
    m <- terra::as.matrix(x, wide = FALSE)
    if(ncol(m) > 1) return(m)
    return(as.vector(m))
  })

  base_names <- names(beta_list)
  ud_list <- lapply(seq_along(beta_list), function(i) {
    nm <- base_names[i]
    item <- beta_list[[i]]

    if (is.list(item) && !is.null(item$type) && item$type == "hierLangevin_population") {
      n_mc <- max(1000, nSims)
      message("   Using Monte Carlo integration (", n_mc, " draws) to compute marginal Population UD...")

      pop_cov <- diag(item$sd^2, nrow = length(item$sd))
      beta_draws <- as.matrix(MASS::mvrnorm(n_mc, item$mu, pop_cov))

      if (!is.null(barrier)) beta_draws <- cbind(beta_draws, 1)

      cpp_res <- simulate_ud_cpp(
        nSims = n_mc,
        n_cells = n_cells,
        n_ud_layers = n_ud_layers,
        n_covs = n_covs_mod,
        beta_draws = beta_draws,
        cov_mats_list = cov_mats_mod,
        show_progress = FALSE
      )

      ud_prob <- spatialCovs[[1]][[1:n_ud_layers]]
      ud_prob <- terra::setValues(ud_prob, cpp_res$mean_pi)

      if(!is.null(maskRast)){
        if(!inherits(maskRast, "SpatRaster")) stop("'maskRast' must be a SpatRaster")
        if (!terra::compareGeom(ud_prob, maskRast, stopOnError = FALSE)) stop("The 'maskRast' raster must share the same projection (CRS), extent, and resolution as the rasters in 'spatialCovs'.")
        m_na <- terra::ifel(maskRast <= 0, NA, 1)
        ud_prob <- terra::mask(ud_prob, m_na, maskvalues = NA)
      }

      layer_sums <- terra::global(ud_prob, "sum", na.rm = TRUE)$sum
      for(k in 1:terra::nlyr(ud_prob)) ud_prob[[k]] <- ud_prob[[k]] / layer_sums[k]

    } else {
      b_vec <- if (is.list(item)) item$vec else item
      if (!is.null(barrier)) b_vec <- c(b_vec, 1)
      ud_prob <- calc_prob_ud_base(b_vec, mod_spatialCovs, maskRast)
    }

    if (log) {
      ud_res <- log(ud_prob)
      pref <- "log_UD_"
    } else {
      ud_res <- ud_prob
      pref <- "UD_"
    }

    if (nm == "UD") {
      names(ud_res) <- rep(if(log) "log_UD" else "UD", terra::nlyr(ud_res))
    } else {
      names(ud_res) <- rep(paste0(pref, nm), terra::nlyr(ud_res))
    }
    return(list(prob = ud_prob, res = ud_res))
  })

  # Simplify beta_list back to standard vectors for downstream SE logic
  beta_list <- lapply(beta_list, function(item) {
    if (is.list(item) && !is.null(item$type) && item$type == "hierLangevin_population") {
      return(item$mu)
    }
    return(item)
  })

  ud_base <- do.call(c, lapply(ud_list, function(x) x$res))

  if (plot) print(plotUD(ud_base, log = log, extent = extent))

  if (!can_calc_se && !run_sim_se) {
    out_rast <- ud_base
    if (!is.null(extent) && !normalize) {
      crop_ext <- tryCatch(terra::ext(extent), error = function(e) NULL)
      if (!is.null(crop_ext)) out_rast <- terra::crop(out_rast, crop_ext)
    }
    return(out_rast)
  }

  out_rast <- ud_base
  template <- ud_list[[1]]$res

  if (can_calc_se && !is_hierLangevin) {
    if (hierarchical_logical) {
      message("   Calculating Delta Method UD uncertainty (Population-level)...")
    } else {
      message("   Calculating Delta Method UD uncertainty...")
    }

    n_covs_orig <- length(spatialCovs)

    cov_mats_orig <- lapply(spatialCovs, function(x) {
      m <- terra::as.matrix(x, wide = FALSE)
      if(ncol(m) > 1) return(m)
      return(as.vector(m))
    })

    for (i in seq_along(beta_list)) {
      nm <- names(beta_list)[i]

      if (hierarchical_logical && nm != "Population") {
        message("     Skipping Delta method for ", nm, " (use nSims > 0 for individual-level uncertainty).")
        next
      }

      ud_prob_mat <- terra::as.matrix(ud_list[[i]]$prob, wide = FALSE)
      delta_se_mat <- matrix(0, nrow = n_cells, ncol = n_ud_layers)

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
      pref_se <- if(nm == "UD") "UD_SE_delta" else paste0("UD_SE_delta_", nm)
      names(se_rast) <- rep(pref_se, n_ud_layers)

      cv_rast <- terra::setValues(template, delta_cv_mat)
      pref_cv <- if(nm == "UD") "UD_CV_delta" else paste0("UD_CV_delta_", nm)
      names(cv_rast) <- rep(pref_cv, n_ud_layers)

      out_rast <- c(out_rast, se_rast, cv_rast)
    }
  } else if (can_calc_se && is_hierLangevin) {
    message("   Skipping Delta Method (mathematically invalid for hierLangevin marginal UDs). Relying on Monte Carlo...")
  }

  if (run_sim_se) {
    if (hierarchical_logical) {
      message("   Simulating ", nSims, " draws to estimate hierarchical UD uncertainties...")
    } else {
      message("   Simulating ", nSims, " draws to estimate UD uncertainty...")
    }

    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package \"MASS\" needed for simulation.")

    draw_data <- if (hierarchical_logical && length(beta_list) > 1) {
      .sample_joint_precision_re(fit, is_hierLangevin, beta_idx_cov, nSims)
    } else {
      list(has_Q = FALSE)
    }

    for (i in seq_along(beta_list)) {
      nm <- names(beta_list)[i]

      skip_standard_sim <- FALSE

      if (hierarchical_logical) {
        if (nm == "Population") {
          if (is_hierLangevin) {
            message("     Running nested Monte Carlo to estimate Population UD uncertainty...")

            pop_params_draws <- MASS::mvrnorm(nSims, c(beta_list[[i]], log(sd_beta_vec)), beta_cov)
            M_inner <- 30

            mean_ud <- numeric(n_cells * n_ud_layers)
            M2_ud <- numeric(n_cells * n_ud_layers)

            if (show_progress) pb <- txtProgressBar(min = 0, max = nSims, style = 3, width=50)

            for (s in 1:nSims) {
              mu_s <- pop_params_draws[s, 1:length(beta_list[[i]])]
              sd_s <- exp(pop_params_draws[s, (length(beta_list[[i]])+1):ncol(pop_params_draws)])

              inner_cov <- diag(sd_s^2, nrow = length(sd_s))
              inner_beta_draws <- as.matrix(MASS::mvrnorm(M_inner, mu_s, inner_cov))
              if (!is.null(barrier)) inner_beta_draws <- cbind(inner_beta_draws, 1)

              cpp_res_s <- simulate_ud_cpp(
                nSims = M_inner, n_cells = n_cells, n_ud_layers = n_ud_layers,
                n_covs = n_covs_mod, beta_draws = inner_beta_draws,
                cov_mats_list = cov_mats_mod, show_progress = FALSE
              )

              val <- cpp_res_s$mean_pi

              delta <- val - mean_ud
              mean_ud <- mean_ud + delta / s
              delta2 <- val - mean_ud
              M2_ud <- M2_ud + delta * delta2

              if (show_progress) setTxtProgressBar(pb, s)
            }
            if (show_progress) close(pb)

            sim_se_mat <- matrix(sqrt(M2_ud / (nSims - 1)), nrow = n_cells, ncol = n_ud_layers)
            sim_cv_mat <- sim_se_mat / matrix(mean_ud, nrow = n_cells, ncol = n_ud_layers)
            skip_standard_sim <- TRUE

          } else {
            beta_draws <- MASS::mvrnorm(nSims, beta_list[[i]], beta_cov)
          }
        } else {
          if (!draw_data$has_Q) next
          beta_draws <- .get_individual_beta_draws(nm, fit, is_hierLangevin, draw_data, beta_idx_cov, nSims)
        }
      } else {
        beta_draws <- MASS::mvrnorm(nSims, beta_list[[i]], beta_cov)
      }

      if (!skip_standard_sim) {
        if (hierarchical_logical && nm != "Population") {
          message("     Simulating UD uncertainty for individual: ", sub("^ID_", "", nm), "...")
        }

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
      }

      sim_se_rast <- terra::setValues(template, sim_se_mat)
      pref_se <- if(nm == "UD") "UD_SE_sim" else paste0("UD_SE_sim_", nm)
      names(sim_se_rast) <- rep(pref_se, n_ud_layers)

      sim_cv_rast <- terra::setValues(template, sim_cv_mat)
      pref_cv <- if(nm == "UD") "UD_CV_sim" else paste0("UD_CV_sim_", nm)
      names(sim_cv_rast) <- rep(pref_cv, n_ud_layers)

      out_rast <- c(out_rast, sim_se_rast, sim_cv_rast)
    }
  }

  dyn_idx <- which(sapply(spatialCovs, terra::nlyr) == n_ud_layers)[1]
  if (!is.na(dyn_idx) && !is.null(terra::time(spatialCovs[[dyn_idx]]))) {
    time_vals <- terra::time(spatialCovs[[dyn_idx]])
    num_blocks <- terra::nlyr(out_rast) / n_ud_layers
    terra::time(out_rast) <- rep(time_vals, num_blocks)
  }

  if (!is.null(extent) && !normalize) {
    crop_ext <- tryCatch(terra::ext(extent), error = function(e) NULL)
    if (!is.null(crop_ext)) out_rast <- terra::crop(out_rast, crop_ext)
  }

  return(out_rast)
}
