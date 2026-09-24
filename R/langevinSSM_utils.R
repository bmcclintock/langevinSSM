#' @importFrom utils globalVariables
utils::globalVariables(c("mu.x", "mu.y", "vel.x", "vel.y", "id", "psi", "tau", "dt", "x", "y", "smaj", "smin", "eor", "x.err", "y.err", "val", "lag", "UD", "type", "theoretical", "uid", "i", "split_flag", "segment","max_seg","mu.x_pr","mu.y_pr",".data","temp_x","temp_y","pr_idx","date_num"))

#' Example Spatial Covariates
#'
#' A list of SpatRaster objects used as example covariates.
#' @name exampleCovs
#' @docType data
#' @export
NULL

.onLoad <- function(libname, pkgname) {
  env <- asNamespace(pkgname)

  makeActiveBinding("exampleCovs", function() {

    cov_names <- c("exampleCov1", "exampleCov2", "exampleCov3", "exampleCov4")
    lyr_names <- c("cov1","cov2","cov3","d2c")
    names(lyr_names) <- cov_names

    cov_list <- lapply(cov_names, function(cov) {

      file_path <- system.file("extdata", paste0(cov, ".tif"), package = "langevinSSM")

      # Safety check
      if (file_path == "") {
        stop(paste("Could not find", paste0(cov, ".tif"), "in the extdata/ folder."))
      }

      r <- terra::rast(file_path)
      names(r) <- lyr_names[cov]
      return(r)
    })

    names(cov_list) <- lyr_names

    return(cov_list)

  }, env)
}

#' Example formatted tracking data
#'
#' A \code{dataLangevin} object containing formatted movement tracks appropriate for \code{\link{fitLangevin}}..
#'
#' @name exampleDat
#' @docType data
NULL

#' Example unformatted tracking Data
#'
#' A data frame containing example movement tracks appropriate for \code{\link{formatData}}.
#'
#' @name unformatDat
#' @docType data
NULL

#' S3 Methods for Langevin model fits
#'
#' Standard S3 methods for extracting information from \code{fitLangevin} objects.
#'
#' @param object A \code{fitLangevin} object.
#' @param type Character string indicating which scale or type of estimates to extract. Options are \code{"natural"} (default), \code{"working"}, or \code{"random"}. Ignored for \code{fitted}.
#' @param parm A specification of which parameters to extract. For fixed parameters (when \code{type} is \code{"natural"} or \code{"working"}), this can be a vector of parameter names, or \code{"fixed"} (default) to return all fixed parameters. For random effects (when \code{type = "random"}, or when using \code{fitted}), this must be either \code{"mu"} (default) or \code{"vel"}.
#' @param level The confidence level required. Default: \code{0.95}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' \itemize{
#'   \item \strong{\code{logLik}:} Returns an object of class \code{logLik} containing the maximized log-likelihood, with attributes for degrees of freedom (\code{df}) and number of observations (\code{nobs}).
#'   \item \strong{\code{coef}:} Returns a named numeric vector of parameter point estimates. If \code{type = "random"}, returns a data frame of the latent states.
#'   \item \strong{\code{vcov}:} Returns the variance-covariance matrix of the estimated parameters. If \code{type = "random"}, returns the sparse joint precision matrix of the random effects.
#'   \item \strong{\code{confint}:} For fixed parameters, returns a numeric matrix with lower and upper bounds. For random effects, returns a wide data frame containing the \code{id}, \code{date}, and the upper and lower confidence bounds for the coordinates.
#'   \item \strong{\code{fitted}:} Returns a data frame containing the \code{id}, \code{date}, and the estimated expected values of the latent states.
#'   \item \strong{\code{summary}:} Returns a \code{summary.fitLangevin} object containing convergence details and matrices of coefficients with standard errors (and Z-test p-values for habitat selection coefficients).
#' }
#'
#' @name langevin_methods
NULL

#' Print a fitLangevin object
#'
#' @param x A \code{fitLangevin} object returned by \code{\link{fitLangevin}}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @rdname langevin_methods
#' @importFrom stats printCoefmat
#' @export
print.fitLangevin <- function(x, ...) {

  cat("\nHabitat-Driven Langevin Diffusion Model\n")
  cat("=======================================\n")

  # Determine model type
  model_type <- ifelse("vel" %in% names(x$estimates$random), "Underdamped", "Overdamped")
  cat("Model type:       ", model_type, "\n")

  # Convergence status
  conv_text <- ifelse(x$convergence == 0, "Successful", paste("Failed (Code", x$convergence, ")"))
  cat("Convergence:      ", conv_text, "\n")
  if (x$convergence != 0 && !is.null(x$message)) {
    cat("Message:          ", x$message, "\n")
  }

  cat("Max Log-Likelihood:", -x$objective, "\n")
  cat("Optimization time: ", round(x$elapsedTime[3], 2), "seconds\n")

  if (!is.null(x$conditions$barrier) && !is.null(x$conditions$lambda)) {
    cat("Barrier penalty:   ", signif(x$conditions$lambda, 4), "\n")
  }
  cat("\n")

  cat("Parameter Estimates (Natural Scale):\n")
  cat("---------------------------------------\n")

  # Clean up the rownames for the natural estimates matrix
  nat_est <- x$estimates$natural

  stats::printCoefmat(nat_est, digits = 4, signif.stars = FALSE, na.print = "NA", ...)

  boundsWarning(x, as_warning = FALSE)

  invisible(x)
}

#' @rdname langevin_methods
#' @export
logLik.fitLangevin <- function(object, ...) {

  boundsWarning(object, as_warning = FALSE)

  # TMB minimizes the negative log-likelihood
  val <- -object$objective

  # The number of estimated parameters (degrees of freedom)
  attr(val, "df") <- length(object$par)

  # The number of observations (required for BIC)
  if (!is.null(object$signatures$data$nrow)) {
    attr(val, "nobs") <- object$signatures$data$nrow
  }

  class(val) <- "logLik"
  return(val)
}

#' @rdname langevin_methods
#' @export
coef.fitLangevin <- function(object, type = c("natural", "working"), ...) {

  type <- match.arg(type)
  boundsWarning(object, as_warning = FALSE)

  if (!type %in% names(object$estimates)) {
    stop("Estimates for type '", type, "' are not available in this model fit.")
  }

  est_vector <- object$estimates[[type]][, "Estimate"]
  names(est_vector) <- rownames(object$estimates[[type]])

  return(est_vector)
}

#' @rdname langevin_methods
#' @export
fitted.fitLangevin <- function(object, parm = c("mu", "vel"), ...) {

  parm <- match.arg(parm)
  boundsWarning(object, as_warning = FALSE)

  if (is.null(object$estimates$random[[parm]])) {
    stop("Latent state '", parm, "' was not estimated in this model.")
  }

  return(object$estimates$random[[parm]]$est)
}

#' @rdname langevin_methods
#' @export
vcov.fitLangevin <- function(object, parm = "fixed", type = c("natural", "working", "random"), ...) {

  if (!missing(parm) && any(parm %in% c("mu", "vel"))) type <- "random"
  type <- match.arg(type)

  boundsWarning(object, as_warning = FALSE)

  if (type %in% c("natural", "working")) {
    if (!type %in% names(object$covariance)) {
      stop("Covariance matrix for type '", type, "' is not available in this model fit.")
    }
    vc <- object$covariance[[type]]

    if (!"fixed" %in% parm) {
      missing_parms <- parm[!parm %in% rownames(vc)]
      if (length(missing_parms) > 0) {
        stop("Parameter(s) not found in model: ", paste(missing_parms, collapse = ", "))
      }
      vc <- vc[parm, parm, drop = FALSE]
    }
    return(vc)

  } else if (type == "random") {
    if (is.null(object$covariance$random$jointPrecision)) {
      stop("Joint precision matrix not available. Re-run model with getJointPrecision = TRUE.")
    }
    return(object$covariance$random$jointPrecision)
  }
}

#' @rdname langevin_methods
#' @export
confint.fitLangevin <- function(object, parm = "fixed", level = 0.95, type = c("natural", "working", "random"), ...) {

  if (!missing(parm) && any(parm %in% c("mu", "vel"))) type <- "random"
  type <- match.arg(type)

  boundsWarning(object, as_warning = FALSE)

  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
  fac <- stats::qnorm(a)

  if (type %in% c("natural", "working")) {
    cf_all <- object$estimates[[type]][, "Estimate"]
    names(cf_all) <- rownames(object$estimates[[type]])
    se_all <- sqrt(diag(object$covariance[[type]]))

    ci <- array(NA, dim = c(length(cf_all), 2L), dimnames = list(names(cf_all), pct))
    ci[, 1] <- cf_all + fac[1] * se_all
    ci[, 2] <- cf_all + fac[2] * se_all

    if (!"fixed" %in% parm) {
      missing_parms <- parm[!parm %in% rownames(ci)]
      if (length(missing_parms) > 0) {
        stop("Parameter(s) not found in model: ", paste(missing_parms, collapse = ", "))
      }
      ci <- ci[parm, , drop = FALSE]
    }
    return(ci)

  } else if (type == "random") {
    if ("fixed" %in% parm) parm <- "mu" # Default
    if (length(parm) > 1) stop("For type='random', please specify a single parm ('mu' or 'vel').")

    if (!parm %in% names(object$estimates$random)) {
      stop("Random effect '", parm, "' is not available in this model fit.")
    }

    est_df <- object$estimates$random[[parm]]$est
    se_df  <- object$estimates$random[[parm]]$se

    x_col <- paste0(parm, ".x")
    y_col <- paste0(parm, ".y")

    pct_clean <- gsub(" ", "", pct)

    ci_df <- data.frame(id = est_df$id)
    if ("date" %in% names(est_df)) {
      ci_df$date <- est_df$date
    }

    ci_df[[paste0(x_col, "_", pct_clean[1])]] <- est_df[[x_col]] + fac[1] * se_df[[x_col]]
    ci_df[[paste0(x_col, "_", pct_clean[2])]] <- est_df[[x_col]] + fac[2] * se_df[[x_col]]

    ci_df[[paste0(y_col, "_", pct_clean[1])]] <- est_df[[y_col]] + fac[1] * se_df[[y_col]]
    ci_df[[paste0(y_col, "_", pct_clean[2])]] <- est_df[[y_col]] + fac[2] * se_df[[y_col]]

    return(ci_df)
  }
}

#' @rdname langevin_methods
#' @export
summary.fitLangevin <- function(object, ...) {

  boundsWarning(object, as_warning = FALSE)

  nat_est <- object$estimates$natural
  beta_idx <- grepl("^beta", rownames(nat_est))

  coef_beta <- nat_est[beta_idx, , drop = FALSE]

  se_zero <- coef_beta[, "Std. Error"] == 0
  z_val <- ifelse(se_zero, NA_real_, coef_beta[, "Estimate"] / coef_beta[, "Std. Error"])
  p_val <- ifelse(se_zero, NA_real_, 2 * stats::pnorm(abs(z_val), lower.tail = FALSE))

  coef_beta$z_value <- z_val
  coef_beta$`Pr(>|z|)` <- p_val

  coef_process <- nat_est[!beta_idx, , drop = FALSE]

  res <- list(
    model_type = ifelse("vel" %in% names(object$estimates$random), "Underdamped", "Overdamped"),
    convergence = object$convergence,
    message = object$message,
    loglik = -object$objective,
    elapsed = object$elapsedTime[3],
    coef_beta = as.matrix(coef_beta), # printCoefmat prefers matrices
    coef_process = as.matrix(coef_process),
    barrier = object$conditions$barrier,
    lambda = object$conditions$lambda
  )

  class(res) <- "summary.fitLangevin"
  return(res)
}

#' Print summary of a fitLangevin object
#'
#' @param x A \code{summary.fitLangevin} object.
#' @param digits Minimal number of significant digits, see \code{\link[stats]{printCoefmat}}.
#' @param signif.stars Logical. See \code{\link[stats]{printCoefmat}}.
#' @param ... further arguments passed to \code{\link{print.default}}.
#'
#' @importFrom stats printCoefmat
#' @export
print.summary.fitLangevin <- function(x, digits = 4, signif.stars = TRUE, ...) {

  cat("\nHabitat-Driven Langevin Diffusion Model\n")
  cat("=======================================\n")
  cat("Model type:       ", x$model_type, "\n")

  conv_text <- ifelse(x$convergence == 0, "Successful", paste("Failed (Code", x$convergence, ")"))
  cat("Convergence:      ", conv_text, "\n")
  if (x$convergence != 0 && !is.null(x$message)) {
    cat("Message:          ", x$message, "\n")
  }

  cat("Max Log-Likelihood:", x$loglik, "\n")
  cat("Optimization time: ", round(x$elapsed, 2), "seconds\n")

  if (!is.null(x$barrier) && !is.null(x$lambda)) {
    cat("Barrier penalty:   ", signif(x$lambda, 4), "\n")
  }
  cat("\n")

  cat("Habitat Selection Coefficients:\n")
  cat("---------------------------------------\n")
  stats::printCoefmat(x$coef_beta, digits = digits, signif.stars = signif.stars,
                      na.print = "NA", has.Pvalue = TRUE, ...)

  cat("\nProcess & Observation Parameters (Natural Scale):\n")
  cat("---------------------------------------\n")
  stats::printCoefmat(x$coef_process, digits = digits, signif.stars = FALSE,
                      na.print = "NA", ...)

  invisible(x)
}

#' Print a resLangevin object
#' @method print resLangevin
#'
#' @param x A \code{resLangevin} object returned by \code{\link{residuals.fitLangevin}}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @export
print.resLangevin <- function(x, ...) {

  cat("\n=== One-Step-Ahead (OSA) Residuals ===\n")

  n_obs <- nrow(x)
  n_tracks <- length(unique(x$id))
  cat("Total observations:", n_obs, "\n")
  cat("Number of tracks:  ", n_tracks, "\n\n")

  tests_df <- attr(x, "tests")
  if(!is.null(tests_df)){
    cat("---- Goodness-of-Fit Tests ----\n")
    print(tests_df, row.names = FALSE)
    cat("-------------------------------\n\n")
  }

  cat("Residual Summary:\n")
  res_summary <- summary(x[, c("residual.x", "residual.y")])
  print(res_summary)

  #cat("\n* Tip: Use plot() on this object to view diagnostic plots.\n")

  invisible(x)
}

#' Print a regLangevin object
#'
#' @param x A \code{regLangevin} object returned by \code{\link{regionProb}}.
#' @param digits Minimal number of significant digits to print. Default: \code{4}.
#' @param header Logical. Indicates whether to print the main descriptive header. Default: \code{TRUE}.
#' @param ... Additional arguments passed to \code{print}.
#' @export
print.regLangevin <- function(x, digits = 4, header = TRUE, ...) {
  # Handle container lists (e.g., when regionProb returns population + individual results)
  if (!("Point_Estimate" %in% names(x)) && is.list(x) && length(x) > 0) {
    item_names <- names(x)
    for (i in seq_along(x)) {
      nm <- if (!is.null(item_names)) item_names[i] else paste("Item", i)
      clean_nm <- sub("^ID_", "Individual: ", nm)
      cat("\n============================================\n")
      cat("Regional Probability -", clean_nm, "\n")
      cat("============================================\n")
      print.regLangevin(x[[i]], digits = digits, header = FALSE, ...)
    }
    return(invisible(x))
  }

  n_layers <- length(x$Point_Estimate)
  conf_level <- x$level * 100

  if (n_layers == 1) {
    if (isTRUE(header)) {
      cat("Regional Probability Estimate\n")
      cat("=============================\n")
    }
    cat(sprintf("Point Estimate: %.*f\n\n", digits, x$Point_Estimate[1]))

    if (!is.null(x$SE_delta) && !all(is.na(x$SE_delta))) {
      cat("Delta Method Approximation:\n")
      cat(sprintf("  Standard Error: %.*f\n", digits, x$SE_delta[1]))
      cat(sprintf("  %.0f%% CI:         [%.*f, %.*f]\n", conf_level, digits, x$CI_delta[1, 1], digits, x$CI_delta[1, 2]))
      cat("\n")
    }

    if (!is.null(x$SE_sim)) {
      cat("Monte Carlo Simulation:\n")
      cat(sprintf("  Standard Error: %.*f\n", digits, x$SE_sim[1]))
      cat(sprintf("  %.0f%% CI:         [%.*f, %.*f]\n", conf_level, digits, x$CI_sim[1, 1], digits, x$CI_sim[1, 2]))
      cat(sprintf("  (Based on %d draws)\n", nrow(x$simulated_draws)))
    }

  } else {
    if (isTRUE(header)) {
      cat("Regional Probability Estimates (Multi-Layer)\n")
      cat("============================================\n")
    }

    time_vals <- tryCatch(terra::time(x$prob_raster), error = function(e) NULL)
    has_time <- !is.null(time_vals) && !all(is.na(time_vals))

    df_base <- data.frame(Layer = seq_len(n_layers))
    if (has_time) df_base$Time <- time_vals

    # 1. Point Estimates
    df_est <- df_base
    df_est$Estimate <- round(x$Point_Estimate, digits)
    cat("\n--- Point Estimates ---\n")
    print(df_est, row.names = FALSE)

    # 2. Delta Method
    if (!is.null(x$SE_delta) && !all(is.na(x$SE_delta))) {
      df_delta <- df_base
      df_delta$SE <- round(x$SE_delta, digits)
      df_delta$CI <- sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), x$CI_delta[, 1], x$CI_delta[, 2])
      names(df_delta)[names(df_delta) == "CI"] <- sprintf("%.0f%%_CI", conf_level)
      cat("\n--- Delta Method Approximation ---\n")
      print(df_delta, row.names = FALSE)
    }

    # 3. Monte Carlo
    if (!is.null(x$SE_sim)) {
      df_sim <- df_base
      df_sim$SE <- round(x$SE_sim, digits)
      df_sim$CI <- sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), x$CI_sim[, 1], x$CI_sim[, 2])
      names(df_sim)[names(df_sim) == "CI"] <- sprintf("%.0f%%_CI", conf_level)

      cat("\n--- Monte Carlo Simulation ---\n")
      print(df_sim, row.names = FALSE)
      cat(sprintf("\n(Based on %d draws)\n", nrow(x$simulated_draws)))
    }
  }

  invisible(x)
}

rasterList <- function (rast)
{
  lim <- as.vector(terra::ext(rast))
  res <- terra::res(rast)
  xgrid <- seq(lim[1] + res[1]/2, lim[2] - res[1]/2, by = res[1])
  ygrid <- seq(lim[3] + res[2]/2, lim[4] - res[2]/2, by = res[2])

  # Add wide = TRUE so terra returns [nrow, ncol] instead of [ncells, nlayers]
  z_mat <- terra::as.matrix(rast, wide = TRUE)
  z <- t(apply(z_mat, 2, rev))

  return(list(x = xgrid, y = ygrid, z = z))
}

## modified from aniMotum version 1.2-15
#' @importFrom sf st_as_sf st_crs
#' @importFrom dplyr tibble
format_data <- function(x, id = "id", date = "date", lc = "lc", coord = c("x", "y"), epar = c("smaj", "smin", "eor"), sderr = c("x.err", "y.err"), tz = "UTC") {

  if (id %in% names(x))
    stopifnot(`id must be a character string` = is.character(id))
  else stop("An 'id' variable must be included in the input data\n")

  stopifnot(`date must be a character string` = is.character(date))
  stopifnot(`lc must be a character string` = is.character(lc))
  stopifnot(`coord must be a character vector with 2 elements` = all(is.character(coord)) && length(coord) == 2)
  stopifnot(`epar must be a character vector with 3 elements` = all(is.character(epar)) && length(epar) == 3)
  stopifnot(`sderr must be a character vector with 2 elements` = all(is.character(sderr)) && length(sderr) == 2)

  stopifnot(`An id variable must be included in the input data` = id %in% names(x))
  stopifnot(`A date/time variable must be included in the input data` = date %in% names(x))
  stopifnot(`Coordinate variables must be included in the input data` = all(coord %in% names(x)))

  # --- DYNAMIC COLUMN INJECTION ---
  if (all(!epar %in% names(x))) {
    x[[epar[1]]] <- as.double(NA)
    x[[epar[2]]] <- as.double(NA)
    x[[epar[3]]] <- as.double(NA)
  }
  if (all(!sderr %in% names(x))) {
    x[[sderr[1]]] <- as.double(NA)
    x[[sderr[2]]] <- as.double(NA)
  }

  xt.vars <- names(x)[!names(x) %in% c(id, date, lc, coord, epar, sderr)]

  # Subset cleanly
  xx <- x[, c(id, date, lc, coord, epar, sderr, xt.vars)]

  # Force standard names
  new_names <- c("id", "date", "lc", coord, "smaj", "smin", "eor", "x.err", "y.err", xt.vars)
  names(xx) <- new_names

  if (is.factor(xx$id)) xx$id <- droplevels(xx$id)
  xx$id <- as.character(xx$id)

  if (!inherits(xx$date, "POSIXt")) {
    xx$date <- try(as.POSIXct(xx$date, tz = tz), silent = TRUE)
    if (inherits(xx$date, "try-error"))
      stop("dates must be in a standard format: YYYY-MM-DD HH:MM:SS")
  }

  xx <- xx[order(xx$date), ]
  return(xx)
}

# modfied from aniMotum version 1.2-15
#' Error multiplication factors
#'
#' A function to generate a data frame of error multiplication factors (EMF) for different location classes, which can be used to account for measurement error for observations that lack error information but have a known location quality class. The default values are based on the EMF values for Argos satellite telemetry data, but users can specify their own EMF values for different location classes as needed. It is a modified version of the \code{emf} function from the \href{https://ianjonsen.github.io/aniMotum/}{aniMotum} package.
#'
#' @param gps A numeric value or a vector of length 2 specifying the error multiplication factor for GPS locations. If a single value is provided, it will be used for both x and y axes. Default is 0.1 (i.e. GPS errors are 10x more accurate than Argos \code{lc} 3.
#' @param emf.x A numeric vector of length 6 specifying the error multiplication factors for the x-axis for each location class (in order: 3, 2, 1, 0, A, B, where Z is assumed equal to B). Default values are based on the EMF values for Argos satellite telemetry data.
#' @param emf.y A numeric vector of length 6 specifying the error multiplication factors for the y-axis for each location class (in order: 3, 2, 1, 0, A, B, where Z is assumed equal to B). Default values are based on the EMF values for Argos satellite telemetry data.
#' @param prob A numeric vector of length 8 specifying the sampling probability for each location quality class (in order: G, 3, 2, 1, 0, A, B, Z). These probabilities are used internally by \code{\link{simLangevin}} and \code{\link{addMeasurementError}} to simulate observed location quality class data. Default values represent a general approximation of typical marine telemetry Argos distributions, with GPS ("G") and invalid ("Z") locations set to 0. \code{probs} must sum to 1.
#' @return A data frame with columns \code{lc}, \code{emf.x}, \code{emf.y}, and \code{prob} containing the error multiplication factors and observation probabilities for each location class. The location classes included are "G" for GPS and "3", "2", "1", "0", "A", "B", and "Z" for Argos satellite telemetry data.
#' @export
getEMF <- function (gps = 0.1,
                    emf.x = c(1, 1.54, 3.72, 13.51, 23.9, 44.22),
                    emf.y = c(1, 1.29, 2.55, 14.99, 22, 32.53),
                    prob = c(0, 0.05, 0.05, 0.10, 0.15, 0.25, 0.40, 0))
{
  if (!length(gps) %in% 1:2)
    stop("GPS emf must be a vector of length 1 or 2")
  if (length(emf.x) != 6)
    stop("Argos emf.x must be a vector of length 6")
  if (length(emf.y) != 6)
    stop("Argos emf.y must be a vector of length 6")
  if (length(prob) != 8)
    stop("prob must be a numeric vector of length 8 corresponding to c('G', '3', '2', '1', '0', 'A', 'B', 'Z')")
  if (abs(sum(prob) - 1) > 1e-6)
    stop("The 'prob' vector must exactly sum to 1")

  if (length(gps) == 1)
    gps <- c(gps, gps)

  data.frame(emf.x = c(gps[1], emf.x, emf.x[6]),
             emf.y = c(gps[2], emf.y, emf.y[6]),
             lc = as.character(c("G", "3", "2", "1", "0", "A", "B", "Z")),
             prob = prob)
}

checkErrorData <- function(data, coord=c("x","y"), measurementError = NULL, knownError = TRUE){
  if(any(is.na(data[,coord[1]]) & !is.na(data[,coord[2]])) | any(!is.na(data[,coord[1]]) & is.na(data[,coord[2]]))) stop("Missing values (NA) in coordinates must be in both x and y columns.")
  if(knownError){
    if(any(is.na(data[,coord[1]]) & (!is.na(data$smaj) | !is.na(data$smin) | !is.na(data$eor) | !is.na(data$x.err) | !is.na(data$y.err)))) stop("Measurement error terms must be NA when there are missing values in the coordinates.")
    if(any(is.na(data$smaj) & (!is.na(data$smin) & !is.na(data$eor)))) stop("When using the error ellipse model, smaj, smin, and eor must all be provided or all be NA.")
    if(any(is.na(data$x.err) & !is.na(data$y.err))) stop("When using the x- and y-axis error model, x.err and y.err must both be provided or both be NA.")
    if(any((!is.na(data$smaj) & !is.na(data$smin) & !is.na(data$eor)) & (!is.na(data$x.err) | !is.na(data$y.err)))) stop("Cannot provide both error ellipse and x- and y-axis error terms.\nIf using the error ellipse, 'smaj', 'smin', and 'eor' must all be provided and 'x.err' and 'y.err' must both be NA.\nIf using the x- and y-axis error model, 'x.err' and 'y.err' must both be provided and 'smaj', 'smin', and 'eor' must all be NA.")
    if(any((!is.na(data$smaj) | !is.na(data$smin) | !is.na(data$eor)) & (!is.na(data$x.err) & !is.na(data$y.err)))) stop("Cannot provide both error ellipse and x- and y-axis error terms.\nIf using the error ellipse, 'smaj', 'smin', and 'eor' must all be provided and 'x.err' and 'y.err' must both be NA.\nIf using the x- and y-axis error model, 'x.err' and 'y.err' must both be provided and 'smaj', 'smin', and 'eor' must all be NA.")
    if(isTRUE(any(data$eor<0 | data$eor > pi, na.rm=TRUE))) stop("Error ellipse orientation (eor) must be between 0 and pi radians.")
  }
  if(!is.null(measurementError)){
    if(knownError) stop("Cannot provide 'measurementError' parameters when the data already contains measurement error information. Please provide either 'measurementError' or appropriate measurement error columns in 'data', but not both.")
    if(!is.list(measurementError)) stop("'measurementError' must be a list.")
    if(!all(c("smaj.sd", "smin.sd") %in% names(measurementError)) && !all(c("x.sd", "y.sd") %in% names(measurementError))) stop("When providing 'measurementError' parameters, you must provide either 'smaj.sd', 'smin.sd', and 'eor.lim' for the error ellipse model, or 'x.sd' and 'y.sd' for the x- and y-axis error model.\nPlease provide the appropriate parameters for your chosen error model.")
    if(all(c("smaj.sd", "smin.sd") %in% names(measurementError)) && all(c("x.sd", "y.sd") %in% names(measurementError))) stop("Cannot provide both error ellipse and x- and y-axis error parameters in 'measurementError'.\nPlease provide either 'smaj.sd', 'smin.sd', and 'eor.lim' for the error ellipse model, or 'x.sd' and 'y.sd' for the x- and y-axis error model, but not both.")
    if(!is.null(measurementError$smaj.sd)){
      if(!is.numeric(measurementError$smaj.sd) || length(measurementError$smaj.sd)>1 || measurementError$smaj.sd <= 0) stop("smaj.sd must be positive numeric of length 1")
      if(!is.numeric(measurementError$smin.sd) || length(measurementError$smin.sd)>1 || measurementError$smin.sd <= 0) stop("smin.sd must be positive numeric of length 1")
      if(!is.null(measurementError$eor.lim)){
        if(!is.numeric(measurementError$eor.lim) || length(measurementError$eor.lim)!=2 || measurementError$eor.lim[1]<0 || measurementError$eor.lim[2]>180 || measurementError$eor.lim[2]<measurementError$eor.lim[1]) stop("eor.lim must be a numeric vector of length 2, where eor.lim[2] >= eor.lim[1].")
      }
    }
  }
}

mapDuplicatedTimes <- function(dat, map, par, re) {
  if (any(dat$dt < 1.e-6)) {
    # extract current maps (or initialize 1:N if NULL)
    mu_map <- if(is.null(map$mu)) 1:length(par$mu) else as.character(map$mu)
    vel_map <- if(is.null(map$vel)) 1:length(par$vel) else as.character(map$vel)

    for(i in 2:ncol(dat$Y)) {
      # if same track and dt is near 0, map current state to previous state
      if(dat$ID[i] == dat$ID[i-1] && dat$dt[i] < 1.e-6) {

        if(dat$dt[i] > 0){
          warning("Extremely small (0 < dt < 1.e-6) time step detected. Mapping states together (i.e. no change in location or velocity) to prevent numerical instability.")
        }

        # matrix indexing in TMB is column-major:
        # col 1 = elements 1 & 2; col 2 = elements 3 & 4
        idx_curr_x <- 2 * i - 1
        idx_curr_y <- 2 * i
        idx_prev_x <- 2 * (i - 1) - 1
        idx_prev_y <- 2 * (i - 1)

        mu_map[idx_curr_x] <- mu_map[idx_prev_x]
        mu_map[idx_curr_y] <- mu_map[idx_prev_y]

        vel_map[idx_curr_x] <- vel_map[idx_prev_x]
        vel_map[idx_curr_y] <- vel_map[idx_prev_y]
      }
    }
    map$mu <- factor(mu_map, levels = unique(mu_map[!is.na(mu_map)]))
    if("vel" %in% re) map$vel <- factor(vel_map, levels = unique(vel_map[!is.na(vel_map)]))
  }
  return(map)
}

get_data_signature <- function(data, coord = c("x", "y")) {
  if (is.null(data)) return(NULL)
  list(
    nrow = nrow(data),
    # Summing coordinates is a fast, unique identifier for the specific track path
    coord_sum = round(sum(data[[coord[1]]], data[[coord[2]]], na.rm = TRUE), 4),
    date_range = as.numeric(range(data$date, na.rm = TRUE))
  )
}

get_covs_signature <- function(spatialCovs) {
  if (is.null(spatialCovs)) return(NULL)

  # 1. Define fixed, deterministic cell locations
  n_cells <- terra::ncell(spatialCovs[[1]])
  first_cell  <- 1
  center_cell <- ceiling(n_cells / 2)
  last_cell   <- n_cells

  target_cells <- c(first_cell, center_cell, last_cell)

  list(
    names = names(spatialCovs),
    extent = round(as.vector(terra::ext(spatialCovs[[1]])), 4),
    nlyr = unname(sapply(spatialCovs, terra::nlyr)),

    # extract values for all target cells across all layers
    val_check = round(as.numeric(unlist(lapply(spatialCovs, function(x) x[target_cells]))), 6)
  )
}

verify_signatures <- function(fit, data = NULL, spatialCovs = NULL) {
  if (!is.null(data) && !is.null(fit$signatures$data)) {
    coord <- if (!is.null(fit$conditions$coord)) fit$conditions$coord else c("x", "y")
    current_data_sig <- get_data_signature(data, coord)

    if (!isTRUE(all.equal(fit$signatures$data, current_data_sig, tolerance = 1e-5))) {
      stop("Safeguard triggered: the 'data' provided does not match the 'data' originally used to fit the model. Did you pass a filtered or otherwise modified dataset?")
    }
  }

  if (!is.null(spatialCovs) && !is.null(fit$signatures$covs)) {
    current_covs_sig <- get_covs_signature(spatialCovs)

    if (!isTRUE(all.equal(fit$signatures$covs, current_covs_sig, tolerance = 1e-5))) {
      stop("Safeguard triggered: the 'spatialCovs' provided do not match the covariates originally used to fit the model. Please ensure you are passing the exact same raster list used to fit the model.")
    }
  }
}

gof_tests <- function(res_df){

  message("   Calculating goodness-of-fit tests...")

  res_x <- res_df$residual.x
  res_y <- res_df$residual.y

  valid_idx <- which(!is.na(res_x) & !is.na(res_y))
  rx <- res_x[valid_idx]
  ry <- res_y[valid_idx]
  mah <- rx^2 + ry^2

  if(length(rx) > 2) {
    # Suppress warnings for ties, which can occasionally happen in large datasets
    ks_x <- stats::ks.test(rx, "pnorm", mean = 0, sd = 1)
    ks_y <- stats::ks.test(ry, "pnorm", mean = 0, sd = 1)
    ks_mah <- stats::ks.test(mah, "pchisq", df = 2)

    # Box-Ljung test for autocorrelation (lag typically defaults to log(N))
    lag_val <- max(1, floor(log(length(rx))))
    lb_x <- stats::Box.test(rx, lag = lag_val, type = "Ljung-Box")
    lb_y <- stats::Box.test(ry, lag = lag_val, type = "Ljung-Box")
    lb_mah <- stats::Box.test(mah, lag = lag_val, type = "Ljung-Box")

    tests_df <- data.frame(
      metric = c("KS_x", "KS_y", "KS_mah", "LB_x", "LB_y", "LB_mah"),
      statistic = unname(c(ks_x$statistic, ks_y$statistic, ks_mah$statistic,
                           lb_x$statistic, lb_y$statistic, lb_mah$statistic)),
      p.value = unname(c(ks_x$p.value, ks_y$p.value, ks_mah$p.value,
                         lb_x$p.value, lb_y$p.value, lb_mah$p.value)),
      stringsAsFactors = FALSE
    )
  } else {
    warning("Not enough valid residuals to calculate GOF tests.")
    return(data.frame())
  }
  return(tests_df)
}

boundsWarning <- function(fit, as_warning = TRUE) {
  if (isTRUE(fit$conditions$out_of_bounds)) {
    msg_text <- "One or more estimated locations fell outside the spatial covariate extent. The raster extent should be expanded and the model refitted."

    if (as_warning) {
      warning("MODEL FIT LIKELY INVALID: ", msg_text, call. = FALSE, immediate. = TRUE)
    } else {
      cat("\n*** WARNING: MODEL FIT LIKELY INVALID ***\n")
      cat(msg_text, "\n")
    }
  }
}

# --- Internal Validation & Processing Helpers ---

.validate_lambda <- function(lambda) {
  if (is.null(lambda)) return(invisible(NULL))

  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("'lambda' must be a single numeric value.")
  }
  if (lambda < 0) {
    stop("'lambda' must be non-negative.")
  }
  return(invisible(TRUE))
}

# --- Shared Internal Helpers for getUD & regionProb ---

.extract_langevin_beta_list <- function(fit, beta, individual) {
  if ((missing(fit) && missing(beta)) || (!missing(fit) && !missing(beta))) {
    stop("Either 'fit' or 'beta' must be provided, but not both.")
  }

  is_hierLangevin <- !missing(fit) && inherits(fit, "hierLangevin")
  hierarchical_logical <- is_hierLangevin || (!missing(fit) && !is.null(fit$conditions$hierarchical) && !isFALSE(fit$conditions$hierarchical))

  if (!missing(beta) && !is.null(individual)) {
    warning("Argument 'individual' is ignored when 'beta' is provided manually.")
  }

  beta_list <- list()
  sd_beta_vec <- NULL

  if (missing(beta)) {
    rn_nat <- rownames(fit$estimates$natural)
    rn_work <- rownames(fit$estimates$working)

    if (is_hierLangevin) {
      beta_idx_mu <- which(grepl("^mu_beta", rn_work))
      beta_idx_sd <- which(grepl("^sd_beta", rn_work))

      mu_beta_vec <- fit$estimates$working[beta_idx_mu, "Estimate"]
      sd_beta_vec <- fit$estimates$working[beta_idx_sd, "Estimate"]
      beta_names <- gsub("^mu_", "", rn_work[beta_idx_mu])
      names(mu_beta_vec) <- beta_names
      names(sd_beta_vec) <- beta_names

      beta_list[["Population"]] <- list(type = "hierLangevin_population", mu = mu_beta_vec, sd = sd_beta_vec)

      if (!is.null(individual)) {
        id_col_name <- names(fit$estimates$random[[beta_names[1]]]$est)[1]
        ind_ids <- as.character(fit$estimates$random[[beta_names[1]]]$est[[id_col_name]])

        if (length(individual) == 1 && individual == "all") {
          individual <- ind_ids
        } else {
          individual <- as.character(individual)
          missing_inds <- setdiff(individual, ind_ids)
          if (length(missing_inds) > 0) {
            stop("The following individuals were not found in the fitted model: ", paste(missing_inds, collapse = ", "))
          }
        }

        for (ind in individual) {
          b_vec <- numeric(length(beta_names))
          names(b_vec) <- beta_names
          for (j in seq_along(beta_names)) {
            df_est <- fit$estimates$random[[beta_names[j]]]$est
            b_vec[j] <- df_est[df_est[[id_col_name]] == ind, "est"]
          }
          beta_list[[paste0("ID_", ind)]] <- b_vec
        }
      }
    } else if (hierarchical_logical) {
      beta_idx <- which(grepl("^mu_beta", rn_nat))
      beta_list[["Population"]] <- fit$estimates$natural[beta_idx, "Estimate"]

      if (!is.null(individual)) {
        beta_ind_df <- fit$estimates$random$beta_ind$est
        ind_ids <- as.character(beta_ind_df[[1]])

        if (length(individual) == 1 && individual == "all") {
          individual <- ind_ids
        } else {
          individual <- as.character(individual)
          missing_inds <- setdiff(individual, ind_ids)
          if (length(missing_inds) > 0) {
            stop("The following individuals were not found in the fitted model: ", paste(missing_inds, collapse = ", "))
          }
        }

        for (i in 1:nrow(beta_ind_df)) {
          ind_id <- as.character(beta_ind_df[i, 1])
          if (ind_id %in% individual) {
            beta_list[[paste0("ID_", ind_id)]] <- as.numeric(beta_ind_df[i, -1])
          }
        }
      }
    } else {
      beta_idx <- which(grepl("^beta", rn_nat))
      beta_list[["UD"]] <- fit$estimates$natural[beta_idx, "Estimate"]
    }
  } else {
    beta_list[["UD"]] <- as.numeric(beta)
  }

  list(
    beta_list = beta_list,
    is_hierLangevin = is_hierLangevin,
    hierarchical_logical = hierarchical_logical,
    sd_beta_vec = sd_beta_vec
  )
}

.extract_langevin_cov <- function(fit, is_hierLangevin, hierarchical_logical) {
  can_calc_se <- !missing(fit) && (!is.null(fit$covariance$natural) || !is.null(fit$covariance$working))
  beta_cov <- NULL
  beta_idx_cov <- NULL

  if (can_calc_se) {
    if (is_hierLangevin) {
      param_names <- gsub("^mu_", "", rownames(fit$estimates$working)[grepl("^mu_", rownames(fit$estimates$working))])
      beta_idx_cov <- grep("^beta", param_names)
      p_total <- length(param_names)

      joint_cov_idx <- c(beta_idx_cov, p_total + beta_idx_cov)
      beta_cov <- fit$covariance$working[joint_cov_idx, joint_cov_idx, drop = FALSE]
    } else if (hierarchical_logical) {
      cov_idx <- which(grepl("^mu_beta", rownames(fit$estimates$natural)))
      beta_cov <- fit$covariance$natural[cov_idx, cov_idx, drop = FALSE]
      beta_idx_cov <- cov_idx
    } else {
      cov_idx <- which(grepl("^beta", rownames(fit$estimates$natural)))
      beta_cov <- fit$covariance$natural[cov_idx, cov_idx, drop = FALSE]
      beta_idx_cov <- cov_idx
    }

    if (any(!is.finite(beta_cov))) {
      warning("The model's covariance matrix contains NaN or infinite values for the habitat selection coefficients. Skipping uncertainty calculations.")
      can_calc_se <- FALSE
      beta_cov <- NULL
    }
  }

  list(
    beta_cov = beta_cov,
    beta_idx_cov = beta_idx_cov,
    can_calc_se = can_calc_se
  )
}

.prep_barrier_raster <- function(spatialCovs, barrier, lambda, scaleFactor) {
  if (is.null(barrier)) {
    return(list(spatialCovs = spatialCovs, mod_spatialCovs = spatialCovs, barrier_sdf = NULL))
  }

  if (!(barrier %in% names(spatialCovs))) {
    stop(sprintf("Barrier raster '%s' not found in spatialCovs.", barrier))
  }

  barrier_sdf <- spatialCovs[[barrier]]
  if (!isTRUE(attr(barrier_sdf, "barLangevin"))) {
    stop("barrier is not a 'barLangevin' object created by prepBarrier.")
  }

  spatialCovs_clean <- spatialCovs
  spatialCovs_clean[[barrier]] <- NULL

  if (is.null(lambda)) {
    stop("To plot a barrier without a fitted model ('fit'), you must manually specify 'lambda'.")
  }
  .validate_lambda(lambda)

  penalty_rast <- terra::app(barrier_sdf, fun = function(x) {
    x_work <- x / scaleFactor
    ifelse(x_work <= 0, -0.5 * lambda * (x_work^2), 0)
  })
  names(penalty_rast) <- "barrier_penalty"
  mod_spatialCovs <- c(spatialCovs_clean, penalty_rast)

  list(spatialCovs = spatialCovs_clean, mod_spatialCovs = mod_spatialCovs, barrier_sdf = barrier_sdf)
}


.sample_joint_precision_re <- function(fit, is_hierLangevin, beta_idx_cov, nSims) {
  if (is.null(fit$covariance$random$jointPrecision)) {
    warning("Joint precision matrix is missing. Skipping individual-level uncertainties.")
    return(list(has_Q = FALSE))
  }

  message("     Drawing from joint precision matrix for individual-level uncertainties...")
  Q <- fit$covariance$random$jointPrecision
  L <- Matrix::Cholesky(Q, super = TRUE)
  z <- matrix(stats::rnorm(ncol(Q) * nSims), nrow = ncol(Q), ncol = nSims)

  step <- as.matrix(Matrix::solve(L, Matrix::solve(L, z, system = "Lt"), system = "P"))

  Q_names <- if (!is.null(colnames(Q))) {
    colnames(Q)
  } else {
    rep(names(fit$tmb_setup$parList), lengths(fit$tmb_setup$parList))
  }

  if (is_hierLangevin) {
    idx_mu <- which(Q_names == "mu")[beta_idx_cov]
    idx_u <- which(Q_names == "u")

    step_mu_beta <- step[idx_mu, , drop = FALSE]
    mu_beta_mle <- fit$par[which(names(fit$par) == "mu")[beta_idx_cov]]
    mu_beta_draws <- mu_beta_mle + step_mu_beta

    u_mle <- as.vector(fit$tmb_setup$parList$u)
    step_u <- step[idx_u, , drop = FALSE]
    u_draws_all <- u_mle + step_u

    list(
      has_Q = TRUE,
      mu_beta_draws = mu_beta_draws,
      u_draws_all = u_draws_all
    )
  } else {
    idx_mu_beta <- which(Q_names == "mu_beta")
    idx_log_sd_beta <- which(Q_names == "log_sd_beta")
    idx_beta_re <- which(Q_names == "beta_re")

    step_mu_beta <- step[idx_mu_beta, , drop = FALSE]
    step_log_sd_beta <- step[idx_log_sd_beta, , drop = FALSE]
    step_beta_re <- step[idx_beta_re, , drop = FALSE]

    mu_beta_mle <- fit$tmb_setup$parList$mu_beta
    log_sd_beta_mle <- fit$tmb_setup$parList$log_sd_beta
    beta_re_mle <- as.vector(fit$tmb_setup$parList$beta_re)

    mu_beta_draws <- mu_beta_mle + step_mu_beta
    log_sd_beta_draws <- log_sd_beta_mle + step_log_sd_beta
    beta_re_draws <- beta_re_mle + step_beta_re

    list(
      has_Q = TRUE,
      mu_beta_draws = mu_beta_draws,
      log_sd_beta_draws = log_sd_beta_draws,
      beta_re_draws = beta_re_draws
    )
  }
}

.get_individual_beta_draws <- function(nm, fit, is_hierLangevin, draw_data, beta_idx_cov, nSims) {
  ind_id <- sub("^ID_", "", nm)

  if (is_hierLangevin) {
    ind_idx <- which(fit$conditions$ids == ind_id)
    n_indiv <- length(fit$conditions$ids)

    my_u_draws <- matrix(0, nrow = length(beta_idx_cov), ncol = nSims)
    for (c in seq_along(beta_idx_cov)) {
      flat_idx <- ind_idx + (beta_idx_cov[c] - 1) * n_indiv
      my_u_draws[c, ] <- draw_data$u_draws_all[flat_idx, ]
    }

    t(draw_data$mu_beta_draws + my_u_draws)
  } else {
    beta_ind_df <- fit$estimates$random$beta_ind$est
    ind_idx <- which(as.character(beta_ind_df[[1]]) == ind_id)
    n_indiv <- nrow(beta_ind_df)

    n_covs <- nrow(draw_data$mu_beta_draws)
    my_beta_re_draws <- matrix(0, nrow = n_covs, ncol = nSims)
    for (c in seq_len(n_covs)) {
      flat_idx <- ind_idx + (c - 1) * n_indiv
      my_beta_re_draws[c, ] <- draw_data$beta_re_draws[flat_idx, ]
    }

    my_beta_ind_draws <- draw_data$mu_beta_draws + exp(draw_data$log_sd_beta_draws) * my_beta_re_draws
    t(my_beta_ind_draws)
  }
}

#' Print a hierLangevin object
#'
#' @param x A \code{hierLangevin} object.
#' @param ... Ignored.
#' @export
print.hierLangevin <- function(x, ...) {
  cat("Stage II (multistage CIHM) langevinSSM fit\n")
  cat("Model:", x$conditions$model, " | Individuals:", length(x$conditions$ids), "\n")
  cat("Convergence:", x$convergence, "-", x$message, "\n")
  cat("Negative log-likelihood:", round(x$objective, 3), "\n\n")
  cat("Population-level estimates (working scale):\n")
  print(x$estimates$working, digits = 4)
  cat("\nPopulation-level estimates (natural scale):\n")
  print(x$estimates$natural, digits = 4)
  invisible(x)
}
