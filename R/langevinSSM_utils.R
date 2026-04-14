#' @importFrom utils globalVariables
utils::globalVariables(c("mu.x", "mu.y", "vel.x", "vel.y", "id", "psi", "tau", "dt", "x", "y", "smaj", "smin", "eor", "x.err", "y.err", "val", "lag", "UD", "type", "theoretical"))

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
#' @param type Character string indicating which scale to extract. Options are \code{"natural"} (default), \code{"working"}, \code{"mu"}, or \code{"vel"} (depending on the method).
#' @param parm A specification of which parameters are to be given confidence intervals.
#' @param level The confidence level required. Default: \code{0.95}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @name langevin_methods
NULL

#' @rdname langevin_methods
#' @export
logLik.fitLangevin <- function(object, ...) {
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
coef.fitLangevin <- function(object, type = "natural", ...) {

  if (!type %in% names(object$estimates)) {
    stop("Estimates for type '", type, "' are not available in this model fit.")
  }

  # Extract just the Estimate column
  est_vector <- object$estimates[[type]][, "Estimate"]

  # Assign the parameter names
  names(est_vector) <- rownames(object$estimates[[type]])

  return(est_vector)
}

#' @rdname langevin_methods
#' @export
vcov.fitLangevin <- function(object, type = "natural", ...) {

  if (!type %in% names(object$covariance)) {
    stop("Covariance matrix for type '", type, "' is not available in this model fit.")
  }

  return(object$covariance[[type]])
}

#' @rdname langevin_methods
#' @export
confint.fitLangevin <- function(object, parm, level = 0.95, type = "natural", ...) {

  # --- Smart catching for common user syntax ---
  # If a user accidentally passes 'mu' or 'vel' to 'parm', seamlessly route it to 'type'
  if (!missing(parm) && length(parm) == 1 && parm %in% c("mu", "vel")) {
    type <- parm
    is_re_parm <- TRUE
  } else {
    is_re_parm <- FALSE
  }

  # Calculate the Z multiplier based on the requested confidence level
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
  fac <- stats::qnorm(a)

  # 1. Standard Fixed Effects ("natural" or "working")
  if (type %in% c("natural", "working")) {
    cf <- coef.fitLangevin(object, type = type)
    vc <- vcov.fitLangevin(object, type = type)
    se <- sqrt(diag(vc))

    ci <- array(NA, dim = c(length(cf), 2L), dimnames = list(names(cf), pct))
    ci[, 1] <- cf + fac[1] * se
    ci[, 2] <- cf + fac[2] * se

    if (!missing(parm) && !is_re_parm) {
      missing_parms <- parm[!parm %in% rownames(ci)]
      if (length(missing_parms) > 0) {
        stop("Parameter(s) not found in model: ", paste(missing_parms, collapse = ", "))
      }
      ci <- ci[parm, , drop = FALSE]
    }

    return(ci)

    # 2. Random Effects Trajectories ("mu" or "vel")
  } else if (type %in% c("mu", "vel")) {

    if (is.null(object$estimates$random) || !type %in% names(object$estimates$random)) {
      stop("Random effect '", type, "' is not available in this model fit.")
    }

    est_df <- object$estimates$random[[type]]$est
    se_df  <- object$estimates$random[[type]]$se

    x_col <- paste0(type, ".x")
    y_col <- paste0(type, ".y")

    # Clean the percentage strings for column names (e.g., "2.5 %" -> "2.5%")
    pct_clean <- gsub(" ", "", pct)

    # Construct the compromised wide-format data frame
    ci_df <- data.frame(
      id = est_df$id,
      # Force integer conversion so tests don't fail against c(1, 2)
      time_step = as.integer(ave(as.character(est_df$id), est_df$id, FUN = seq_along))
    )

    # X-axis Estimates and Bounds
    ci_df[[x_col]] <- est_df[[x_col]]
    ci_df[[paste0(x_col, "_", pct_clean[1])]] <- est_df[[x_col]] + fac[1] * se_df[[x_col]]
    ci_df[[paste0(x_col, "_", pct_clean[2])]] <- est_df[[x_col]] + fac[2] * se_df[[x_col]]

    # Y-axis Estimates and Bounds
    ci_df[[y_col]] <- est_df[[y_col]]
    ci_df[[paste0(y_col, "_", pct_clean[1])]] <- est_df[[y_col]] + fac[1] * se_df[[y_col]]
    ci_df[[paste0(y_col, "_", pct_clean[2])]] <- est_df[[y_col]] + fac[2] * se_df[[y_col]]

    return(ci_df)

  } else {
    stop("The 'type' argument must be one of: 'natural', 'working', 'mu', or 'vel'.")
  }
}

#' @rdname langevin_methods
#' @export
residuals.fitLangevin <- function(object, ...) {

  if (is.null(object$residuals)) {
    stop("Residuals were not calculated during model fitting. Please re-run fitLangevin with 'calcResiduals = TRUE'.")
  }

  res <- object$residuals

  return(res)
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
#' @return A data frame with columns \code{lc}, \code{emf.x}, and \code{emf.y} containing the error multiplication factors for each location class. The location classes included are "G" for GPS and "3", "2", "1", "0", "A", "B", and "Z" for Argos satellite telemetry data.
#' @export
getEMF <- function (gps = 0.1, emf.x = c(1, 1.54, 3.72, 13.51, 23.9, 44.22),
                     emf.y = c(1, 1.29, 2.55, 14.99, 22, 32.53))
{
  if (!length(gps) %in% 1:2)
    stop("GPS emf must be a vector of length 1 or 2")
  if (length(emf.x) != 6)
    stop("Argos emf.x must be a vector of length 6")
  if (length(emf.y) != 6)
    stop("Argos emf.y must be a vector of length 6")
  if (length(gps) == 1)
    gps <- c(gps, gps)
  data.frame(emf.x = c(gps[1], emf.x, emf.x[6]), emf.y = c(gps[2],
                                                           emf.y, emf.y[6]), lc = as.character(c("G", "3", "2",
                                                                                                 "1", "0", "A", "B", "Z")))
}

checkErrorData <- function(data, coord=c("x","y"), measurementError = NULL, knownError = TRUE){
  if(any(is.na(data[,coord[1]]) & !is.na(data[,coord[2]])) | any(!is.na(data[,coord[1]]) & is.na(data[,coord[2]]))) stop("Missing values (NA) in coordinates must be in both x and y columns.")
  if(knownError){
    if(any(is.na(data[,coord[1]]) & (!is.na(data$smaj) | !is.na(data$smin) | !is.na(data$eor) | !is.na(data$x.err) | !is.na(data$y.err)))) stop("Measurement error terms must be NA when there are missing values in the coordinates.")
    if(any(is.na(data$smaj) & (!is.na(data$smin) & !is.na(data$eor)))) stop("When using the error ellipse model, smaj, smin, and eor must all be provided or all be NA.")
    if(any(is.na(data$x.err) & !is.na(data$y.err))) stop("When using the x- and y-axis error model, x.err and y.err must both be provided or both be NA.")
    if(any((!is.na(data$smaj) & !is.na(data$smin) & !is.na(data$eor)) & (!is.na(data$x.err) | !is.na(data$y.err)))) stop("Cannot provide both error ellipse and x- and y-axis error terms.\nIf using the error ellipse, 'smaj', 'smin', and 'eor' must all be provided and 'x.err' and 'y.err' must both be NA.\nIf using the x- and y-axis error model, 'x.err' and 'y.err' must both be provided and 'smaj', 'smin', and 'eor' must all be NA.")
    if(any((!is.na(data$smaj) | !is.na(data$smin) | !is.na(data$eor)) & (!is.na(data$x.err) & !is.na(data$y.err)))) stop("Cannot provide both error ellipse and x- and y-axis error terms.\nIf using the error ellipse, 'smaj', 'smin', and 'eor' must all be provided and 'x.err' and 'y.err' must both be NA.\nIf using the x- and y-axis error model, 'x.err' and 'y.err' must both be provided and 'smaj', 'smin', and 'eor' must all be NA.")
    if(isTRUE(any(data$eor<0 | data$eor > pi))) stop("Error ellipse orientation (eor) must be between 0 and pi radians.")
  }
  if(!is.null(measurementError)){
    if(knownError) stop("Cannot provide 'measurementError' parameters when the data already contains measurement error information. Please provide either 'measurementError' or appropriate measurement error columns in 'data', but not both.")
    if(!is.list(measurementError)) stop("'measurementError' must be a list.")
    if(!all(c("smaj.sd", "smin.sd") %in% names(measurementError)) && !all(c("x.sd", "y.sd") %in% names(measurementError))) stop("When providing 'measurementError' parameters, you must provide either 'smaj.sd', 'smin.sd', and 'eor.lim' for the error ellipse model, or 'x.sd' and 'y.sd' for the x- and y-axis error model.\nPlease provide the appropriate parameters for your chosen error model.")
    if(all(c("smaj.sd", "smin.sd") %in% names(measurementError)) && all(c("x.sd", "y.sd") %in% names(measurementError))) stop("Cannot provide both error ellipse and x- and y-axis error parameters in 'measurementError'.\nPlease provide either 'smaj.sd', 'smin.sd', and 'eor.lim' for the error ellipse model, or 'x.sd' and 'y.sd' for the x- and y-axis error model, but not both.")
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
    warning("Not enough valid residuals to calculate quantitative GOF tests.")
  }
  return(tests_df)
}
