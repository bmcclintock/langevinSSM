#' Calculate one-step-ahead (OSA) residuals post-fit
#'
#' Reconstructs the TMB objective function from a fitted model to calculate OSA residuals using user-specified methods, without needing to refit the model.
#'
#' @param object A \code{fitLangevin} object returned by \code{\link{fitLangevin}}.
#' @param data The \code{dataLangevin} object originally used to fit the model.
#' @param spatialCovs The list of \code{SpatRaster} objects originally used to fit the model.
#' @param method Character string specifying the OSA method. Default is \code{"oneStepGaussianOffMode"}. See \code{\link[TMB]{oneStepPredict}}.
#' @param trace Logical; Trace progress? See \code{\link[TMB]{oneStepPredict}}. Default: \code{FALSE}.
#' @param run_tests Logical; calculate quantitative goodness-of-fit tests? (Kolmogorov-Smirnov for normality/chi-square, and Ljung-Box for autocorrelation). The results are attached as a data frame to the \code{"tests"} attribute of the output. Default: \code{TRUE}.
#' @param ncores Integer; Number of cores to use for parallel processing of the independent tracks in \code{data}. Default is \code{1} (sequential).
#' @param ... Additional arguments passed to \code{\link[TMB]{oneStepPredict}}.
#' @return A \code{resLangevin} data frame containing the OSA residuals. If \code{run_tests = TRUE}, the data frame will have an attribute \code{"tests"} containing a data frame of goodness-of-fit statistics and p-values.
#' @examples
#' par <- list(beta = c(-4, 6, 5, -0.1), sigma = 5, gamma = 0.5)
#' measurementError <- list(smaj.sd = 1.5, smin.sd = 0.75, eor.lim = c(0,180))
#'
#' set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")
#' # ridiculously small dataset for example purposes
#' smallDat <- simLangevin(par = par,
#'                           spatialCovs = exampleCovs,
#'                           nbAnimals = 3,
#'                           obsPerAnimal = 50,
#'                           measurementError = measurementError)
#'
#' fit <- fitLangevin(model = "underdamped",
#'                    data = smallDat,
#'                    spatialCovs = exampleCovs,
#'                    silent = TRUE,
#'                    control = list(trace = 1))
#'
#' res <- residuals(fit, data = smallDat, spatialCovs = exampleCovs)
#' print(res)
#'
#' @importFrom TMB MakeADFun oneStepPredict
#' @importFrom stats ks.test Box.test residuals
#' @export
residuals.fitLangevin <- function(object, data, spatialCovs, method = "oneStepGaussianOffMode", trace = FALSE, run_tests = TRUE, ncores = 1, ...) {

  if(!inherits(object, "fitLangevin")) stop("'object' must be a fitLangevin object.")
  if(is.null(object$tmb_setup)) stop("The provided fit object does not contain a 'tmb_setup' blueprint.")
  if(!inherits(data,"dataLangevin")) stop("'data' must be a dataLangevin object (as returned by formatData or simLangevin.")

  if (!requireNamespace("foreach", quietly = TRUE))
    stop("Package 'foreach' is required for calculating residuals. Please install it.")

  verify_signatures(object, data = data, spatialCovs = spatialCovs)

  cond <- object$conditions
  boundsWarning(object)

  time.unit <- attr(data, "time.unit")
  coord <- cond$coord

  raster_data <- prepareRaster(spatialCovs, scaleFactor = cond$scaleFactor, time.unit = time.unit, data = data, coord = coord)

  if (inherits(data$date, "POSIXt") || inherits(data$date, "Date")) {
    track_times <- as.numeric(difftime(data$date, as.POSIXct("1970-01-01 00:00:00", tz = "UTC"), units = time.unit))
  } else {
    track_times <- as.numeric(data$date)
  }

  dat <- list(
    process_model = ifelse(cond$model == "underdamped", 1, 0),
    Y = t(data[, coord]) / cond$scaleFactor,
    times = track_times,
    dt = data$dt
  )

  dat$skip_step <- as.integer(dat$dt < 1.e-6)
  dat$smaj <- data$smaj / cond$scaleFactor
  dat$smin <- data$smin / cond$scaleFactor
  dat$eor <- data$eor
  dat$K <- as.matrix(data[, c("x.err", "y.err")] / cond$scaleFactor)
  dat$isd <- as.numeric(!is.na(dat$Y[1,]) & ((!is.na(dat$K[,1]) & !is.na(dat$K[,2])) | (!is.na(dat$smaj) & !is.na(dat$smin) & !is.na(dat$eor))))
  dat$obs_mod <- rep(NA, ncol(dat$Y))
  dat$obs_mod[dat$isd == 1 & (!is.na(dat$K[,1]) & !is.na(dat$K[,2]))] <- 0
  dat$obs_mod[dat$isd == 1 & (!is.na(dat$smaj) & !is.na(dat$smin) & !is.na(dat$eor))] <- 1
  dat$ID <- data$id
  dat$nbObs <- rep(1, ncol(dat$Y))
  dat$scale_factor <- cond$scaleFactor
  dat <- c(dat, raster_data)
  dat$smoothGradient <- ifelse(cond$smoothGradient, 1, 0)
  dat$weights <- c(cond$curweight, rep((1 - cond$curweight) / cond$npoints, cond$npoints))
  dat$zetaScale <- cond$zetaScale

  subset_track_data <- function(dat, parList, mapList, uid, model_type) {
    track_idx <- which(dat$ID == uid)

    t_dat <- dat
    t_dat$Y <- t_dat$Y[, track_idx, drop = FALSE]
    t_dat$times <- t_dat$times[track_idx]
    t_dat$dt <- t_dat$dt[track_idx]
    t_dat$skip_step <- t_dat$skip_step[track_idx]
    t_dat$smaj <- t_dat$smaj[track_idx]
    t_dat$smin <- t_dat$smin[track_idx]
    t_dat$eor <- t_dat$eor[track_idx]
    t_dat$K <- t_dat$K[track_idx, , drop = FALSE]
    t_dat$isd <- t_dat$isd[track_idx]
    t_dat$obs_mod <- t_dat$obs_mod[track_idx]
    t_dat$ID <- t_dat$ID[track_idx]
    t_dat$nbObs <- t_dat$nbObs[track_idx]

    t_pars <- parList
    t_pars$mu <- t_pars$mu[, track_idx, drop = FALSE]
    if (model_type == "underdamped") {
      t_pars$vel <- t_pars$vel[, track_idx, drop = FALSE]
    }

    t_map <- mapList
    if (!is.null(t_map$mu)) t_map$mu <- as.factor(matrix(t_map$mu, nrow = 2)[, track_idx, drop = FALSE])
    if (!is.null(t_map$vel) && model_type == "underdamped") t_map$vel <- as.factor(matrix(t_map$vel, nrow = 2)[, track_idx, drop = FALSE])

    list(dat = t_dat, pars = t_pars, map = t_map, track_idx = track_idx)
  }

  unique_ids <- unique(data$id)

  if (ncores > 1) {
    if (!requireNamespace("doFuture", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("Packages 'future' and 'doFuture' are required for multicore processing. Please install them.")
    } else {
      oldDoPar <- doFuture::registerDoFuture()
      on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)
      future::plan(future::multisession, workers = ncores)
      `%loop%` <- foreach::`%dopar%`
    }
  } else {
    `%loop%` <- foreach::`%do%`
  }

  message("   Calculating OSA residuals using method: '", method, "' (ncores = ", ncores, ")...")

  results_list <- foreach::foreach(uid = unique_ids, .packages = c("TMB"), .errorhandling = "pass") %loop% {

    sub_info <- subset_track_data(dat, object$tmb_setup$parList, object$tmb_setup$map, uid, cond$model)

    valid_cols <- which(sub_info$dat$isd == 1)

    if (length(valid_cols) <= 1) {
      return(list(uid = uid, error = "Not enough valid observations for OSA calculation."))
    }

    eval_cols <- valid_cols[-1]
    cond_cols <- valid_cols[1]

    track_elements <- sort(c(2 * eval_cols - 1, 2 * eval_cols))
    cond_elements <- sort(c(2 * cond_cols - 1, 2 * cond_cols))

    obj_track <- try({
      TMB::MakeADFun(
        data = c(model = "langevinSSM", sub_info$dat),
        par = sub_info$pars,
        map = sub_info$map,
        random = object$tmb_setup$random,
        DLL = "langevinSSM_TMBExports",
        silent = TRUE
      )
    }, silent = TRUE)

    if (inherits(obj_track, "try-error")) return(list(uid = uid, error = paste("MakeADFun failed:", attr(obj_track, "condition")$message)))

    obj_track$fn(object$par)

    message("      Processing track ID: ", uid, "...")
    track_res <- try({
      TMB::oneStepPredict(
        obj = obj_track,
        observation.name = "Y",
        data.term.indicator = "keep",
        method = method,
        trace = trace,
        discrete = FALSE,
        subset = track_elements,
        conditional = cond_elements,
        ...
      )
    }, silent = TRUE)

    res_x_vec <- rep(NA_real_, length(sub_info$track_idx))
    res_y_vec <- rep(NA_real_, length(sub_info$track_idx))

    if (!inherits(track_res, "try-error") && !is.null(track_res)) {
      idx_x <- seq(1, nrow(track_res), by = 2)
      idx_y <- seq(2, nrow(track_res), by = 2)
      res_x_vec[eval_cols] <- track_res$residual[idx_x]
      res_y_vec[eval_cols] <- track_res$residual[idx_y]
      error_msg <- NULL
    } else {
      error_msg <- if(inherits(track_res, "try-error")) attr(track_res, "condition")$message else "oneStepPredict returned NULL"
    }

    list(uid = uid, track_idx = sub_info$track_idx, res_x = res_x_vec, res_y = res_y_vec, error = error_msg)
  }

  res_x <- rep(NA_real_, ncol(dat$Y))
  res_y <- rep(NA_real_, ncol(dat$Y))

  for (res in results_list) {
    if (!is.null(res$error)) {
      warning("OSA calculation failed for track ID ", res$uid, ". Error: ", res$error)
    } else {
      res_x[res$track_idx] <- res$res_x
      res_y[res$track_idx] <- res$res_y
    }
  }

  res_df <- data.frame(
    id = data$id,
    date = data$date,
    residual.x = res_x,
    residual.y = res_y
  )

  class(res_df) <- c("resLangevin", "data.frame")

  if(run_tests) {
    tests_df <- gof_tests(res_df)
    attr(res_df,"tests") <- tests_df
  }

  return(res_df)
}
