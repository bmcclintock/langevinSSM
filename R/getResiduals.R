#' Calculate one-step-ahead (OSA) residuals post-fit
#'
#' Reconstructs the TMB objective function from a fitted model to calculate OSA residuals using user-specified methods, without needing to refit the model.
#'
#' @param fit A \code{fitLangevin} object returned by \code{\link{fitLangevin}}.
#' @param data The \code{dataLangevin} object originally used to fit the model.
#' @param spatialCovs The list of \code{SpatRaster} objects originally used to fit the model.
#' @param method Character string specifying the OSA method. Default is \code{"oneStepGaussianOffMode"}. See \code{\link[TMB]{oneStepPredict}}.
#' @param trace Logical; Trace progress? See \code{\link[TMB]{oneStepPredict}}. Default: \code{FALSE}.
#' @param run_tests Logical; calculate quantitative goodness-of-fit tests? (Kolmogorov-Smirnov for normality/chi-square, and Ljung-Box for autocorrelation). The results are printed to the console and attached as a data frame to the \code{"tests"} attribute of the output. Default: \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[TMB]{oneStepPredict}}.
#' @return An \code{resLangevin} data frame containing the OSA residuals. If \code{run_tests = TRUE}, the data frame will have an attribute \code{"tests"} containing a data frame of goodness-of-fit statistics and p-values.
#' @examples
#' par <- list(beta = c(-4, 6, 5, -0.1), sigma = 5, gamma = 0.5)
#' measurementError <- list(smaj.sd = 1.5, smin.sd = 0.75, eor = c(0,180))
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
#' fit$residuals <- getResiduals(fit, data = smallDat, spatialCovs = exampleCovs, run_tests = TRUE)
#'
#' @importFrom TMB MakeADFun oneStepPredict
#' @importFrom stats ks.test Box.test
#' @export
getResiduals <- function(fit, data, spatialCovs, method = "oneStepGaussianOffMode", trace = FALSE, run_tests = FALSE, ...) {

  if(!inherits(fit, "fitLangevin")) stop("'fit' must be a fitLangevin object.")
  if(is.null(fit$tmb_setup)) stop("The provided fit object does not contain a 'tmb_setup' blueprint.")

  if(!inherits(data,"dataLangevin")) stop("'data' must be a dataLangevin object (as returned by formatData or simLangevin.")

  verify_signatures(fit, data = data, spatialCovs = spatialCovs)

  cond <- fit$conditions
  time.unit <- attr(data, "time.unit")
  coord <- cond$coord

  raster_data <-  prepareRaster(spatialCovs, scaleFactor=cond$scaleFactor, time.unit=time.unit, data = data, coord = coord)

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
  dat$K <- as.matrix(data[, c("x.sd", "y.sd")] / cond$scaleFactor)
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

  obj <- try({
    TMB::MakeADFun(
      data = c(model = "langevinSSM", dat),
      par = fit$tmb_setup$parList,
      map = fit$tmb_setup$map,
      random = fit$tmb_setup$random,
      DLL = "langevinSSM_TMBExports",
      hessian = cond$hessian,
      method = cond$method,
      silent = cond$silent,
      inner.control =  cond$inner.control
    )
  }, silent = TRUE)

  if (inherits(obj, "try-error")) stop("Failed to reconstruct TMB object: ", attr(obj, "condition")$message)

  obj$fn(fit$par)

  message("   Calculating OSA residuals using method: '", method, "'...")

  res_x <- rep(NA_real_, ncol(dat$Y))
  res_y <- rep(NA_real_, ncol(dat$Y))
  unique_ids <- unique(data$id)
  all_valid_cols <- which(dat$isd == 1)

  for (uid in unique_ids) {
    track_cols <- which(dat$ID == uid)
    valid_track_cols <- intersect(track_cols, all_valid_cols)

    if (length(valid_track_cols) > 1) {
      eval_track_cols <- valid_track_cols[-1]
    } else {
      warning("Track ID ", uid, " does not have enough valid observations for OSA calculation.")
      next
    }

    track_elements <- sort(c(2 * eval_track_cols - 1, 2 * eval_track_cols))
    other_valid_cols <- setdiff(all_valid_cols, eval_track_cols)
    cond_elements <- sort(c(2 * other_valid_cols - 1, 2 * other_valid_cols))

    message("      Processing track ID: ", uid, "...")
    track_res <- tryCatch({
      TMB::oneStepPredict(
        obj = obj,
        observation.name = "Y",
        data.term.indicator = "keep",
        method = method,
        trace = trace,
        discrete = FALSE,
        subset = track_elements,
        conditional = cond_elements,
        ...
      )
    }, error = function(e) {
      warning("OSA residual calculation failed for track ID ", uid, ". Error: ", e$message)
      return(NULL)
    })

    if (!is.null(track_res)) {
      idx_x <- seq(1, nrow(track_res), by = 2)
      idx_y <- seq(2, nrow(track_res), by = 2)

      if (length(eval_track_cols) == length(idx_x)) {
        res_x[eval_track_cols] <- track_res$residual[idx_x]
        res_y[eval_track_cols] <- track_res$residual[idx_y]
      }
    }
  }

  res_df <- data.frame(
    id = data$id,
    date = data$date,
    residual.x = res_x,
    residual.y = res_y
  )

  class(res_df) <- c("resLangevin", "data.frame")

  #    --- Goodness-of-Fit Tests ---
  if(run_tests) {

    tests_df <- gof_tests(res_df)

    attr(res_df,"tests") <- tests_df

    message("\n   --- OSA Goodness-of-Fit Results ---")

    # Capture the pretty-printed table as text
    formatted_table <- utils::capture.output(print(tests_df, row.names = FALSE))

    # Collapse the lines and send as a single message
    message(paste("   ",formatted_table, collapse = "\n"))

    message("-----------------------------------")
  }

  return(res_df)
}
