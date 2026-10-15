#' Calculate one-step-ahead (OSA) residuals post-fit
#'
#' Reconstructs the TMB objective function from a fitted model to calculate OSA residuals using user-specified methods, without needing to refit the model.
#'
#' @details
#' The function automatically selects the appropriate residual calculation method based on the data:
#'
#' \strong{Data with Measurement Error:}
#' Calculates standard One-Step-Ahead (OSA) observation residuals using TMB's \code{\link[TMB]{oneStepPredict}}.
#' These residuals evaluate the joint process and observation model by predicting the raw observation
#' \eqn{Y_t} conditioned on all previous observations \eqn{Y_{1:t-1}}.
#'
#' \strong{Data without Measurement Error (Known Locations):}
#' Adds true locations (\code{mu}) to the random effects vector and injects a
#' negligible pseudo-observation error (\eqn{10^{-4} \times \sigma \sqrt{\text{median}(\Delta t)}}) to the observations.
#' TMB's OSA algorithm is then used to perform one-step-ahead prediction.
#'
#' \strong{Windowing:} To keep the calculation stable and fast for long tracks, each track is split into windows of \code{window} observations.
#' Within each window, the OSA residuals are conditioned on the preceding \code{burnin} observations rather than on the entire past of the track.
#' If \code{oneStepPredict} still returns non-finite values for a window, the window is split in half and recomputed.
#'
#' @param object A \code{fitLangevin} object returned by \code{\link{fitLangevin}}.
#' @param data The \code{dataLangevin} object originally used to fit the model.
#' @param spatialCovs The list of \code{SpatRaster} objects originally used to fit the model.
#' @param method Character string specifying the OSA method for noisy data. See \code{\link[TMB]{oneStepPredict}}. If \code{NULL} (the default), the method depends on the observation distribution of the fitted model: \code{"oneStepGaussianOffMode"} for \code{obs_dist = "Gaussian"}, and \code{"oneStepGeneric"} for \code{obs_dist = "t"} or \code{"Huber"}.
#' The Gaussian methods (\code{"oneStepGaussianOffMode"}, \code{"oneStepGaussian"}, \code{"fullGaussian"}) assume that the one-step predictive distribution of each observation is Gaussian. This is not true for t-distributed measurement error: the residuals then have inflated variance and the wrong tail shape, so goodness-of-fit tests reject even a correctly specified model. A warning is issued if one of these methods is requested for a t model.
#' Method \code{"cdf"} is also supported for all three observation distributions. For \code{obs_dist = "Huber"} the one-step CDF has no closed form; it is evaluated by a fixed 24-point Gauss-Hermite quadrature over the normal variance-mixture representation (relative error of the tail probabilities below \eqn{2 \times 10^{-5}}). It is about 3 times faster than \code{"oneStepGeneric"} for Gaussian and t errors, but only about 1.2 times faster for \code{"Huber"}, but it applies the Laplace approximation to the one-step \emph{CDF} (rather than the density). This systematically inflates the absolute residuals when the predictive uncertainty of the latent location is not small compared with the measurement error, as for the \eqn{y | x} component. Use it for quick checks only, and use \code{"oneStepGeneric"} for formal goodness-of-fit testing.
#' @param trace Logical; Trace progress? See \code{\link[TMB]{oneStepPredict}}. Default: \code{FALSE}.
#' @param run_tests Logical; calculate quantitative goodness-of-fit tests? The results are attached as a data frame to the \code{"tests"} attribute of the output, with columns \code{metric}, \code{statistic}, \code{p.value} and \code{effect}. Default: \code{TRUE}. The tests are:
#' \itemize{
#'   \item \code{KS_x}, \code{KS_y}: Kolmogorov-Smirnov tests of the x and y residuals against N(0,1). \code{effect} is the KS distance (the largest difference between the empirical and the normal CDF).
#'   \item \code{KS_mah_shape}: scale-free test of the shape of the squared Mahalanobis distance \eqn{D^2 = r_x^2 + r_y^2}. \eqn{D^2} is rescaled by its mean and compared with the exponential distribution (\eqn{D^2/2 \sim} Exp(1) under the model); the p-value is obtained by Monte Carlo (Lilliefors-type, \code{nsim = 1000} with a fixed internal seed; the user's random number stream is not affected). \code{effect} is the KS distance.
#'   \item \code{Disp_mah}: dispersion test. The statistic is \eqn{\overline{D^2}/2}, the average residual variance (1 under the model); the two-sided p-value uses \eqn{\sum D^2 \sim \chi^2_{2n}}. \code{effect} is the relative variance error \eqn{\overline{D^2}/2 - 1}.
#'   \item \code{LB_x}, \code{LB_y}, \code{LB_mah}: Ljung-Box tests for autocorrelation of the x residuals, the y residuals and \eqn{D^2}. They are computed within each track (lag \eqn{\lfloor \log n_i \rfloor} for a track with \eqn{n_i} residuals) and pooled by summing the Ljung-Box statistics and their degrees of freedom, so the end of one track is never treated as the predecessor of the start of the next. \code{effect} is the largest absolute autocorrelation over the tested lags (averaged over tracks, weighted by \eqn{n_i}).
#' }
#' Together, \code{KS_mah_shape} and \code{Disp_mah} replace a single KS test of \eqn{D^2} against \eqn{\chi^2_2}, which mixes scale and shape errors. \code{Disp_mah} rejecting while \code{KS_mah_shape} does not indicates residuals that are correctly shaped but too small (\code{statistic < 1}) or too large (\code{statistic > 1}), typically because of mis-estimated process parameters rather than a misspecified distribution. The Mahalanobis-based tests are omitted for \code{obs_dist = "t"} and \code{"Huber"}.
#'
#' With long tracks, these tests can detect departures far too small to matter. In particular, the OSA residuals are only approximately independent (Laplace approximation, windowing, estimated parameters), so with tens of thousands of residuals the Ljung-Box tests can reject on autocorrelations of about 0.02 even when the model is correct. Judge the \code{effect} column, not just the p-value: for example, autocorrelations below about 0.05 or a variance error of a few percent are rarely of practical concern. Up to seven tests are reported, so a single p-value below 0.05 is expected fairly often under a correct model.
#' @param ncores Integer; Number of cores to use for parallel processing of independent tracks. Default is \code{1}.
#' @param window Integer; maximum number of observations whose OSA residuals are calculated in a single call to \code{\link[TMB]{oneStepPredict}}. Long tracks are processed in consecutive windows of this size. Short windows keep the unobserved tail of latent states short, which stops the inner optimisation from diverging (which otherwise gives \code{NaN} residuals). \code{window} mainly affects run time. Within a window, the not-yet-observed future latent states are integrated out by the Laplace approximation, which is not exact when the drift is nonlinear. This adds a small amount of noise to the residuals (typically about 0.003, at most about 0.07, with no bias) that levels off after about 5 future observations. Shorter windows therefore never lose accuracy. Each window pays the cost of its \code{burnin} observations once, which favours longer windows. But every OSA step re-optimises all the latent states in the window, which favours shorter ones. The fastest window is roughly \eqn{\sqrt{10 \times burnin}} for \code{"oneStepGaussianOffMode"}, and 25--50 for \code{"oneStepGeneric"}. Default: \code{50}.
#' @param burnin Integer; number of valid observations immediately before each window that are conditioned on (treated as known) when calculating the OSA residuals within the window. This approximates conditioning on all past observations. The approximation error decays quickly once \code{burnin} exceeds the velocity autocorrelation time \eqn{1/(\gamma \, \text{median}(\Delta t))} (in observations): about \eqn{10^{-5}} at 1 autocorrelation time and below \eqn{10^{-10}} at 2. Default (\code{NULL}): 1.5 autocorrelation times (\eqn{1.5/(\gamma \, \text{median}(\Delta t))}) for the underdamped model, constrained to [20, 500], and 20 for the overdamped model. The overdamped model has no velocity state, so the effect of the past dies out within about 10 observations.
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
#'                    control = list(trace = 1))
#'
#' res <- residuals(fit, data = smallDat, spatialCovs = exampleCovs)
#' print(res)
#'
#' @importFrom TMB MakeADFun oneStepPredict
#' @importFrom stats ks.test Box.test residuals median
#' @export
residuals.fitLangevin <- function(object, data, spatialCovs, method = NULL, trace = FALSE, run_tests = TRUE, ncores = 1, window = 50, burnin = NULL, ...) {

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
  scaleFactor <- cond$scaleFactor

  dat <- .fit_tmb_data(object, data, spatialCovs)

  # Check for valid measurement error observations
  valid_obs_counts <- sapply(unique(dat$ID), function(id) sum(dat$isd[dat$ID == id] == 1))
  has_measurement_error <- max(valid_obs_counts) > 1

  tmb_pars <- object$tmb_setup$parList
  tmb_map  <- object$tmb_setup$map
  par_to_pass <- object$par
  random_effs <- object$tmb_setup$random

  if (!has_measurement_error) {
    message("   No measurement error detected. Restoring latent states for process OSA residuals...")

    # Force known locations to be treated as valid observations
    valid_locs <- !is.na(dat$Y[1, ]) & !is.nan(dat$Y[1, ])
    dat$isd[valid_locs] <- 1L
    dat$obs_mod[valid_locs] <- 0L

    dat$K <- matrix(1.0, nrow = length(dat$ID), ncol = 2)

    # Scale pseudo-error to typical step's process standard deviation
    sigma_est <- unname(exp(object$par["log_sigma"]))
    median_dt <- median(dat$delta_t[dat$delta_t > 0], na.rm = TRUE)
    typical_process_sd <- sigma_est * sqrt(median_dt)
    pseudo_err <- typical_process_sd * 1e-4

    tmb_pars$l_tau <- c(log(pseudo_err), log(pseudo_err))
    tmb_pars$l_rho_o <- 0

    if (is.null(tmb_map)) tmb_map <- list()
    tmb_map$l_tau <- factor(c(NA, NA))
    tmb_map$l_rho_o <- factor(NA)

    # Unlock mu and vel to act as dynamic random effects for forward filtering
    tmb_map$mu <- NULL
    if (!"mu" %in% random_effs) {
      random_effs <- c(random_effs, "mu")
    }

    if (cond$model == "underdamped") {
      tmb_map$vel <- NULL
      if (!"vel" %in% random_effs) {
        random_effs <- c(random_effs, "vel")
      }
    }

    par_to_pass <- object$par[!(names(object$par) %in% c("l_tau", "l_rho_o"))]
  }
  # ---- choose / check the OSA method ----
  is_t <- !is.null(cond$obs_dist) && cond$obs_dist %in% c("t", "Huber") && has_measurement_error
  gaussian_methods <- c("oneStepGaussianOffMode", "oneStepGaussian", "fullGaussian")
  if (is.null(method)) {
    method <- if (is_t) "oneStepGeneric" else "oneStepGaussianOffMode"
  } else {
    method <- match.arg(method, c(gaussian_methods, "oneStepGeneric", "cdf"))
    if (is_t && method %in% gaussian_methods)
      warning("method = '", method, "' assumes a Gaussian one-step predictive distribution, but the model ",
              "has non-Gaussian measurement error (obs_dist = '", cond$obs_dist, "'). The residuals will not be N(0,1) even ",
              "if the model is correct (inflated variance, wrong tails), and goodness-of-fit tests will be ",
              "unreliable. Use method = 'oneStepGeneric' or 'cdf' instead.", call. = FALSE)
  }

  # =========================================================================
  # ONE-STEP-AHEAD RESIDUAL CALCULATION (SUBSETTED BY TRACK AND WINDOW)
  # =========================================================================
  #
  # Why windows?  oneStepGaussianOffMode (TMB default reverse = TRUE) starts with
  # every observation "on" and then switches them off one at a time, from the last
  # observation backwards, warm-starting each inner (Laplace) optimisation at the
  # previous solution (random.start = last.par[random]). After k steps, the last k
  # latent states (mu, vel) are no longer tied to any data. For long tracks
  # (thousands of locations) that unobserved tail is very long. Its conditional mode
  # under the nonlinear, raster-driven drift wanders off and leaves the covariate
  # extent, and eventually the inner Newton optimiser fails. TMB then returns NaN,
  # and last.par is left holding the failed Newton iterate. That iterate becomes the
  # starting value for the next step, so every remaining (earlier) observation in the
  # track also comes back NaN.
  #
  # Fix: run OSA on short overlapping windows. Within a window:
  #   * future observations beyond the window are dropped. This is exact, because
  #     unobserved future states integrate to one;
  #   * the 'burnin' valid observations just before the window are treated as
  #     'conditional'. This approximates conditioning on the full past, and the
  #     approximation is very good once burnin covers a few velocity
  #     autocorrelation times (1/gamma).
  # Every unobserved tail is then at most 'window' observations long. Each window is
  # also much smaller than the whole track, so every OSA step is far cheaper.
  # =========================================================================
  subset_cols_data <- function(dat, parList, mapList, cols, model_type) {
    track_idx <- cols

    t_dat <- dat
    t_dat$Y <- t_dat$Y[, track_idx, drop = FALSE]
    t_dat$times <- t_dat$times[track_idx]
    t_dat$delta_t <- t_dat$delta_t[track_idx]
    t_dat$skip_step <- t_dat$skip_step[track_idx]
    t_dat$smaj <- t_dat$smaj[track_idx]
    t_dat$smin <- t_dat$smin[track_idx]
    t_dat$eor <- t_dat$eor[track_idx]
    t_dat$K <- t_dat$K[track_idx, , drop = FALSE]
    t_dat$isd <- t_dat$isd[track_idx]
    t_dat$obs_mod <- t_dat$obs_mod[track_idx]
    t_dat$ID <- t_dat$ID[track_idx]
    t_dat$nbObs <- t_dat$nbObs[track_idx]

    global_cols <- track_idx - 1L

    if (t_dat$has_prior_mu == 1L) {
      c_global_mu <- t_dat$prior_idx_mu %/% 2L
      keep_mu <- c_global_mu %in% global_cols
      if (any(keep_mu)) {
        c_local_mu <- match(c_global_mu[keep_mu], global_cols) - 1L
        t_dat$prior_idx_mu <- as.integer((t_dat$prior_idx_mu[keep_mu] %% 2L) + 2L * c_local_mu)
        t_dat$prior_mean_mu_val <- t_dat$prior_mean_mu_val[keep_mu]
        t_dat$prior_sd_mu_val <- t_dat$prior_sd_mu_val[keep_mu]
      } else {
        t_dat$has_prior_mu <- 0L
        t_dat$prior_idx_mu <- integer(0)
        t_dat$prior_mean_mu_val <- numeric(0)
        t_dat$prior_sd_mu_val <- numeric(0)
      }
    }

    if (t_dat$has_prior_vel == 1L) {
      c_global_vel <- t_dat$prior_idx_vel %/% 2L
      keep_vel <- c_global_vel %in% global_cols
      if (any(keep_vel)) {
        c_local_vel <- match(c_global_vel[keep_vel], global_cols) - 1L
        t_dat$prior_idx_vel <- as.integer((t_dat$prior_idx_vel[keep_vel] %% 2L) + 2L * c_local_vel)
        t_dat$prior_mean_vel_val <- t_dat$prior_mean_vel_val[keep_vel]
        t_dat$prior_sd_vel_val <- t_dat$prior_sd_vel_val[keep_vel]
      } else {
        t_dat$has_prior_vel <- 0L
        t_dat$prior_idx_vel <- integer(0)
        t_dat$prior_mean_vel_val <- numeric(0)
        t_dat$prior_sd_vel_val <- numeric(0)
      }
    }

    t_pars <- parList
    t_pars$mu <- t_pars$mu[, track_idx, drop = FALSE]
    if (model_type == "underdamped") {
      t_pars$vel <- t_pars$vel[, track_idx, drop = FALSE]
    }

    t_map <- mapList
    if (!is.null(t_map$mu)) t_map$mu <- factor(matrix(t_map$mu, nrow = 2)[, track_idx, drop = FALSE])
    if (!is.null(t_map$vel) && model_type == "underdamped") t_map$vel <- factor(matrix(t_map$vel, nrow = 2)[, track_idx, drop = FALSE])

    list(dat = t_dat, pars = t_pars, map = t_map, track_idx = track_idx)
  }

  # ---- window settings ----
  if (!is.numeric(window) || length(window) != 1 || window < 1) stop("'window' must be a positive integer")
  window <- as.integer(window)
  if (is.null(burnin)) {
    med_dt <- median(dat$delta_t[dat$delta_t > 0], na.rm = TRUE)
    if (cond$model == "underdamped") {
      # ~5 velocity autocorrelation times (1/gamma), expressed in observations
      # (from the full parameter list: log_gamma is absent from object$par when it is fixed via map)
      gamma_est <- unname(exp(object$tmb_setup$parList$log_gamma))
      if (!length(gamma_est) || !is.finite(gamma_est)) gamma_est <- unname(exp(object$par["log_gamma"]))
      # 1.5 velocity autocorrelation times (1/gamma), expressed in observations.
      # Truncation error: ~1e-2 at 0.5, ~1e-5 at 1 and <1e-10 at 2 autocorrelation times,
      # versus ~1e-3 Laplace noise that OSA has anyway, so 1.5 leaves ample margin.
      burnin <- ceiling(1.5 / (gamma_est * med_dt))
      if (!is.finite(burnin)) burnin <- 20L
    } else {
      # no velocity state: the dependence on the past dies out within a few observations
      burnin <- 20L
    }
    burnin <- as.integer(min(max(burnin, 20L), 500L))
  }
  if (!is.numeric(burnin) || length(burnin) != 1 || burnin < 1) stop("'burnin' must be a positive integer")
  burnin <- as.integer(burnin)
  min_window <- 5L

  unique_ids <- unique(data$id)

  if (ncores > 1) {
    if (!requireNamespace("doFuture", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("Packages 'future' and 'doFuture' are required for multicore processing. Please install them.")
    } else {
      oldDoPar <- doFuture::registerDoFuture()
      on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)
      oldPlan <- future::plan(future::multisession, workers = ncores)
      on.exit(future::plan(oldPlan), add = TRUE)
      `%loop%` <- foreach::`%dopar%`
    }
  } else {
    `%loop%` <- foreach::`%do%`
  }

  message("   Calculating OSA residuals using method: '", method, "' (ncores = ", ncores,
          ", window = ", window, ", burnin = ", burnin, ")...")

  # OSA for the observations in 'eval_local' (track-local column indices), conditioning on
  # up to 'burnin' preceding valid observations; returns a 2 x length(eval_local) matrix
  osa_window <- function(track_idx, valid_cols, eval_local) {
    s <- eval_local[1]
    e <- eval_local[length(eval_local)]
    prev_valid <- valid_cols[valid_cols < s]
    cond_local <- utils::tail(prev_valid, burnin)
    w0 <- if (length(cond_local)) min(cond_local) else s
    win_local <- w0:e
    sub_info <- subset_cols_data(dat, tmb_pars, tmb_map, track_idx[win_local], cond$model)

    ev <- eval_local - w0 + 1L
    cv <- cond_local - w0 + 1L
    track_elements <- sort(c(2L * ev - 1L, 2L * ev))
    cond_elements <- if (length(cv)) sort(c(2L * cv - 1L, 2L * cv)) else NULL

    out <- matrix(NaN, 2, length(eval_local))

    obj_w <- try({
      TMB::MakeADFun(
        data = c(model = "langevinSSM", sub_info$dat),
        par = sub_info$pars,
        map = sub_info$map,
        random = random_effs,
        DLL = "langevinSSM_TMBExports",
        silent = TRUE
      )
    }, silent = TRUE)
    if (inherits(obj_w, "try-error")) return(out)

    obj_w$fn(par_to_pass)

    osa <- try({
      TMB::oneStepPredict(
        obj = obj_w,
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

    if (!inherits(osa, "try-error") && !is.null(osa)) {
      out[1, ] <- osa$residual[seq(1, nrow(osa), by = 2)]
      out[2, ] <- osa$residual[seq(2, nrow(osa), by = 2)]
    }
    out
  }

  # Recursive driver: if a window still yields non-finite residuals, split it and retry
  osa_recursive <- function(track_idx, valid_cols, eval_local) {
    r <- osa_window(track_idx, valid_cols, eval_local)
    if (all(is.finite(r)) || length(eval_local) <= min_window) return(r)
    h <- ceiling(length(eval_local) / 2)
    cbind(osa_recursive(track_idx, valid_cols, eval_local[1:h]),
          osa_recursive(track_idx, valid_cols, eval_local[(h + 1):length(eval_local)]))
  }

  results_list <- foreach::foreach(uid = unique_ids, .packages = c("TMB", "langevinSSM"), .errorhandling = "pass") %loop% {

    track_idx <- which(dat$ID == uid)
    valid_cols <- which(dat$isd[track_idx] == 1)

    if (length(valid_cols) <= 1) {
      return(list(uid = uid, error = "Not enough valid observations for OSA calculation."))
    }

    message("      Processing track ID: ", uid, "...")

    eval_all <- valid_cols[-1]
    chunks <- split(eval_all, ceiling(seq_along(eval_all) / window))

    res_x_vec <- rep(NA_real_, length(track_idx))
    res_y_vec <- rep(NA_real_, length(track_idx))

    for (ch in chunks) {
      r <- osa_recursive(track_idx, valid_cols, ch)
      res_x_vec[ch] <- r[1, ]
      res_y_vec[ch] <- r[2, ]
    }

    n_bad <- sum(!is.finite(res_x_vec[eval_all]) | !is.finite(res_y_vec[eval_all]))
    error_msg <- NULL
    if (n_bad == length(eval_all)) error_msg <- "oneStepPredict failed for all observations"

    list(uid = uid, track_idx = track_idx, res_x = res_x_vec, res_y = res_y_vec,
         n_bad = n_bad, error = error_msg)
  }

  res_x <- rep(NA_real_, ncol(dat$Y))
  res_y <- rep(NA_real_, ncol(dat$Y))

  for (res in results_list) {
    if (inherits(res, "error")) {
      warning("OSA calculation failed: ", conditionMessage(res))
    } else if (!is.null(res$error)) {
      warning("OSA calculation failed for track ID ", res$uid, ". Error: ", res$error)
    } else {
      if (isTRUE(res$n_bad > 0)) warning("OSA residuals could not be calculated for ", res$n_bad,
                                         " observation(s) of track ID ", res$uid, " (returned as NaN).")
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
  attr(res_df, "obs_dist") <- cond$obs_dist

  if(run_tests) {
    tests_df <- gof_tests(res_df)

    # Exclude Mahalanobis distance tests for the t and Huber models (D^2 is not chi^2_2 there)
    if (!is.null(cond$obs_dist) && cond$obs_dist %in% c("t", "Huber")) {
      tests_df <- tests_df[!grepl("mah", tests_df$metric, ignore.case = TRUE), ]
    }

    attr(res_df,"tests") <- tests_df
  }

  return(res_df)
}
