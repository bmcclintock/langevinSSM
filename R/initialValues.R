#' Generate initial parameter values
#'
#' This function generates initial parameter values for fitting the Langevin movement model. It uses empirical estimates from the data to provide reasonable starting values for the optimization algorithm. The user can also specify initial values for any parameters, which will override the empirical estimates.
#'
#' @param data A data frame of class \code{dataLangevin} containing the formatted tracking data. See \code{\link{formatData}}.
#' @param model Character string specifying the movement model to be fitted. Must be either "underdamped" or "overdamped". Default: "underdamped".
#' @param obs_dist The measurement error distribution: ``Gaussian'' (default), ``t'', or ``Huber'' (pseudo-Huber) for all observations, or a named character vector by error type (\code{"epar"}, \code{"sderr"}) and/or location class, e.g. \code{c(epar = "Huber", sderr = "Gaussian")}. See \code{\link{fitLangevin}}.
#' @param par A list of initial parameter values. The names of the list should be a subset of c("beta","sigma","gamma","lambda","mu","vel","psi","psi_maj","b_maj","b_min","tau","rho_o","df","k"). If a parameter is not included in the list, an empirical estimate will be used as the initial value. See Details. Default: NULL.
#' @param spatialCovs A list of \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates to be included in the model. The order of the covariates in the list should match the order of the coefficients in \code{par$beta}.
#' @param coord Character vector of length 2 specifying the column names for the coordinates in the \code{data} data frame. Default: c("x", "y").
#' @param type Character string indicating the initialization strategy for the movement parameters \code{sigma} and \code{gamma}. Can be "variogram" (default), "neutral", or "empirical". See Details.
#' @param lambda Barrier strength setting passed to \code{\link{fitLangevin}}. Only \code{lambda = "estimate"} affects the result: the barrier strength then becomes a parameter and \code{par$lambda} is included in the output (see Details). For any other value (default: \code{NULL}), \code{par$lambda} cannot be specified.
#' @return A list of initial parameter values, with names corresponding to the parameters used in the model. The list will include the following parameters:
#' \item{beta}{Numeric vector of initial values for the coefficients of the spatial covariates. Length should match the number of spatial covariates.}
#' \item{sigma}{Numeric value for the initial estimate of the diffusion (or speed) parameter}
#' \item{gamma}{Numeric value for the initial estimate of the friction parameter (only for the underdamped model).}
#' \item{lambda}{Numeric value for the initial barrier strength on the natural scale (only when \code{lambda = "estimate"}).}
#' \item{mu}{Numeric matrix of initial values for the true locations of the animals. Should have the same number of rows as the number of observations in the data and 2 columns for the x and y coordinates. If \code{mu} includes missing values (\code{NA}), these are filled in using linear interpolation.}
#' \item{vel}{Numeric matrix of initial values for the true velocities of the animals (only for the underdamped model). Should have the same number of rows as the number of observations in the data and 2 columns for the x and y velocity components.} If \code{vel} includes missing values (\code{NA}), these are filled in using linear interpolation.
#' And, if provided in \code{par}, the following observation process parameters:
#' \item{psi}{Numeric value for the scaling factor of the Argos KF error ellipse model.}
#' \item{psi_maj, b_maj, b_min}{Numeric values for the scale of the semi-major axis and the exponents of the semi-major and semi-minor axes in the error ellipse power law (see \code{\link{fitLangevin}}).}
#' \item{tau}{Numeric vector of length 2 for the scaling factors of the x and y standard deviations in the LS/GPS error model, or a 2-row matrix (or named list) with one column per location class (and/or \code{"sderr"}).}
#' \item{rho_o}{Numeric value for the correlation parameter in the LS/GPS error model. Must be >=0 and <1.}
#' \item{df}{Numeric values for the degrees of freedom parameter for the t-distributed measurement error observation model.}
#' @details
#' \strong{Movement Process Parameters:}
#' If not provided in \code{par}, the initial values for the movement process parameters are generated as follows:
#' \itemize{
#'  \item \code{beta}: Initialized as a vector of zeros with length equal to the number of spatial covariates.
#'  \item \code{sigma} and \code{gamma} (``underdamped'' model only):
#'  \itemize{
#'   \item \code{type="variogram"} (default): a measurement-error-corrected variogram estimator. The mean squared displacement between locations is computed as a function of time lag (pooled across tracks and binned on a log scale). The habitat-driven process is approximated by the same Langevin process in a quadratic potential whose stationary variance matches the error-corrected variance of the observed locations (a damped harmonic oscillator for the underdamped model; an Ornstein-Uhlenbeck process for the overdamped model). The model variogram plus a measurement error term, \eqn{\psi(e_i + e_j)} where \eqn{e_i} is the trace of the error covariance matrix of location \eqn{i} implied by the error ellipse or standard errors, is fitted by weighted least squares with \eqn{\psi} free. If no error information is available, a constant nugget is estimated instead. This removes the dominant bias of the naive estimators when the measurement error is large compared with the movement between consecutive locations. Because real habitat potentials are not quadratic, \code{gamma} tends to be overestimated (by a factor of about 2-4 in simulations), but it is typically within an order of magnitude of the maximum likelihood estimate. If the variogram fit fails (e.g., too few locations), the ``neutral'' values are used with a message.
#'   \item \code{type="neutral"}: \code{sigma} is initialized using a heuristic assuming standard diffusive scaling, \eqn{\sqrt{\text{mean}(R^2 / (2 \Delta))}}, where \eqn{R^2} is the squared displacement between consecutive non-missing observations and \eqn{\Delta} is the time step; \code{gamma} is initialized as \eqn{1 / \text{median}(\Delta)}. These ignore measurement error and can be orders of magnitude too large when measurement error dominates the step lengths.
#'   \item \code{type="empirical"} (underdamped only): \code{gamma} is initialized from the lag-1 autocorrelation of the finite-differenced velocities, and \code{sigma} by a method-of-moments estimator. Measurement error makes the velocity autocorrelation negative, in which case this reverts to ``neutral''.
#'  }
#'  \item \code{lambda} (\code{lambda = "estimate"} only): Initialized at the largest numerically stable barrier strength for the initial \code{sigma} and \code{gamma} at the median time step, \eqn{\lambda_{max}(\sigma, \gamma, \text{median}(\Delta))} (the value used by \code{lambda = "tied"}; see \code{\link{fitLangevin}}). \code{lambda} is on the natural scale (units of 1/distance^2 in the units of the coordinates), i.e. it does not depend on \code{scaleFactor}.
#'  \item \code{mu}: Initialized at the locations in \code{data} as a matrix with 2 columns corresponding to the x and y coordinates. If there are missing values (\code{NA}) in \code{data[,coord]}, these are filled in using linear interpolation separately for each track (see \code{\link[stats]{approx}}).
#'  \item \code{vel} (``underdamped'' model only): Initialized as a matrix of zeros with the same number of rows as \code{data} and 2 columns corresponding to the x and y velocity components.
#' }
#' \strong{Observation Process Parameters:}
#' If not provided in \code{par}, the initial values for the observation process parameters are set to default values as follows:
#' \itemize{
#'  \item \code{psi}: Initialized to 1, which means no scaling of the error ellipse.
#'  \item \code{tau}: Initialized to c(1, 1), which means no scaling of the x and y standard deviations.
#'  \item \code{rho_o}: Initialized to 0, which means no correlation between the x and y errors.
#'  \item \code{df}: Initialized to 5, which provides a robust, heavy-tailed default.
#' }
#' @importFrom stats median approx ave cor
#' @export
initialValues <- function(data, model=c("underdamped","overdamped"), obs_dist = c("Gaussian", "t", "Huber"), par, spatialCovs, coord = c("x","y"), type = c("variogram", "neutral", "empirical"), lambda = NULL){

  model <- match.arg(model)
  obs_dist <- .norm_obs_dist(obs_dist)
  type <- match.arg(type)
  estimate_lambda <- identical(lambda, "estimate")

  if(!inherits(data,"dataLangevin")) stop("'data' is not formatted as a 'dataLangevin' object. See ?formatData")
  if(!missing(par)){
    if(!is.list(par)) stop("par must be a list.")
    else if(!all(names(par) %in% .par_names_all)) stop("names(par) is limited to c(", paste0("'", .par_names_all, "'", collapse = ","), ")")

    if (!estimate_lambda && "lambda" %in% names(par)) {
      stop("Cannot specify par$lambda unless a barrier is specified and lambda = \"estimate\".")
    }

    has_ee <- any(!is.na(data$smaj))
    has_ls <- any(!is.na(data$x.err))

    if (!has_ee && any(c("psi", .ell_par_names) %in% names(par))) {
      stop("Cannot specify par$psi, par$psi_maj, par$b_maj or par$b_min because the data do not contain error ellipse observations ('smaj', 'smin', 'eor').")
    }
    if (!has_ls && any(c("tau", "rho_o") %in% names(par))) {
      stop("Cannot specify par$tau or par$rho_o because the data do not contain standard error observations ('x.err', 'y.err').")
    }
    if (model == "overdamped" && any(c("gamma", "vel") %in% names(par))) {
      stop("Cannot specify par$gamma or par$vel when model = 'overdamped'.")
    }
    if (!"t" %in% obs_dist && (("df" %in% names(par)))) {
      stop("Cannot specify 'df' in 'par' unless obs_dist = 't' (for at least some observations).")
    }
    if (!"Huber" %in% obs_dist && (("k" %in% names(par)))) {
      stop("Cannot specify 'k' in 'par' unless obs_dist = 'Huber' (for at least some observations).")
    }
  } else par <- list()

  # cumulative absolute time to correctly span NA gaps
  abs_time <- stats::ave(data$dt, data$id, FUN = cumsum)

  # isolate only valid, non-NA observations
  valid_idx <- which(!is.na(data[[coord[1]]]) & !is.na(data[[coord[2]]]))

  valid_id <- data$id[valid_idx]
  valid_x <- data[[coord[1]]][valid_idx]
  valid_y <- data[[coord[2]]][valid_idx]
  valid_t <- abs_time[valid_idx]

  dx <- diff(valid_x)
  dy <- diff(valid_y)
  dt_valid <- diff(valid_t)
  R_squared <- dx^2 + dy^2

  # Only keep steps within the same track with positive time differences
  idx_keep <- which((valid_id[-1] == valid_id[-length(valid_id)]) & dt_valid > 0)

  # Baseline ("neutral") heuristic guesses
  empirical_sigma <- sqrt(mean(R_squared[idx_keep] / (2 * dt_valid[idx_keep])))
  empirical_gamma <- 1 / stats::median(dt_valid[idx_keep], na.rm = TRUE)

  need_sigma <- is.null(par$sigma)
  need_gamma <- model == "underdamped" && is.null(par$gamma)

  # ---------------------------------------------------------------------
  # Measurement-error-corrected variogram estimator (type == "variogram")
  # ---------------------------------------------------------------------
  if (type == "variogram" && (need_sigma || need_gamma)) {
    vg_est <- tryCatch(variogram_init(data, model, coord), error = function(e) NULL)
    if (is.null(vg_est)) {
      stop("   Variogram-based initialization failed (e.g., too few locations).")
    } else {
      empirical_sigma <- vg_est[["sigma"]]
      if (model == "underdamped") empirical_gamma <- vg_est[["gamma"]]
    }
  }

  # ---------------------------------------------------------------------
  # Method-of-moments (MoM) refinement of empirical_sigma/empirical_gamma
  # for the underdamped model when type == "empirical".
  # ---------------------------------------------------------------------
  if (type == "empirical" && model == "underdamped" && (is.null(par$sigma) || is.null(par$gamma))) {
    empirical_sigma_success <- FALSE
    empirical_gamma_success <- FALSE

    if (length(idx_keep) > 20) {
      vx <- dx / dt_valid
      vy <- dy / dt_valid

      same_pair <- (valid_id[-1] == valid_id[-length(valid_id)]) & dt_valid > 0
      pair_idx <- which(same_pair[-length(same_pair)] & same_pair[-1])

      if (length(pair_idx) > 20) {
        v1 <- c(vx[pair_idx], vy[pair_idx])
        v2 <- c(vx[pair_idx + 1], vy[pair_idx + 1])
        rho_v <- suppressWarnings(stats::cor(v1, v2))
        dt_pair <- stats::median(dt_valid[pair_idx], na.rm = TRUE)

        if (is.finite(rho_v) && rho_v > 1e-4 && rho_v < 1 - 1e-8 &&
            is.finite(dt_pair) && dt_pair > 0) {
          gamma_mom <- -log(rho_v) / dt_pair
          empirical_gamma <- gamma_mom
          empirical_gamma_success <- TRUE

          gdt <- gamma_mom * dt_valid[idx_keep]
          g_gdt <- 2 * gdt + 4 * expm1(-gdt) - expm1(-2 * gdt)
          Kfac <- g_gdt / gamma_mom^2 + (expm1(-gdt))^2 / (2 * gamma_mom^3)
          s2_mom <- stats::median(R_squared[idx_keep] / (2 * Kfac), na.rm = TRUE)

          if (is.finite(s2_mom) && s2_mom > 0) {
            empirical_sigma <- sqrt(s2_mom)
            empirical_sigma_success <- TRUE
          }
        }
      }
    }

    failed_sigma <- is.null(par$sigma) && !empirical_sigma_success
    failed_gamma <- is.null(par$gamma) && !empirical_gamma_success

    if (failed_sigma || failed_gamma) {
      stop("   Empirical initialization failed (e.g., due to sparse data or high measurement noise).")
    }
  }

  if(is.null(par$beta)) par$beta <- rep(0,length(spatialCovs))

  if(is.null(par$sigma)) {
    par$sigma <- empirical_sigma
  }

  if(model == "underdamped"){
    if(is.null(par$gamma)){
      par$gamma <- empirical_gamma
    }
  }

  # barrier strength (natural scale): start at the stability limit of the initial sigma and gamma
  # at the median time step, i.e. the value lambda = "tied" would use
  if (estimate_lambda && is.null(par$lambda)) {
    pos_dt <- data$dt[!is.na(data$dt) & data$dt > 1.e-6]
    dt0 <- if (length(pos_dt)) stats::median(pos_dt) else 1
    par$lambda <- .stability_limit(model, par$sigma, if (model == "underdamped") par$gamma else 1, dt0)
  }

  if(model == "underdamped"){
    if(is.null(par$vel)){
      par$vel <- matrix(0,nrow(data),2)
    } else {
      if (!is.matrix(par$vel) || nrow(par$vel) != nrow(data) || ncol(par$vel) != 2) {
        stop("'par$vel' must be a matrix with the same number of rows as 'data' and 2 columns corresponding to the x and y velocity components.")
      }
      if (any(is.na(par$vel))) {
        vel_x <- par$vel[, 1]
        vel_y <- par$vel[, 2]
        for(uid in unique(data$id)) {
          trk_idx <- which(data$id == uid)
          t_num <- as.numeric(data$date[trk_idx])

          valid_x <- !is.na(vel_x[trk_idx])
          if(any(!valid_x)) {
            vel_x[trk_idx] <- stats::approx(x = t_num[valid_x],
                                            y = vel_x[trk_idx][valid_x],
                                            xout = t_num,
                                            rule = 2)$y
          }

          valid_y <- !is.na(vel_y[trk_idx])
          if(any(!valid_y)) {
            vel_y[trk_idx] <- stats::approx(x = t_num[valid_y],
                                            y = vel_y[trk_idx][valid_y],
                                            xout = t_num,
                                            rule = 2)$y
          }
        }
        par$vel <- unname(cbind(vel_x, vel_y))
      }
    }
  }

  if(is.null(par$mu)) {
    mu_x <- data[[coord[1]]]
    mu_y <- data[[coord[2]]]

    if(any(is.na(mu_x) | is.na(mu_y))) {
      for(uid in unique(data$id)) {
        trk_idx <- which(data$id == uid)
        t_num <- as.numeric(data$date[trk_idx])

        valid_x <- !is.na(mu_x[trk_idx])
        if(any(!valid_x)) {
          mu_x[trk_idx] <- stats::approx(x = t_num[valid_x],
                                         y = mu_x[trk_idx][valid_x],
                                         xout = t_num,
                                         rule = 2)$y
        }

        valid_y <- !is.na(mu_y[trk_idx])
        if(any(!valid_y)) {
          mu_y[trk_idx] <- stats::approx(x = t_num[valid_y],
                                         y = mu_y[trk_idx][valid_y],
                                         xout = t_num,
                                         rule = 2)$y
        }
      }
    }
    par$mu <- unname(cbind(mu_x, mu_y))
  } else {
    if (!is.matrix(par$mu) || nrow(par$mu) != nrow(data) || ncol(par$mu) != 2) {
      stop("'par$mu' must be a matrix with the same number of rows as 'data' and 2 columns.")
    }
    if (any(is.na(par$mu))) {
      mu_x <- par$mu[, 1]
      mu_y <- par$mu[, 2]
      for(uid in unique(data$id)) {
        trk_idx <- which(data$id == uid)
        t_num <- as.numeric(data$date[trk_idx])

        valid_x <- !is.na(mu_x[trk_idx])
        if(any(!valid_x)) {
          mu_x[trk_idx] <- stats::approx(x = t_num[valid_x],
                                         y = mu_x[trk_idx][valid_x],
                                         xout = t_num,
                                         rule = 2)$y
        }

        valid_y <- !is.na(mu_y[trk_idx])
        if(any(!valid_y)) {
          mu_y[trk_idx] <- stats::approx(x = t_num[valid_y],
                                         y = mu_y[trk_idx][valid_y],
                                         xout = t_num,
                                         rule = 2)$y
        }
      }
      par$mu <- unname(cbind(mu_x, mu_y))
    }
  }
  return(par)
}

# ---------------------------------------------------------------------------
# Internal helpers for the measurement-error-corrected variogram initializer
# ---------------------------------------------------------------------------

# Trace of the measurement error covariance matrix for each observation, using
# the same parameterization as the TMB template (psi = tau = 1, rho_o = 0):
#   Argos KF ellipse:  (smaj^2 + smin^2) / 2
#   LS/GPS:            x.err^2 + y.err^2
# Locations without error information get 0.
#' @noRd
obs_error_trace <- function(data) {
  e <- rep(0, nrow(data))
  has_col <- function(v) !is.null(data[[v]])
  kf <- if (all(vapply(c("smaj", "smin", "eor"), has_col, logical(1))))
    !is.na(data$smaj) & !is.na(data$smin) & !is.na(data$eor) else rep(FALSE, nrow(data))
  ls <- if (all(vapply(c("x.err", "y.err"), has_col, logical(1))))
    !kf & !is.na(data$x.err) & !is.na(data$y.err) else rep(FALSE, nrow(data))
  e[kf] <- (data$smaj[kf]^2 + data$smin[kf]^2) / 2
  e[ls] <- data$x.err[ls]^2 + data$y.err[ls]^2
  e
}

# Binned empirical variogram (mean squared displacement summed over both
# coordinates) as a function of time lag, pooled across tracks. For each pair of
# locations the expected measurement error contribution e_i + e_j is also stored.
#' @noRd
empirical_variogram <- function(data, coord, nlag = 40, maxfrac = 0.25, minpairs = 10) {
  t_abs <- stats::ave(data$dt, data$id, FUN = cumsum)
  e <- obs_error_trace(data)
  ok <- !is.na(data[[coord[1]]]) & !is.na(data[[coord[2]]])
  out <- list()
  for (u in unique(data$id[ok])) {
    i <- which(ok & data$id == u)
    n <- length(i)
    if (n < 3) next
    lags <- unique(round(exp(seq(0, log(max(1, floor(n * maxfrac))), length.out = nlag))))
    for (k in lags) {
      a <- i[1:(n - k)]
      b <- i[(1 + k):n]
      out[[length(out) + 1]] <- cbind(
        tau = t_abs[b] - t_abs[a],
        r2  = (data[[coord[1]]][b] - data[[coord[1]]][a])^2 + (data[[coord[2]]][b] - data[[coord[2]]][a])^2,
        e   = e[a] + e[b])
    }
  }
  if (!length(out)) return(NULL)
  P <- as.data.frame(do.call(rbind, out))
  P <- P[P$tau > 0 & is.finite(P$r2), ]
  if (nrow(P) < 5 * minpairs) return(NULL)
  br <- unique(stats::quantile(log(P$tau), seq(0, 1, length.out = nlag + 1)))
  grp <- cut(log(P$tau), br, include.lowest = TRUE)
  vg <- do.call(rbind, lapply(split(P, grp), function(z) {
    if (nrow(z) >= minpairs)
      data.frame(tau = mean(z$tau), r2 = mean(z$r2), e = mean(z$e),
                 se = max(stats::sd(z$r2) / sqrt(nrow(z)), 1e-12))
  }))
  if (is.null(vg) || nrow(vg) < 5) return(NULL)
  vg
}

# Autocorrelation of position for a damped harmonic oscillator (the underdamped
# Langevin model in a quadratic potential) with friction g and spring constant k
#' @noRd
dho_acf <- function(tau, g, k) {
  disc <- g^2 / 4 - k
  if (disc > 1e-10 * g^2) {
    r <- sqrt(disc); l1 <- g / 2 + r; l2 <- g / 2 - r
    (l1 * exp(-l2 * tau) - l2 * exp(-l1 * tau)) / (l1 - l2)
  } else if (disc < -1e-10 * g^2) {
    w <- sqrt(-disc)
    exp(-g * tau / 2) * (cos(w * tau) + g / (2 * w) * sin(w * tau))
  } else {
    exp(-g * tau / 2) * (1 + g * tau / 2)
  }
}

# Measurement-error-corrected variogram estimates of sigma (and gamma).
#
# The habitat-driven Langevin process is approximated by the same process in a
# quadratic potential whose stationary variance matches the (error-corrected)
# variance of the observed locations. For the underdamped model this is a damped
# harmonic oscillator with stationary velocity variance sigma^2 per coordinate, so
#   E|X(t+tau) - X(t)|^2 = 4 (sigma^2 / k) (1 - acf(tau; gamma, k)),  k = 2 sigma^2 / s2,
# and for the overdamped model an OU process with
#   E|X(t+tau) - X(t)|^2 = 2 s2 (1 - exp(-sigma^2 tau / s2)),
# where s2 is the summed (x + y) location variance. The observed mean squared
# displacement is modelled as the above plus psi * (e_i + e_j), where the measurement
# error contribution comes from the error ellipses / standard errors and psi is a free
# multiplier (absorbing, e.g., the larger variance of t-distributed errors). If there
# is no error information, a constant nugget is estimated instead. The fit is weighted
# least squares on the binned variogram.
#' @noRd
variogram_init <- function(data, model, coord) {
  data <- as.data.frame(data)
  vg <- empirical_variogram(data, coord)
  if (is.null(vg)) return(NULL)

  # error-corrected stationary variance of locations (x + y), weighted across tracks
  e <- obs_error_trace(data)
  ok <- !is.na(data[[coord[1]]]) & !is.na(data[[coord[2]]])
  idx <- split(which(ok), data$id[ok])
  idx <- idx[vapply(idx, length, integer(1)) > 2]
  if (!length(idx)) return(NULL)
  v <- vapply(idx, function(i) stats::var(data[[coord[1]]][i]) + stats::var(data[[coord[2]]][i]) - mean(e[i]), numeric(1))
  n_i <- vapply(idx, length, integer(1))
  s2 <- sum(v * n_i) / sum(n_i)
  if (!is.finite(s2) || s2 <= 0) s2 <- max(vg$r2)

  use_err <- any(vg$e > 0)
  w <- 1 / vg$se^2
  tau_min <- min(vg$tau); tau_max <- max(vg$tau)

  nugget <- function(p) if (use_err) exp(p) * vg$e else exp(p)
  nug0 <- if (use_err) 0 else log(max(vg$r2[1] / 2, 1e-12))
  sig0 <- sqrt(max(vg$r2[nrow(vg)] - if (use_err) vg$e[nrow(vg)] else 0, vg$r2[nrow(vg)] / 10) / tau_max)  # rough L/T scale

  wls <- function(m) { if (any(!is.finite(m))) return(1e300); sum(w * (vg$r2 - m)^2) }
  best <- NULL
  try_fit <- function(p0, f) {
    o <- try(stats::optim(p0, f, control = list(maxit = 5000, reltol = 1e-12)), silent = TRUE)
    if (!inherits(o, "try-error") && is.finite(o$value) && (is.null(best) || o$value < best$value)) best <<- o
  }

  if (model == "underdamped") {
    f <- function(p) {
      s2v <- exp(2 * p[1]); g <- exp(p[2]); k <- 2 * s2v / s2
      wls(4 * (s2v / k) * (1 - dho_acf(vg$tau, g, k)) + nugget(p[3]))
    }
    for (lg in log(c(0.01, 0.1, 1, 10) / tau_min))
      for (ls in log(sig0 * c(0.1, 1, 10)))
        try_fit(c(ls, lg, nug0), f)
    if (is.null(best)) return(NULL)
    sigma <- exp(best$par[1]); gamma <- exp(best$par[2])
    # keep gamma within the range the data can inform
    gamma <- min(max(gamma, 0.1 / tau_max), 10 / tau_min)
    if (!is.finite(sigma) || sigma <= 0) return(NULL)
    return(c(sigma = sigma, gamma = gamma))
  } else {
    f <- function(p) {
      s2d <- exp(2 * p[1])
      wls(2 * s2 * (-expm1(-s2d * vg$tau / s2)) + nugget(p[2]))
    }
    for (ls in log(sqrt(s2 / tau_max) * c(0.1, 1, 10, 100)))
      try_fit(c(ls, nug0), f)
    if (is.null(best)) return(NULL)
    sigma <- exp(best$par[1])
    if (!is.finite(sigma) || sigma <= 0) return(NULL)
    return(c(sigma = sigma))
  }
}
