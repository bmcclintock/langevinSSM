#' Generate initial parameter values
#'
#' This function generates initial parameter values for fitting the Langevin movement model. It uses empirical estimates from the data to provide reasonable starting values for the optimization algorithm. The user can also specify initial values for any parameters, which will override the empirical estimates.
#'
#' @param data A data frame of class \code{dataLangevin} containing the formatted tracking data. See \code{\link{formatData}}.
#' @param model Character string specifying the movement model to be fitted. Must be either "underdamped" or "overdamped". Default: "underdamped".
#' @param par A list of initial parameter values. The names of the list should be a subset of c("beta","sigma","gamma","mu","vel","psi","tau","rho_o"). If a parameter is not included in the list, an empirical estimate will be used as the initial value. See Details. Default: NULL.
#' @param spatialCovs A list of \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates to be included in the model. The order of the covariates in the list should match the order of the coefficients in \code{par$beta}.
#' @param coord Character vector of length 2 specifying the column names for the coordinates in the \code{data} data frame. Default: c("x", "y").
#' @param type Character string indicating the initialization strategy for the movement parameters \code{sigma} and \code{gamma}. Can be "neutral" (default) or "empirical". See Details.
#' @return A list of initial parameter values, with names corresponding to the parameters used in the model. The list will include the following parameters:
#' \item{beta}{Numeric vector of initial values for the coefficients of the spatial covariates. Length should match the number of spatial covariates.}
#' \item{sigma}{Numeric value for the initial estimate of the diffusion (or speed) parameter}
#' \item{gamma}{Numeric value for the initial estimate of the friction parameter (only for the underdamped model).}
#' \item{mu}{Numeric matrix of initial values for the true locations of the animals. Should have the same number of rows as the number of observations in the data and 2 columns for the x and y coordinates. If \code{mu} includes missing values (\code{NA}), these are filled in using linear interpolation.}
#' \item{vel}{Numeric matrix of initial values for the true velocities of the animals (only for the underdamped model). Should have the same number of rows as the number of observations in the data and 2 columns for the x and y velocity components.} If \code{vel} includes missing values (\code{NA}), these are filled in using linear interpolation.
#' And, if provided in \code{par}, the following observation process parameters:
#' \item{psi}{Numeric value for the scaling factor of the Argos KF error ellipse model.}
#' \item{tau}{Numeric vector of length 2 for the scaling factors of the x and y standard deviations in the LS/GPS error model.}
#' \item{rho_o}{Numeric value for the correlation parameter in the LS/GPS error model. Must be >=0 and <1.}
#' @details
#' \strong{Movement Process Parameters:}
#' If not provided in \code{par}, the initial values for the movement process parameters are generated as follows:
#' \itemize{
#'  \item \code{beta}: Initialized as a vector of zeros with length equal to the number of spatial covariates.
#'  \item \code{sigma}: Under \code{type="empirical"}, for the underdamped model, it is initialized using a Method-of-Moments estimator based on the stationary variance of the OU process. Under \code{type="neutral"}, it is initialized using a heuristic assuming standard diffusive scaling: \eqn{\sqrt{\text{mean}(R^2 / (2 * \Delta))}}, where \eqn{R^2} is the squared displacement between consecutive non-missing observations and \eqn{\Delta} is the time step.
#'  \item \code{gamma} (``underdamped'' model only): Under \code{type="empirical"}, it is initialized based on the lag-1 autocorrelation of the finite-differenced velocities. Under \code{type="neutral"}, initialized as \eqn{1 / \text{median}(\Delta)}.
#'  \item \code{mu}: Initialized at the locations in \code{data} as a matrix with 2 columns corresponding to the x and y coordinates. If there are missing values (\code{NA}) in \code{data[,coord]}, these are filled in using linear interpolation separately for each track (see \code{\link[stats]{approx}}).
#'  \item \code{vel} (``underdamped'' model only): Initialized as a matrix of zeros with the same number of rows as \code{data} and 2 columns corresponding to the x and y velocity components.
#' }
#' \strong{Observation Process Parameters:}
#' If not provided in \code{par}, the initial values for the observation process parameters are set to default values as follows:
#' \itemize{
#'  \item \code{psi}: Initialized to 1, which means no scaling of the error ellipse.
#'  \item \code{tau}: Initialized to c(1, 1), which means no scaling of the x and y standard deviations.
#'  \item \code{rho_o}: Initialized to 0, which means no correlation between the x and y errors.
#' }
#' @importFrom stats median approx ave cor
#' @export
initialValues <- function(data, model=c("underdamped","overdamped"), par, spatialCovs, coord = c("x","y"), type = c("neutral", "empirical")){

  model <- match.arg(model)
  type <- match.arg(type)

  if(!inherits(data,"dataLangevin")) stop("'data' is not formatted as a 'dataLangevin' object. See ?formatData")
  if(!missing(par)){
    if(!is.list(par)) stop("par must be a list.")
    else if(!all(names(par) %in% c("beta","sigma","gamma","mu","vel","psi","tau","rho_o"))) stop("names(par) is limited to c('beta','sigma','gamma','mu','vel','psi','tau','rho_o')")

    has_ee <- any(!is.na(data$smaj))
    has_ls <- any(!is.na(data$x.err))

    if (!has_ee && "psi" %in% names(par)) {
      stop("Cannot specify par$psi because the data do not contain error ellipse observations ('smaj', 'smin', 'eor').")
    }
    if (!has_ls && any(c("tau", "rho_o") %in% names(par))) {
      stop("Cannot specify par$tau or par$rho_o because the data do not contain standard error observations ('x.err', 'y.err').")
    }
    if (model == "overdamped" && any(c("gamma", "vel") %in% names(par))) {
      stop("Cannot specify par$gamma or par$vel when model = 'overdamped'.")
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
      message("   Empirical initialization failed (e.g., due to sparse data or high measurement noise). Reverting to 'neutral' initialization.")
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
