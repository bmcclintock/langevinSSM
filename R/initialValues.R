#' Generate initial parameter values
#'
#' This function generates initial parameter values for fitting the Langevin movement model. It uses empirical estimates from the data to provide reasonable starting values for the optimization algorithm. The user can also specify initial values for any parameters, which will override the empirical estimates.
#'
#' @param data A data frame of class \code{dataLangevin} containing the formatted tracking data. See \code{\link{formatData}}.
#' @param model Character string specifying the movement model to be fitted. Must be either "underdamped" or "overdamped". Default: "underdamped".
#' @param par A list of initial parameter values. The names of the list should be a subset of c("beta","sigma","gamma","mu","vel","psi","tau","rho_o"). If a parameter is not included in the list, an empirical estimate will be used as the initial value. See Details. Default: NULL.
#' @param spatialCovs A list of \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates to be included in the model. The order of the covariates in the list should match the order of the coefficients in \code{par$beta}.
#' @param coord Character vector of length 2 specifying the column names for the coordinates in the \code{data} data frame. Default: c("x", "y").
#' @return A list of initial parameter values, with names corresponding to the parameters used in the model. The list will include the following parameters:
#' \item{beta}{Numeric vector of initial values for the coefficients of the spatial covariates. Length should match the number of spatial covariates.}
#' \item{sigma}{Numeric value for the initial estimate of the diffusion (or speed) parameter}
#' \item{gamma}{Numeric value for the initial estimate of the friction parameter (only for the underdamped model).}
#' \item{mu}{Numeric matrix of initial values for the true locations of the animals. Should have the same number of rows as the number of observations in the data and 2 columns for the x and y coordinates.}
#' \item{vel}{Numeric matrix of initial values for the true velocities of the animals (only for the underdamped model). Should have the same number of rows as the number of observations in the data and 2 columns for the x and y velocity components.}
#' And, if provided in \code{par}, the following observation process parameters:
#' \item{psi}{Numeric value for the scaling factor of the Argos KF error ellipse model.}
#' \item{tau}{Numeric vector of length 2 for the scaling factors of the x and y standard deviations in the LS/GPS error model.}
#' \item{rho_o}{Numeric value for the correlation parameter in the LS/GPS error model. Must be >=0 and <1.}
#' @details
#' \strong{Movement Process Parameters:}
#' If not provided in \code{par}, the initial values for the movement process parameters are generated as follows:
#' \itemize{
#'  \item \code{beta}: Initialized as a vector of zeros with length equal to the number of spatial covariates.
#'  \item \code{sigma}: Initialized as a ``neutral'' empirical estimate of the diffusion (or speed) parameter, calculated as \eqn{\sqrt{\text{mean}(R^2 / (2 * \Delta))}}, where \eqn{R^2} is the squared displacement between consecutive non-missing observations and \eqn{\Delta} is the time step between those observations.
#'  \item \code{gamma} (``underdamped'' model only): Initialized as a ``neutral'' empirical estimate of the friction parameter, calculated as \eqn{1 / \text{median}(\Delta)}, where \eqn{\Delta} is the time step between consecutive non-missing observations.
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
#' @importFrom stats median approx ave
#' @export
initialValues <- function(data, model=c("underdamped","overdamped"), par, spatialCovs, coord = c("x","y")){

  model <- match.arg(model)

  if(!inherits(data,"dataLangevin")) stop("'data' is not formatted as a 'dataLangevin' object. See ?formatData")
  if(!missing(par)){
    if(!is.list(par)) stop("par must be a list.")
    else if(!all(names(par) %in% c("beta","sigma","gamma","mu","vel","psi","tau","rho_o"))) stop("names(par) is limited to c('beta','sigma','gamma','mu','vel','psi','tau','rho_o')")

    has_ee <- any(!is.na(data$smaj))
    has_ls <- any(!is.na(data$x.sd))

    if (!has_ee && "psi" %in% names(par)) {
      stop("Cannot specify par$psi because the data do not contain error ellipse observations ('smaj', 'smin', 'eor').")
    }
    if (!has_ls && any(c("tau", "rho_o") %in% names(par))) {
      stop("Cannot specify par$tau or par$rho_o because the data do not contain standard deviation observations ('x.sd', 'y.sd'.")
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

  # guess for sigma: sqrt(R_squared/(2*dt))
  dx <- diff(valid_x)
  dy <- diff(valid_y)
  dt_valid <- diff(valid_t)
  R_squared <- dx^2 + dy^2

  # guess for gamma: inverse of the median time step
  # Only keep steps within the same track with positive time differences
  idx_keep <- which((valid_id[-1] == valid_id[-length(valid_id)]) & dt_valid > 0)

  empirical_sigma <- sqrt(mean(R_squared[idx_keep] / (2 * dt_valid[idx_keep])))
  empirical_gamma <- 1 / stats::median(dt_valid[idx_keep], na.rm = TRUE)

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
        stop("'par$vel' cannot contain missing values (NAs).")
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
      stop("'par$mu' cannot contain missing values (NAs). Please interpolate missing locations.")
    }
  }
  return(par)
}
