#' Generate initial parameter values
#'
#' This function generates initial parameter values for fitting the Langevin movement model. It uses empirical estimates from the data to provide reasonable starting values for the optimization algorithm. The user can also specify initial values for any parameters, which will override the empirical estimates.
#'
#' @param data A data frame of class \code{dataLangevin} containing the formatted tracking data. See \code{\link{formatData}}.
#' @param model Character string specifying the movement model to be fitted. Must be either "underdamped" or "overdamped". Default: "underdamped".
#' @param par A list of initial parameter values. The names of the list should be a subset of c("beta","sigma","gamma","mu","v_mu","psi","tau","rho_o"). If a parameter is not included in the list, an empirical estimate will be used as the initial value. Default: NULL.
#' @param spatialCovs A list of \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates to be included in the model. The order of the covariates in the list should match the order of the coefficients in \code{par$beta}.
#' @param coord Character vector of length 2 specifying the column names for the coordinates in the \code{data} data frame. Default: c("x", "y").
#' @return A list of initial parameter values, with names corresponding to the parameters used in the model. The list will include the following parameters:
#' \item{beta}{Numeric vector of initial values for the coefficients of the spatial covariates. Length should match the number of spatial covariates.}
#' \item{sigma}{Numeric value for the initial estimate of the diffusion coefficient.}
#' \item{gamma}{Numeric value for the initial estimate of the friction coefficient (only for the underdamped model).}
#' \item{mu}{Numeric matrix of initial values for the true locations of the animals. Should have the same number of rows as the number of observations in the data and 2 columns for the x and y coordinates.}
#' \item{v_mu}{Numeric matrix of initial values for the true velocities of the animals (only for the underdamped model). Should have the same number of rows as the number of observations in the data and 2 columns for the x and y velocity components.}
#' And, if provided in \code{par}, the following observation process parameters:
#' \item{psi}{Numeric value for the scaling factor of the Argos KF error ellipse model.}
#' \item{tau}{Numeric vector of length 2 for the scaling factors of the x and y standard deviations in the LS/GPS error model.}
#' \item{rho_o}{Numeric value for the correlation parameter in the LS/GPS error model. Must be >=0 and <1.}
#' @export
initialValues <- function(data, model=c("underdamped","overdamped"), par, spatialCovs, coord = c("x","y")){

  model <- match.arg(model)

  if(!inherits(data,"dataLangevin")) stop("'data' is not formatted as a 'dataLangevin' object. See ?formatData")
  if(!missing(par)){
    if(!is.list(par)) stop("par must be a list.")
    else if(!all(names(par) %in% c("beta","sigma","gamma","mu","v_mu","psi","tau","rho_o"))) stop("names(par) is limited to c('beta','sigma','gamma','mu','v_mu','psi','tau','rho_o')")
  } else par <- list()

  # Calculate sequential differences
  dx <- diff(data$x)
  dy <- diff(data$y)
  R_squared <- dx^2 + dy^2

  idx <- which((data$id[-1] == data$id[-nrow(data)]) & !is.na(R_squared) & data$dt[-1] > 0)

  # guess for sigma: sqrt(R_squared/(2*dt))
  empirical_sigma <- sqrt(mean(R_squared[idx] / (2 * data$dt[-1][idx])))

  # guess for gamma: inverse of the median time step
  empirical_gamma <- 1 / stats::median(data$dt[-1][idx], na.rm = TRUE)

  if(is.null(par$beta)) par$beta <- rep(0,length(spatialCovs))

  if(is.null(par$sigma)) {
    par$sigma <- empirical_sigma
  }

  if(model == "underdamped"){
    if(is.null(par$gamma)){
      par$gamma <- empirical_gamma
    }
    if(is.null(par$v_mu)){
      par$v_mu <- matrix(0,nrow(data),2)
    }
  }

  if(is.null(par$mu)) par$mu <- cbind(data[,coord[1]],data[,coord[2]])

  return(par)
}
