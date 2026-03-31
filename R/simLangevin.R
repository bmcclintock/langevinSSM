#' Simulate trajectories from the habitat-driven Langevin diffusion
#' @param model Langevin model to simulate: "underdamped" (default) or "overdamped".
#' @param par List of parameters. For the "underdamped" model, this must include \code{beta} (a vector of length equal to the number of covariates), \code{sigma} (speed parameter), and \code{gamma} (friction parameter). For the "overdamped" model, this must include \code{beta} and \code{sigma}.
#' If \code{measurementError} is specified, optional observation process parameters include a scaling factor to account for uncertainty in the Argos error ellipse (\code{psi}), a 2-vector scaling factor for uncertainty in the x- and y-axis errors (\code{tau}), and a correlation term betwwen the x- and y-axis errors (\code{rho_o}).
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates. The covariates must be on the same spatial grid and have the same spatial extent.
#' @param nbAnimals Number of animals to simulate. Default: 1.
#' @param obsPerAnimal Number of observations to simulate per animal. Default: 500.
#' @param timeStep Time step to use for the simulation. The smaller the time step, the more accurate the approximation. Default: 0.01.
#' @param initialPosition Initial position(s) for the simulation. This can be a 2-vector providing the x- and y-coordinates of the initial position for all animals. Alternatively, initialPosition can be specified as a list of length \code{nbAnimals} with each element a 2-vector providing the x- and y-coordinates of the initial position for each individual. If \code{NULL} (default), initial positions are randomly generated within the spatial extent of the covariates, with a preference for areas of higher habitat quality (i.e., higher values of the utilization distribution).
#' @param measurementError List of specifications to add measurement error to the simulated trajectories.
#' Measurement error can be added in the form of the Argos Kalman Filter error ellipse (i.e., semi-major axis, semi-minor axis, and error ellipse orientation) or in the form of Argos least squares or GPS x- and y-axis errors.
#' For the error ellipse, the specifications are \code{M} (the standard deviation of the semi-major axis), \code{m} (the standard deviation of the semi-minor axis), and \code{c} (a 2-vector providing the range for the error ellipse orientation in degrees).
#' For x- and y-axis errors, the specifications are \code{x.sd} and \code{y.sd} (the standard deviation of the x- and y-axis errors, respectively). Optional observation process parameters that can be included in \code{par} are \code{psi} (a scaling factor to account for uncertainty in the error ellipse), \code{tau} (a 2-vector scaling factor to account for uncertainty in the x- and y-axis errors), and rho_o (a correlation term between the x- and y-axis errors).
#' Any of the observation process parameters not specified in \code{par} are set to their default values \code{psi = 1}, \code{tau = c(1, 1)}, and \code{rho_o = 0}.
#' @return A data frame of class \code{dataLangevin} containing the simulated trajectories. The data frame contains the following columns:
#' \item{id}{Animal ID}
#' \item{date}{Date of observation}
#' \item{dt}{Time step between observations}
#' \item{x}{Observed x-coordinate of the location}
#' \item{y}{Observed y-coordinate of the location}
#' \item{smaj}{Semi-major axis of the error ellipse}
#' \item{smin}{Semi-minor axis of the error ellipse}
#' \item{eor}{Error ellipse orientation}
#' \item{x.sd}{Standard deviation of the x-axis error}
#' \item{y.sd}{Standard deviation of the y-axis error}
#' \item{mu.x}{True x-coordinate of the location}
#' \item{mu.y}{True y-coordinate of the location}
#' \item{vel.x}{True x-velocity of the location (if \code{model="underdamped"})}
#' \item{vel.y}{True y-velocity of the location (if \code{model="underdamped"})}
#' @examples
#' # underdamped model with measurement error
#'
#' par <- list(beta = c(-4, 6, 5, -0.1),
#'             sigma = 5,
#'             gamma = 0.5)
#'
#' # error ellipse
#' # exampleCovs included in package; see ?exampleCovs for details
#' simDat_ee <- simLangevin(par = par,
#'                       spatialCovs = exampleCovs,
#'                       measurementError = list(M = 1.5,
#'                                               m = 0.75,
#'                                               c = c(0,180)))
#'
#' # x- and y-axis errors
#' # exampleCovs included in package; see ?exampleCovs for details
#' simDat_xy <- simLangevin(par = par,
#'                       spatialCovs = exampleCovs,
#'                       measurementError = list(x.sd = 1.5,
#'                                               y.sd = 1.5))
#' @references
#' Dupont F, McClintock BT, Fischer J-O, Marcoux M, Hussey N, Auger-Methe M. 2025. Inferring resource selection and utilization distributions from irregular and error-prone animal tracking data using the habitat-driven Langevin diffusion.
#'
#' Michelot T, Gloaguen P, Blackwell PG, Etienne M-P. 2019. The Langevin diffusion as a continuous-time model of animal movement and habitat selection. Methods Ecol Evol 10:1894–1907. doi: 10.1111/2041-210X.13275.
#'
#' Michelot T, Hanks E. 2025. Multiscale modelling of animal movement with persistent dynamics. arXiv doi: 10.48550/arXiv.2406.15195.
#'
#' @useDynLib langevinSSM
#' @export
simLangevin <- function(model = c("underdamped","overdamped"),
                        par,
                        spatialCovs,
                        nbAnimals = 1,
                        obsPerAnimal = 500,
                        timeStep = 0.01,
                        initialPosition,
                        measurementError = NULL){

  model <- match.arg(model)

  # Prepare raster data
  raster_data <- prepareRaster(spatialCovs)

  if(!is.finite(nbAnimals) || nbAnimals<1)
    stop("nbAnimals should be at least 1.")

  if(!is.finite(obsPerAnimal) || obsPerAnimal<1)
    stop("obsPerAnimal should be at least 1.")

  if(!is.finite(timeStep) || timeStep<=0)
    stop("timeStep should be greater than zero.")

  par <- checkPar(par, model, spatialCovs = spatialCovs)$par

  gamma <- exp(par$log_gamma)
  sigma <- exp(par$log_sigma)
  beta <- par$beta

  initialPosition <- getInitialPosition(nbAnimals,initialPosition,spatialCovs,beta)

  message("   Simulating ",model," Langevin model...")
  out <- simulate_langevin_cpp(
    model = ifelse(model=="underdamped",1,0),
    nbAnimals = nbAnimals,
    obsPerAnimal = obsPerAnimal,
    timeStep = timeStep,
    gamma = gamma,
    sigma = sigma,
    beta = beta,
    raster_data = raster_data,
    initialPosition = initialPosition
  )

  out <- addMeasurementError(model, out, par, measurementError)

  out$id <- as.factor(out$id)
  if(model=="underdamped") out <- out %>% dplyr::select(id,date,dt,x,y,smaj,smin,eor,x.sd,y.sd,mu.x,mu.y,vel.x,vel.y)
  else if(model=="overdamped") out <- out %>% dplyr::select(id,date,dt,x,y,smaj,smin,eor,x.sd,y.sd,mu.x,mu.y)
  #out$eor <- out$eor * 180 / pi # convert error ellipse orientation from radians to degrees
  class(out) <- append("dataLangevin",class(out))
  return(out)
}
