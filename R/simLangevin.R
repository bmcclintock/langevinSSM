#' @export
simLangevin <- function(model = c("underdamped","overdamped"), 
                        par,
                        spatialCovs, 
                        nbAnimals = 1,
                        obsPerAnimal = 500,
                        timeStep = 0.01,
                        initialPosition,
                        measurementError = NULL){
  
  match.arg(model)
  
  # Prepare raster data
  raster_data <- prepareRaster(spatialCovs)
  
  if(!is.finite(nbAnimals) || nbAnimals<1)
    stop("nbAnimals should be at least 1.")
  
  if(!is.finite(obsPerAnimal) || obsPerAnimal<1)
    stop("obsPerAnimal should be at least 1.")
  
  if(!is.finite(timeStep) || timeStep<=0)
    stop("timeStep should be greater than zero.")
  
  par <- checkPar(par, model)$par
  
  gamma <- exp(par$log_gamma)
  sigma <- exp(par$log_sigma)
  beta <- par$beta
  
  if(length(beta)!=length(spatialCovs)) stop("par$beta is of length ",length(beta),", but spatialCovs is of length ",length(spatialCovs),". They must be the same length.")
  
  initialPosition <- getInitialPosition(nbAnimals,initialPosition,spatialCovs)
  
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
  
  if(model=="underdamped") out <- out %>% dplyr::select(id,date,dt,x,y,smaj,smin,eor,sd.x,sd.y,mu.x,mu.y,v.x,v.y)
  else if(model=="overdamped") out <- out %>% dplyr::select(id,date,dt,x,y,smaj,smin,eor,sd.x,sd.y,mu.x,mu.y)
  class(out) <- append(class(out),"dataLangevin")
  return(out)
}
