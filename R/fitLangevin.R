#' @export
fitLangevin <- function(data, model = c("underdamped","overdamped"), spatialCovs, par, map=NULL, coord = c("x", "y"), scaleFactor = 1, smoothGradient = FALSE, npoints = 4, curweight = 0.5, zetaScale = 1, hessian=FALSE, silent=FALSE, inner.control=list(maxit=1000), control = list(trace=0,iter.max=1000,eval.max=1000)){
  
  if(!inherits(data,"dataLangevin")) stop("'data' is not formatted as a 'dataLangevin' object. See ?formatData")
  match.arg(model)
  
  time.unit <- attr(data,"time.unit")
  
  # Prepare raster data
  raster_data <- prepareRaster(spatialCovs,data=data,scaleFactor=scaleFactor,time.unit=time.unit)
  
  if(!all(coord %in% colnames(data))) stop("coord not found in data.")
  
  if (inherits(data$date, "POSIXt") || inherits(data$date, "Date")) {
    track_times <- as.numeric(difftime(data$date, as.POSIXct("1970-01-01 00:00:00", tz = "UTC"), units = time.unit))
  } else {
    track_times <- as.numeric(data$date)
  }
  
  if(any(!coord %in% colnames(data))) stop("'coord' not found in data.")
  dat <-  list(model=ifelse(model=="underdamped",1,0),
               Y=t(data[,coord])/scaleFactor,
               times=track_times,
               dt=data$dt)
  
  dat$M <- data$smaj / scaleFactor
  dat$m <- data$smin / scaleFactor
  dat$c <- data$eor
  dat$K <- as.matrix(data[,c("sd.x","sd.y")] / scaleFactor)
  dat$isd <- as.numeric(!is.na(dat$Y[1,]) & ((!is.na(dat$K[,1]) & !is.na(dat$K[,2])) | (!is.na(dat$M) & !is.na(dat$m) & !is.na(dat$c))))
  dat$obs_mod <- rep(NA,ncol(dat$Y))
  dat$obs_mod[dat$isd==1 & (!is.na(dat$K[,1]) & !is.na(dat$K[,2]))] <- 0
  dat$obs_mod[dat$isd==1 & (!is.na(dat$M) & !is.na(dat$m) & !is.na(dat$c))] <- 1
  dat$ID <- data$id
  dat$nbObs <- rep(1,ncol(dat$Y))
  dat$scale_factor <- scaleFactor
  dat <- c(dat,raster_data)
  dat$smoothGradient <- ifelse(smoothGradient,1,0)
  dat$weights <- c(curweight,rep((1-curweight)/npoints,npoints))
  dat$zetaScale <- zetaScale
  
  cp <- checkPar(par, model, map, dat=dat)
  par <- cp$par
  map <- cp$map
  re <- cp$re
  
  # scale parameters
  par$log_sigma <- par$log_sigma - log(scaleFactor)
  par$mu <- par$mu / scaleFactor
  par$v_mu <- par$v_mu / scaleFactor
  
  message("   Fitting ",model," Langevin model...")
  obj2 <-
    MakeADFun(
      dat,
      par,
      map = map,
      random = re,
      DLL = "fitLangevin",
      hessian = hessian,
      silent = silent,
      inner.control =  inner.control
    )
  
  start <- proc.time()
  fit <- do.call(stats::nlminb,args = list(
    start = obj2$par,
    objective = obj2$fn,
    gradient = obj2$gr,
    control=control))
  fit$elapsedTime <- proc.time()-start
  
  message("   Calculating SEs...")
  #fit$report <- obj2$report()
  sdreport <- sdreport(obj2)
  fit$estimates <- list()
  fit$estimates$natural <- summary(sdreport,"report") # get SEs 
  fit$estimates$working <- summary(sdreport,"fixed") # get SEs 
  if(length(re)) {
    ran_est <- summary(sdreport, "random")
    fit$estimates$random <- list()
    for(i in 1:length(re)){
      fit$estimates$random[[re[i]]] <- cbind(ran_est[(i-1)*2*ncol(dat$Y)+1:(2*ncol(dat$Y)),1],ran_est[(i-1)*2*ncol(dat$Y)+1:(2*ncol(dat$Y)),2]) * scaleFactor
      colnames(fit$estimates$random[[re[i]]]) <- c("Estimate","Std. Error")
    }
  }
  class(fit) <- append(class(fit),"fitLangevin")
  return(fit)
}
