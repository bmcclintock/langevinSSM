options("rgdal_show_exportToProj4_warnings"="none") # suppress annoying warnings
library(tidyverse)
library(raster)
library(rasterVis)
library(viridis)
library(ggplot2)
library(future)
library(furrr)
library(doRNG)
library(RStoolbox)
library(Rcpp)
library(aniMotum)
library(remotes)
if(!requireNamespace("geoR",quietly=TRUE) || packageVersion("geoR")!="1.8.1"){
  remotes::install_version("geoR", version = "1.8-1", repos = "http://cran.us.r-project.org")
}
library(geoR)
if(!requireNamespace("RandomFieldsUtils",quietly=TRUE)){
  remotes::install_version("RandomFieldsUtils", version = "1.2-5", repos = "http://cran.us.r-project.org")
}
if(!requireNamespace("RandomFields",quietly=TRUE)){
  remotes::install_version("RandomFields", version = "3.3-14", repos = "http://cran.us.r-project.org")
}
remotes::install_github("papayoun/Rhabit@31ddf44",dependencies = TRUE) # last commit before RandomFields was removed from dependencies
library(Rhabit)
if(!requireNamespace("momentuHMM",quietly=TRUE) || as.numeric(substr(packageVersion("momentuHMM"),1,1))<2){
  remotes::install_github("bmcclintock/momentuHMM@develop",dependencies = TRUE) # requires momentuHMM version >= 2.0.0
}
library(momentuHMM)

library(TMB)
## specify, compile, and load TMB model
compile("src/fitLangevin.cpp")
dyn.load(dynlib(paste0("src/fitLangevin")))
source("R/helper_functions.R")

model <- 1 # underdamped Langevin (model = 1) or overdamped Lagevin (model = 0)

nsims <- 100 # number of simulations
nbAnimals <- 5 # number of tracks
obsPerAnimal <- 5000 # number of simulated locations per track
timeStep <- 0.01 # time scale of simulation (should be small to help prevent discretization error)

beta <- c(-4, 6, 5, -0.1) # resource selection coefficients for the spatial covariates (cov_1, cov_2, ... cov_ncov, d2c)
ncov <- length(beta) - 1 # number of spatial covariates
sigma <- 5 # speed parameter 
gamma <- 0.5 # friction parameter (smaller value -> more directional persistence); ignored unless model=1
psi <- 1 # error SD scaling parameter

## sampling rate, missing data, and measurement error
samplingRate <- 1 # for subsampling observations from true continuous-time model (e.g. if samplingRate = 2 then data are roughly thinned by 2); must be >= 1; note bias increases with samplingRate, but relationships largely preserved
propMissing <- 0 # proportion of missing observations; passed as NA observations to TMB (so corresponding true locations treated as random effects to be estimated)

M <- 1.5 # SD for semi-major error ellipse axis ~ abs(Normal(0,M))
m <- M/2 # SD for semi-minor error ellipse axis ~ abs(Normal(0,M/2))
r <- c(0,180) # range for error ellipse orientation (degrees)

## specify scale, resolution, and spatial autocorrelation for covariates
sca <- 200 # bounding box scale
lim <- c(-1, 1, -1, 1)*sca
cropExtent <- extent(lim)
resol <- 1 # cell resolution 
covRange <- c(0.1,0.5) # lower and upper bounds for covariate spatial range parameter (lower has less spatial autocorrelation)

## smooth gradient specifications (based on adjacent cells)
## slows model fitting, but can reduce bias if \Delta_t is large (i.e. samplingRate > 1) relative to timeStep (Blackwell & Matthiopoulos 2024, https://doi.org/10.1002/ecy.4457)
## can it also reduce bias if there is non-negligible measurement error?
npoints <- 0 # number of smoothing points around current cell (0 = none, 4 = diagonal, 8 = queen neighborhood);
curweight <- 1/2 # smoothing weight of current cell location; ignored if npoints=0
weights <- c(curweight,rep((1-curweight)/npoints,npoints)) # smoothing weights (must be of length 5 or 9 and sum to 1); ignored if npoints=0
zetaScale <- 1 # scale factor for smooth gradient neighborhood (>1 increases, <1 decreases neighborhood); underdamped Langevin neighborhood = zetaScale * sqrt(2*pi) * sigma; overdamped langevin neighborhood = zetaScale * sigma / sqrt(2); ignored if npoints=0


tracepar <- FALSE # TRUE = trace parameters during each iteration of optimization

langSim <- langTMB <- list()
covlist <- list()
parMat <- matrix(NA,nrow=nsims,7+model)
if(model==0){
  colnames(parMat) <- c("sigma",paste0("beta",0:(ncov+1)),"UDcor")
} else colnames(parMat) <- c("sigma",paste0("beta",0:(ncov+1)),"gamma","UDcor")

doFuture::registerDoFuture()
set.seed(1,kind="Mersenne-Twister",normal.kind = "Inversion")
for(isim in 1:nsims){
  
  cat("Simulation",isim,"\n")
  #######################
  ## Define covariates ##
  #######################
  # Generate ncov spatial covariates
  message("   Generating covariates...")
  covlist[[isim]] <- list()
  
  # Include squared distance to origin as covariate
  xgrid <- seq(lim[1], lim[2], by=resol)
  ygrid <- seq(lim[3], lim[4], by=resol)
  xygrid <- expand.grid(xgrid,ygrid)
  dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/sca
  covlist[[isim]][[4]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))
  
  for(i in 1:ncov) {
    irange <- runif(1,covRange[1],covRange[2])
    covlist[[isim]][[i]] <- rasterToRhabit(raster(geoR::grf(NULL, grid="reg", nx=2*sca+1,
                                                            ny=2*sca+1, xlims=c(-sca-0.5,sca+0.5),ylims=c(-sca-0.5,sca+0.5),
                                                            cov.pars=c(0.1, irange * sca), messages=FALSE)))
  }
  
  names(covlist[[isim]]) <- c("cov1","cov2","cov3","d2c")
  spatialCovs <- lapply(lapply(covlist[[isim]],rhabitToRaster),function(x) {proj4string(x) <- CRS("+init=epsg:3416");return(x)})
  names(spatialCovs) <- c("cov1","cov2","cov3","d2c")
  for(i in names(spatialCovs)){
    names(spatialCovs[[i]]) <- i
  }
  
  # Compute utilization distribution for states 1 and 2
  UDrhabit <- Rhabit::getUD(covariates=covlist[[isim]], beta=beta,log=TRUE)
  UD <- rhabitToRaster(UDrhabit)
  
  # Plot covariates
  ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                   legend.title = element_text(size=12), legend.text = element_text(size=12))
  c1plot <- Rhabit::plotRaster(rhabitToRaster(covlist[[isim]][[1]]), scale.name = expression(c[1])) + ggtheme
  c2plot <- Rhabit::plotRaster(rhabitToRaster(covlist[[isim]][[2]]), scale.name = expression(c[2])) + ggtheme
  c3plot <- Rhabit::plotRaster(rhabitToRaster(covlist[[isim]][[3]]), scale.name = expression(c[3])) + ggtheme
  UDplot <- Rhabit::plotRaster(UD, scale.name = expression(pi)) + ggtheme
  #UDplot
  
  raster_stack <- stack(spatialCovs)
  
  # Prepare raster data
  raster_data <- list(
    raster_vals = array(as.numeric(values(raster_stack)), 
                        dim = c(nlayers(raster_stack),
                                nrow(raster_stack),
                                ncol(raster_stack))),
    raster_coords = coordinates(raster_stack),
    raster_resolution = res(raster_stack),
    raster_extent = c(xmin(raster_stack), xmax(raster_stack),
                      ymin(raster_stack), ymax(raster_stack)),
    n_covs = nlayers(raster_stack)
  )
  
  # simulate "high resolution" tracks; this can take a while...
  message("   Simulating tracks...")
  langSim[[isim]] <- tryCatch(stop(),error=function(e) e)
  while(inherits(langSim[[isim]],"error")){
    # randomly start tracks in locations based on UD
    initPos <- matrix(sample(ncell(UD),nbAnimals,replace=FALSE,prob=exp(getValues(UD))/sum(exp(getValues(UD)))),
                      1,nbAnimals,byrow=TRUE)
    initialPosition <- t(mapply(function(x) xyFromCell(UD,initPos[,x]),1:nbAnimals,SIMPLIFY = TRUE))
    langSim[[isim]] <- tryCatch(simulate_langevin_cpp(
      model = model,
      nbAnimals = nbAnimals,
      obsPerAnimal = obsPerAnimal,
      timeStep = timeStep,
      gamma = gamma,
      sigma = sigma,
      beta = beta,
      raster_data = raster_data,
      initialPosition = initialPosition
    ),error=function(e) e)
    if(!inherits(langSim[[isim]],"error") && (any(M>0) | any(m>0))){
      langSim[[isim]] <- measurementError_rcpp(langSim[[isim]],M=M,m=m,c=r,psi=psi,model=model)
    } else {
      message("    Retrying Simulation ",isim,": ",langSim[[isim]]$message)
    }
  }
  
  scale_factor <- 1  # sca # max(abs(langSim[[isim]][,c("mu.x","mu.y")]))
  
  # Prepare raster data
  sca_raster_data <- list(
    raster_vals = array(as.numeric(values(raster_stack)), 
                        dim = c(nlayers(raster_stack),
                                nrow(raster_stack),
                                ncol(raster_stack))),
    raster_coords = coordinates(raster_stack) / scale_factor,
    raster_resolution = res(raster_stack) / scale_factor,
    raster_extent = c(xmin(raster_stack), xmax(raster_stack),
                      ymin(raster_stack), ymax(raster_stack)) / scale_factor,
    n_covs = nlayers(raster_stack)
  )
  
  # subsample data
  subDat <-   langSim[[isim]][sort(sample.int(nrow(langSim[[isim]]),ceiling(nrow(langSim[[isim]])/max(samplingRate,1)),replace=FALSE)),]
  subDat$dt <- do.call(c,mapply(function(x) c(0,diff(subDat$time[which(subDat$ID==x)])),1:nbAnimals,SIMPLIFY = FALSE))
  
  data <- list(model=model,
               Y=t(subDat[,c("mu.x","mu.y")])/scale_factor,
               dt=subDat$dt)
  
  # add missing observations
  data$Y[,sample.int(nrow(subDat),nrow(subDat)*propMissing,replace=FALSE)] <- NA
  
  data$isd <- as.numeric(!is.na(data$Y[1,]))
  data$obs_mod <- rep(NA,ncol(data$Y))
  data$obs_mod[data$isd==1] <- 1
  data$ID <- subDat$ID
  data$nbStates <- 1
  #data$nbSteps <- ncol(data$Y)-1
  data$nbObs <- rep(1,ncol(data$Y))
  data$M <- subDat$error_semimajor_axis / scale_factor
  data$m <- subDat$error_semiminor_axis / scale_factor
  data$c <- subDat$error_ellipse_orientation
  data$K <- matrix(NA,ncol(data$Y),2)
  data$scale_factor <- scale_factor
  data <- c(data,sca_raster_data)
  data$smoothGradient <- ifelse(npoints>0,1,0)
  data$weights <- weights
  data$zetaScale <- zetaScale
  
  # generate missing data using aniMotum
  if(propMissing>0 && any(is.na(data$Y[1,]))){
    message("   Fitting aniMotum separately to each track...")
    notNA <- which(!is.na(data$Y[1,]))
    plan(multisession,workers=min(nbAnimals,parallel::detectCores()-1))  # workers sets number of cores to use
    init.mu <- init.mu_aniMotum(subDat[notNA,],model="rw",timeSteps=data.frame(id=as.character(subDat$ID),date=as.POSIXlt(subDat$time * 1/mean(subDat$dt)))) # fit_ssm doesn't like very small \Delta_t
  } else {
    #init.mu <- t(subDat[,c("mux","muy")]) # true values
    init.mu <- t(subDat[,c("mu.x","mu.y")]) # with measurement error
  }
  
  if(model==0){
    re <- "mu"
    init.v_mu <- matrix(0,2,nrow(subDat)) 
  } else {
    re <- c("mu","v_mu")
    #init.v_mu <- t(subDat[,c("v_mux","v_muy")]) # true values
    init.v_mu <- matrix(0,2,nrow(subDat)) 
    for(i in 1:nbAnimals){
      aInd <- which(data$ID==i)
      init.v_mu[1,aInd[-1]] <- diff(init.mu[1,aInd])/(data$dt[aInd[-1]]*exp(-gamma*data$dt[aInd[-1]]))
      init.v_mu[2,aInd[-1]] <- diff(init.mu[2,aInd])/(data$dt[aInd[-1]]*exp(-gamma*data$dt[aInd[-1]]))
    }
  } 
  
  parm <- list(log_sigma=log(sigma)-log(scale_factor),beta=matrix(c(0,beta),1,length(beta)+1),mu=init.mu / scale_factor, v_mu = init.v_mu / scale_factor, log_gamma = log(gamma),
               l_delta=0,l_gamma=matrix(0,1,2),l_psi=log(psi),l_tau=c(0,0),l_rho_o=0)
  
  map <- list(l_delta=factor(NA),l_gamma=factor(c(NA,NA)),l_rho_o=factor(NA),l_tau=factor(c(NA,NA)),l_psi=factor(NA)) #,beta=factor(NA,1:(ncov+1))),
  if(model==0){
    map$v_mu <- factor(rep(NA,length(init.v_mu)))
    map$log_gamma <- factor(NA)
  }
  
  message("   Fitting model...")
  obj2 <-
    MakeADFun(
      data,
      parm,
      map = map,
      random = re,
      DLL = "fitLangevin",
      hessian = TRUE,
      silent = TRUE,
      inner.control =  list(maxit = 10000, trace=0)
    )
  
  obj2$env$tracepar <- tracepar
  langTMB[[isim]] <- do.call("nlminb",args = list(
    start = obj2$par,
    objective = obj2$fn,
    gradient = obj2$gr,
    control=list(trace=10,iter.max=10000,eval.max=10000)))
  
  message("   Calculating SEs...")
  rep <- obj2$report()
  sdrep <- sdreport(obj2)
  langTMB[[isim]]$sdreport <- summary(sdrep,"report") # get SEs 
  
  parMat[isim,1:(6+model)] <- langTMB[[isim]]$par
  parMat[isim,"sigma"] <- rep$sigma
  if(model==1) parMat[isim,"gamma"] <- rep$gamma
  
  par(mfrow=c(2,2))
  plot(UD,main=paste0("Simulation ",isim),col=viridis::viridis(100))
  points(data$Y[1,]*scale_factor,data$Y[2,]*scale_factor,col=2) # observed locations
  for(i in 1:nbAnimals){
    lines(langSim[[isim]]$mux[langSim[[isim]]$ID==i],langSim[[isim]]$muy[langSim[[isim]]$ID==i]) # true path
  }
  plot(data$Y[1,]*scale_factor,data$Y[2,]*scale_factor,col=2,asp=1) # observed locations
  for(i in 1:nbAnimals){
    lines(rep$mu[1,data$ID==i]*scale_factor,rep$mu[2,data$ID==i]*scale_factor,type="o",pch=20,col=3) # estimated path
  } 
  for(i in 1:nbAnimals){
    lines(langSim[[isim]]$mux[langSim[[isim]]$ID==i],langSim[[isim]]$muy[langSim[[isim]]$ID==i]) # true path
  }
  plot(UD,main="True UD",col=viridis::viridis(100),xlim=c(-100,100),ylim=c(-100,100))
  estUD <- langTMB[[isim]]$par[3]*spatialCovs$cov1 + langTMB[[isim]]$par[4] * spatialCovs$cov2 + langTMB[[isim]]$par[5] * spatialCovs$cov3 + langTMB[[isim]]$par[6] * spatialCovs$d2c
  plot(estUD,main="Estimated UD",col=viridis::viridis(100),xlim=c(-100,100),ylim=c(-100,100))
  
  parMat[isim,"UDcor"] <- cor(values(UD),values(estUD)) # correlation between true and estimated UD (could alternatively use raster::corLocal)
  
  print(paste0("            ",paste0(colnames(parMat),collapse="    ")))
  print(paste0("current ",paste0(round(parMat[isim,],6),collapse=" ")))
  print(paste0("overall ",paste0(round(apply(parMat[1:isim,,drop=FALSE],2,mean),6),collapse=" ")))
  
  #hist(sqrt((langSim[[isim]]$mu.x-langSim[[isim]]$mux)^2+(langSim[[isim]]$mu.y-langSim[[isim]]$muy)^2),main="measurement error",xlim=c(0,1.2),breaks=seq(0,2,length=20))
  #hist(sqrt((rep$mu[1,]-langSim[[isim]]$mux)^2+(rep$mu[2,]-langSim[[isim]]$muy)^2),main="loc estimation error",xlim=c(0,1.2),breaks=seq(0,2,length=20))
} 

save(parMat,beta,sigma,gamma,obsPerAnimal,langSim,langTMB,covlist,psi,timeStep,samplingRate,propMissing,M,scale_factor,resol,sca,npoints,file=paste0("data/",ifelse(model==1,"underdamped","overdamped"),"Sim_nbAnimals_",nbAnimals,"_obsPerAnimal_",obsPerAnimal,"_timeStep_",timeStep,"_sigma_",sigma,"_gamma_",gamma,"_beta_",paste0(beta,collapse="_"),"_samplingRate_",samplingRate,"_propMissing_",propMissing,"_M_",M,"_npoints_",npoints,"_psi_",psi,".RData"))
