options("rgdal_show_exportToProj4_warnings"="none") # suppress annoying warnings
library(tidyverse)
library(terra)
library(viridis)
library(ggplot2)
library(Rcpp)
library(aniMotum)
library(fields)
library(TMB)
library(future)
library(furrr)

## specify, compile, and load TMB model
compile("src/fitLangevin.cpp")
dyn.load(dynlib(paste0("src/fitLangevin")))
source("R/helper_functions.R")
source("R/simLangevin.R")
source("R/fitLangevin.R")

model <- c("underdamped","overdamped")[1]

nsims <- 100 # number of simulations
nbAnimals <- 5 # number of tracks
obsPerAnimal <- 5000 # number of simulated locations per track
timeStep <- 0.01 # time scale of simulation (should be small to help prevent discretization error)

beta <- c(-4, 6, 5, -0.1) # resource selection coefficients for the spatial covariates (cov_1, cov_2, ... cov_ncov, d2c)
ncov <- length(beta) - 1 # number of spatial covariates
sigma <- 5 # speed parameter 
gamma <- 0.5 # friction parameter (smaller value -> more directional persistence); ignored unless model=="underdamped"
psi <- 1 # error SD scaling parameter

## sampling rate, missing data, and measurement error
samplingRate <- 1 # for subsampling observations from true continuous-time model (e.g. if samplingRate = 2 then data are roughly thinned by 2); must be >= 1; note bias increases with samplingRate, but relationships largely preserved
propMissing <- 0 # proportion of missing observations; passed as NA observations to TMB (so corresponding true locations treated as random effects to be estimated)

M <- 1.5 # SD for semi-major error ellipse axis ~ abs(Normal(0,M))
m <- M/2 # SD for semi-minor error ellipse axis ~ abs(Normal(0,M/2))
r <- c(0,180) # range for error ellipse orientation (degrees)
measurementError <- list(M=M,m=m,c=r) # setting measurementError <- NULL adds no measurement error to observations

## specify scale and spatial autocorrelation for covariates
sca <- 200 # bounding box scale
covRange <- c(0.1,0.5) # lower and upper bounds for covariate spatial range parameter (lower has less spatial autocorrelation)

## smooth gradient specifications (based on adjacent cells)
## slows model fitting, but can potentially reduce bias if \Delta_t is large (i.e. samplingRate > 1) relative to timeStep (Blackwell & Matthiopoulos 2024, https://doi.org/10.1002/ecy.4457) or if measurement error is large relative to raster cell resolution
npoints <- 0 # number of smoothing points around current cell (0 = none, 4 = diagonal, 8 = queen neighborhood);
curweight <- 1/2 # smoothing weight of current cell location; ignored if npoints=0
weights <- c(curweight,rep((1-curweight)/npoints,npoints)) # smoothing weights (must be of length 5 or 9 and sum to 1); ignored if npoints=0
zetaScale <- 2 # scale factor for smooth gradient neighborhood (>1 increases, <1 decreases neighborhood); underdamped Langevin neighborhood = zetaScale * sqrt(2*pi) * sigma; overdamped langevin neighborhood = zetaScale * sigma / sqrt(2); ignored if npoints=0

langSim <- langTMB <- list()
spatialCovs <- list()
parMat <- matrix(NA,nrow=nsims,6+ifelse(model=="underdamped",1,0))
if(model=="overdamped"){
  colnames(parMat) <- c(paste0("beta",1:(ncov+1)),"sigma","UDcor")
} else colnames(parMat) <- c(paste0("beta",1:(ncov+1)),"sigma","gamma","UDcor")

set.seed(1,kind="Mersenne-Twister",normal.kind = "Inversion")

for(isim in 1:nsims){
  
  cat("Simulation",isim,"\n")
  #######################
  ## Define covariates ##
  #######################
  # Generate ncov spatial covariates
  message("   Generating covariates...")
  spatialCovs[[isim]] <- list()
  
  for(i in 1:ncov) {
    irange <- runif(1,covRange[1],covRange[2])
    spatialCovs[[isim]][[i]] <- simCov(sca = sca, irange=irange, sigma2 = 0.1, kappa = 0.5, M = 2048, N = 2048)
    # terra::crs(spatialCovs[[isim]][[i]]) <- "epsg:3416"
  }
  
  coords <- terra::crds(spatialCovs[[isim]][[1]])
  dist2 <- (coords[, "x"]^2 + coords[, "y"]^2) / sca
  
  spatialCovs[[isim]][[4]] <- spatialCovs[[isim]][[1]]
  terra::values(spatialCovs[[isim]][[4]]) <- dist2
  
  names(spatialCovs[[isim]]) <- c("cov1","cov2","cov3","d2c")
  
  # Compute utilization distribution 
  UD <- getUD(spatialCovs[[isim]], beta=beta,log=TRUE)
  
  # Plot covariates
  ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                   legend.title = element_text(size=12), legend.text = element_text(size=12))
  c1plot <- plotRaster(spatialCovs[[isim]][[1]], scale.name = expression(c[1])) + ggtheme
  c2plot <- plotRaster(spatialCovs[[isim]][[2]], scale.name = expression(c[2])) + ggtheme
  c3plot <- plotRaster(spatialCovs[[isim]][[3]], scale.name = expression(c[3])) + ggtheme
  UDplot <- plotRaster(UD, scale.name = expression(pi)) + ggtheme
  #UDplot
  
  if(model=="underdamped"){
    par <- list(beta=beta,sigma=sigma,gamma=gamma,psi=psi)
  } else {
    par <- list(beta=beta,sigma=sigma,psi=psi)
  }
  
  # simulate "high resolution" tracks; this can take a while...
  langSim[[isim]] <- tryCatch(stop(),error=function(e) e)
  while(inherits(langSim[[isim]],"error")){
    langSim[[isim]] <- tryCatch(simLangevin(
      model = model,
      nbAnimals = nbAnimals,
      obsPerAnimal = obsPerAnimal,
      timeStep = timeStep,
      par=par,
      spatialCovs = spatialCovs[[isim]],
      measurementError = measurementError
    ),error=function(e) e)
    if(inherits(langSim[[isim]],"error")){
      message("    Retrying Simulation ",isim,": ",langSim[[isim]]$message)
    }
  }
  
  #UDplot + geom_point(mapping=aes(x=x,y=y),data=langSim[[isim]])
  
  # subsample data
  probs <- rep(1,nrow(langSim[[isim]]))
  probs[cumsum(c(1,table(langSim[[isim]]$id)[1:(nbAnimals-1)]))] <- 1.e+10 # ensure first observation is sampled
  subDat <-   langSim[[isim]][sort(sample.int(nrow(langSim[[isim]]),ceiling(nrow(langSim[[isim]])/max(samplingRate,1)),prob=probs,replace=FALSE)),]
  subDat$dt <- do.call(c,mapply(function(x) c(0,diff(subDat$date[which(subDat$id==x)])),1:nbAnimals,SIMPLIFY = FALSE))
  
  # add missing observations
  probs <- rep(1,nrow(subDat))
  probs[cumsum(c(1,table(subDat$id)[1:(nbAnimals-1)]))] <- 0 # don't let first observation be missing
  subDat[sample.int(nrow(subDat),nrow(subDat)*propMissing,prob=probs,replace=FALSE),c("x","y","smaj","smin","eor")] <- NA
  
  # generate missing data using aniMotum
  if(propMissing>0 && any(is.na(subDat$x))){
    message("   Fitting aniMotum separately to each track...")
    notNA <- which(!is.na(subDat$x))
    future::plan(future::multisession,workers=min(nbAnimals,parallel::detectCores()-1))  # workers sets number of cores to use
    init.mu <- init.mu_aniMotum(subDat[notNA,],model="rw",timeSteps=data.frame(id=as.character(subDat$id),date=as.POSIXlt(subDat$date * 1000))) # fit_ssm doesn't like very small \Delta_t
  } else {
    #init.mu <- subDat[,c("mu.x","mu.y")] # true values
    init.mu <- subDat[,c("x","y")] # with measurement error
  }
  
  if(model=="underdamped"){
    #init.v_mu <- subDat[,c("v.x","v.y")] # true values
    init.v_mu <- matrix(0,nrow(subDat),2) 
    for(i in 1:nbAnimals){
      aInd <- which(subDat$id==i)
      init.v_mu[aInd[-1],1] <- diff(init.mu[aInd,1])/(subDat$dt[aInd[-1]]*exp(-gamma*subDat$dt[aInd[-1]]))
      init.v_mu[aInd[-1],2] <- diff(init.mu[aInd,2])/(subDat$dt[aInd[-1]]*exp(-gamma*subDat$dt[aInd[-1]]))
    }
  } 
  
  par$mu=init.mu
  if(model=="underdamped") par$v_mu = init.v_mu

  langTMB[[isim]] <- fitLangevin(subDat,model=model,par=par,spatialCovs=spatialCovs[[isim]],
                                 map=list(psi=factor(NA)),
                                 smoothGradient = ifelse(npoints>0,TRUE,FALSE), npoints = npoints, curweight = curweight, zetaScale = zetaScale, 
                                 silent=TRUE, control=list(trace=0))
  
  parMat[isim,1:(5+ifelse(model=="underdamped",1,0))] <- langTMB[[isim]]$par
  parMat[isim,"sigma"] <- langTMB[[isim]]$estimates$natural["sigma",1]
  if(model=="underdamped") parMat[isim,"gamma"] <- langTMB[[isim]]$estimates$natural["gamma",1]
  
  zoom_ext <- terra::ext(-100, 100, -100, 100)
  
  par(mfrow=c(2,2))
  plot(UD[[1]], main = paste0("Simulation ", isim), col = viridis::viridis(100))
  points(subDat$x, subDat$y, col = 2) # observed locations
  for(i in 1:nbAnimals) {
    lines(langSim[[isim]]$mu.x[langSim[[isim]]$id==i], langSim[[isim]]$mu.y[langSim[[isim]]$id==i]) # true path
  }
  plot(subDat$x,subDat$y,col=2,asp=1) # observed locations
  for(i in 1:nbAnimals){
    lines(langTMB[[isim]]$report$mu[1,subDat$id==i],langTMB[[isim]]$report$mu[2,subDat$id==i],type="o",pch=20,col=3) # estimated path
  } 
  for(i in 1:nbAnimals){
    lines(langSim[[isim]]$mu.x[langSim[[isim]]$id==i],langSim[[isim]]$mu.y[langSim[[isim]]$id==i]) # true path
  }
  plot(terra::crop(UD,zoom_ext),main="True UD",col=viridis::viridis(100))
  estUD <- getUD(spatialCovs[[isim]],langTMB[[isim]]$par[1:4],log=TRUE)
  plot(terra::crop(estUD,zoom_ext),main="Estimated UD",col=viridis::viridis(100))
  
  parMat[isim,"UDcor"] <- cor(values(UD),values(estUD)) # correlation between true and estimated UD (could alternatively use raster::corLocal)
  
  print(paste0("            ",paste0(colnames(parMat),collapse="    ")))
  print(paste0("current ",paste0(round(parMat[isim,],6),collapse=" ")))
  print(paste0("overall ",paste0(round(apply(parMat[1:isim,,drop=FALSE],2,mean),6),collapse=" ")))
  
} 

save(parMat,beta,sigma,gamma,obsPerAnimal,langSim,langTMB,covlist,psi,timeStep,samplingRate,propMissing,M,sca,npoints,covRange,file=paste0("data/",model,"Sim_nbAnimals_",nbAnimals,"_obsPerAnimal_",obsPerAnimal,"_timeStep_",timeStep,"_covRange_",covRange,"_sigma_",sigma,"_gamma_",gamma,"_beta_",paste0(beta,collapse="_"),"_samplingRate_",samplingRate,"_propMissing_",propMissing,"_M_",M,"_npoints_",npoints,"_psi_",psi,".RData"))
