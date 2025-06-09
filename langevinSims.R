options("rgdal_show_exportToProj4_warnings"="none") # suppress annoying warnings
library(tidyverse)
library(raster)
library(rasterVis)
library(viridis)
library(ggplot2)
library(doFuture)
library(doRNG)
library(RStoolbox)
library(Rcpp)
library(aniMotum)
if(!requireNamespace("RandomFieldsUtils",quietly=TRUE)){
  install.packages("http://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/RandomFieldsUtils_1.2.5.tar.gz", repos = NULL, type = "source") # most recent archived version; required by RandomFields
}
if(!requireNamespace("RandomFields",quietly=TRUE)){
  install.packages("http://cran.r-project.org/src/contrib/Archive/RandomFields/RandomFields_3.3.14.tar.gz", repos = NULL, type = "source") # most recent archived version; required by Rhabit
}
remotes::install_github("papayoun/Rhabit@31ddf44",dependencies = TRUE) # last commit before RandomFields was removed from dependencies
library(Rhabit)
if(!requireNamespace("momentuHMM",quietly=TRUE) || as.numeric(substr(packageVersion("momentuHMM"),1,1))<2){
  remotes::install_github("bmcclintock/momentuHMM@develop",dependencies = TRUE) # requires momentuHMM version >= 2.0.0
}
library(momentuHMM)

library(TMB)

## specify, compile, and load TMB model (can be "langevin" or "underdampedLangevin")
model <- "underdampedLangevin" 

compile(paste0("src/",model,".cpp"))
dyn.load(dynlib(paste0("src/",model)))

source("R/helper_functions.R")

nsims <- 100 # number of simulations
nbAnimals <- 5 # number of tracks
obsPerAnimal <- 5000 # number of simulated locations per track
lambda <- 1 # observation rate (1/lambda is expected time between successive observations)

tracepar <- FALSE # TRUE = trace parameters during optimization

## specify scale and number of raster covariates
sca <- 200 # bounding box scale
lim <- c(-1, 1, -1, 1)*sca
cropExtent <- extent(lim)
resol <- 1 # cell resolution 
ncov <- 3 # number of spatial covariates

beta <- c(-4*resol,6*resol,5*resol,-0.1*resol) # resource selection coefficients for the spatial covariates (cov_1, cov_2, ... cov_ncov, d2c)
sigma <- 1 * resol / 2 # speed parameter scaled by resol/2 
gamma <- 0.25 # friction parameter (smaller value -> more directional persistence); ignored unless model="underdampedLangevin"
psi <- 1 # error SD scaling parameter

## missing data and measurement error

propMissing <- 0.1 # proportion of missing observations

errorProp <- 0.2 # maximum error as proportion of sigma
M <- c(0,errorProp * sigma) # range for semi-major error ellipse axis
m <- c(0,errorProp/2 * sigma) # range for semi-minor error ellipse axis
r <- c(0,180) # range for ellipse orientation (degrees)


langSim <- langTMB <- list()
covlist <- list()
parMat <- matrix(NA,nrow=nsims,7+ifelse(model=="langevin",0,1))
if(model=="langevin"){
  colnames(parMat) <- c("sigma",paste0("beta",0:(ncov+1)),"UDcor")
} else colnames(parMat) <- c("sigma",paste0("beta",0:(ncov+1)),"gamma","UDcor")

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
    covlist[[isim]][[i]] <- Rhabit::simSpatialCov(lim = lim, nu = 0.6, rho = 50, sigma2 = 0.1, 
                                                  resol = resol, raster_like = TRUE)
    ## expand extent of covariates 
    #covlist[[isim]][[i]]$x <- covlist[[isim]][[4]]$x
    #covlist[[isim]][[i]]$y <- covlist[[isim]][[4]]$y
    #covlist[[isim]][[i]]$z <- rbind(matrix(0,sca/2,2*sca+1),cbind(matrix(0,sca+1,sca/2),covlist[[isim]][[i]]$z,matrix(0,sca+1,sca/2)),matrix(0,sca/2,2*sca+1))
    #covlist[[isim]][[i]]$z <- matrix(scale(c(covlist[[isim]][[i]]$z)),2*sca+1,2*sca+1) # scale covariates
  }
  
  names(covlist[[isim]]) <- c("cov1","cov2","cov3","d2c")
  spatialCovs <- lapply(lapply(covlist[[isim]],rhabitToRaster),function(x) {proj4string(x) <- CRS("+init=epsg:3416");return(x)})
  
  # orthogonalize covariates
  #pca <- RStoolbox::rasterPCA(stack(spatialCovs[1:ncov]))
  
  #pcovlist <- list()
  #for(i in 1:ncov){
  #  pcovlist[[i]] <- rasterToRhabit(raster(pca$map[[i]]))
  #  spatialCovs[[i]] <- raster(pca$map[[i]])
  #}
  #pcovlist[[ncov+1]] <- covlist[[isim]][[ncov+1]]
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
        model = ifelse(model=="langevin",0,1),
        nbAnimals = nbAnimals,
        obsPerAnimal = obsPerAnimal,
        lambda = lambda,
        gamma = gamma,
        sigma = sigma,
        beta = beta,
        raster_data = raster_data,
        initialPosition = initialPosition
      ),error=function(e) e)
    if(!inherits(langSim[[isim]],"error") && (any(M>0) | any(m>0))){
      langSim[[isim]] <- measurementError(langSim[[isim]],M=M,m=m,c=r,psi=psi)
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
  
  data <- list(Y=t(langSim[[isim]][,c("mu.x","mu.y")])/scale_factor,dt=langSim[[isim]]$dt)
  # add missing observations
  data$Y[,sample.int(nbAnimals*obsPerAnimal,nbAnimals*obsPerAnimal*propMissing,replace=FALSE)] <- NA
  data$isd <- as.numeric(!is.na(data$Y[1,]))
  data$obs_mod <- rep(NA,ncol(data$Y))
  data$obs_mod[data$isd==1] <- 1
  data$ID <- langSim[[isim]]$ID
  data$nbStates <- 1
  #data$nbSteps <- ncol(data$Y)-1
  data$nbObs <- rep(1,ncol(data$Y))
  data$M <- langSim[[isim]]$error_semimajor_axis / scale_factor
  data$m <- langSim[[isim]]$error_semiminor_axis / scale_factor
  data$c <- langSim[[isim]]$error_ellipse_orientation
  data$K <- matrix(NA,ncol(data$Y),2)
  data$scale_factor <- scale_factor
  data <- c(data,sca_raster_data)
  
  # generate missing data using aniMotum
  #if(propMissing>0 && any(is.na(data$Y[1,]))){
  #  aniDat <- data.frame(id=as.character(langSim[[isim]]$ID[which(!is.na(data$Y[1,]))]),date=as.POSIXlt(langSim[[isim]]$time[which(!is.na(data$Y[1,]))]),x=data$Y[1,which(!is.na(data$Y[1,]))]*scale_factor,y=data$Y[2,which(!is.na(data$Y[1,]))]*scale_factor,lc=3,smaj=langSim[[isim]]$error_semimajor_axis[which(!is.na(data$Y[1,]))],smin=langSim[[isim]]$error_semiminor_axis[which(!is.na(data$Y[1,]))],eor=langSim[[isim]]$error_ellipse_orientation[which(!is.na(data$Y[1,]))],x.sd=NA,y.sd=NA)
  #  aniDat <- sf::st_as_sf(
  #    x = aniDat,
  #    coords = c("x", "y"),
  #    crs = 3416,  
  #    remove = FALSE
  #  )
  #  aniFit <- aniMotum::fit_ssm(aniDat,spdf = TRUE,time.step=data.frame(id=as.character(langSim[[isim]]$ID),date=as.POSIXlt(langSim[[isim]]$time)),map = list(psi = factor(NA)))
  #  init.mu <- do.call(cbind,mapply(function(x) matrix(unlist(aniFit$ssm[[x]]$predicted$geometry),nrow=2),1:nbAnimals,SIMPLIFY = FALSE))
  #} else {
    init.mu <- t(langSim[[isim]][,c("mux","muy")]) #t(langSim[[isim]][,c("mu.x","mu.y")]) # 
  #}
        
  if(model=="langevin"){
    re <- "mu"
    init.v_mu <- matrix(0,2,nbAnimals*obsPerAnimal)
  } else {
    re <- c("mu","v_mu")
    init.v_mu <- t(langSim[[isim]][,c("v_mux","v_muy")]) # matrix(0,2,nbAnimals*obsPerAnimal) # 
  } 

  parm <- list(log_sigma=log(sigma)-log(scale_factor),beta=matrix(c(0,beta),1,length(beta)+1),mu=init.mu / scale_factor, v_mu = init.v_mu / scale_factor, log_gamma = log(gamma),
               l_delta=0,l_gamma=matrix(0,1,2),l_psi=log(psi),l_tau=c(0,0),l_rho_o=0)
  
  message("   Fitting model...")
  obj2 <-
    MakeADFun(
      data,
      parm,
      map = list(l_delta=factor(NA),l_gamma=factor(c(NA,NA)),l_rho_o=factor(NA),l_tau=factor(c(NA,NA)),l_psi=factor(NA)),#,beta=factor(NA,1:(n))),
      random = re,
      DLL = model,
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
  
  parMat[isim,1:(6+ifelse(model=="langevin",0,1))] <- langTMB[[isim]]$par
  parMat[isim,"sigma"] <- exp(parMat[isim,"sigma"] + log(scale_factor))
  if(model=="underdampedLangevin") parMat[isim,"gamma"] <- exp(parMat[isim,"gamma"])
  
  par(mfrow=c(2,2))
  plot(UD,main=paste0("Simulation ",isim),col=viridis::viridis(100))
  points(data$Y[1,]*scale_factor,data$Y[2,]*scale_factor,col=2) # observed locations
  for(i in 1:nbAnimals){
    lines(langSim[[isim]]$mux[data$ID==i],langSim[[isim]]$muy[data$ID==i]) # true path
  }
  plot(data$Y[1,]*scale_factor,data$Y[2,]*scale_factor,col=2,asp=1) # observed locations
  rep <- obj2$report()
  for(i in 1:nbAnimals){
    lines(rep$mu[1,data$ID==i]*scale_factor,rep$mu[2,data$ID==i]*scale_factor,type="o",pch=20,col=3) # estimated path
    lines(langSim[[isim]]$mux[data$ID==i],langSim[[isim]]$muy[data$ID==i]) # true path
  }
  plot(UD,main="True UD",col=viridis::viridis(100),xlim=c(-100,100),ylim=c(-100,100))
  estUD <- langTMB[[isim]]$par[3]*spatialCovs$cov1 + langTMB[[isim]]$par[4] * spatialCovs$cov2 + langTMB[[isim]]$par[5] * spatialCovs$cov3 + langTMB[[isim]]$par[6] * spatialCovs$d2c
  plot(estUD,main="Estimated UD",col=viridis::viridis(100),xlim=c(-100,100),ylim=c(-100,100))
  
  parMat[isim,"UDcor"] <- cor(values(UD),values(estUD))
  
  print(paste0("            ",paste0(colnames(parMat),collapse="    ")))
  print(paste0("current ",paste0(round(parMat[isim,],6),collapse=" ")))
  print(paste0("overall ",paste0(round(apply(parMat[1:isim,,drop=FALSE],2,mean),6),collapse=" ")))
  
  #hist(sqrt((langSim[[isim]]$mu.x-langSim[[isim]]$mux)^2+(langSim[[isim]]$mu.y-langSim[[isim]]$muy)^2),main="measurement error",xlim=c(0,1.2),breaks=seq(0,2,length=20))
  #hist(sqrt((rep$mu[1,]-langSim[[isim]]$mux)^2+(rep$mu[2,]-langSim[[isim]]$muy)^2),main="loc estimation error",xlim=c(0,1.2),breaks=seq(0,2,length=20))
}

save(parMat,beta,sigma,gamma,obsPerAnimal,langSim,langTMB,covlist,psi,lambda,propMissing,errorProp,scale_factor,resol,sca,file=paste0("data/",model,"Sim_nbAnimals_",nbAnimals,"_obsPerAnimal_",obsPerAnimal,"_lambda_",lambda,"_sd_",sigma,"_gamma_",gamma,"_beta_",paste0(beta,collapse="_"),"_propMissing_",propMissing,"_errorProp_",errorProp,"_psi_",psi,".RData"))
