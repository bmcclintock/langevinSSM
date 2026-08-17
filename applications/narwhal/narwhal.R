#remotes::install_github('bmcclintock/langevinSSM')
library(langevinSSM)
library(terra)
library(readr)
library(dplyr)
library(ggplot2)

set.seed(1,kind="Mersenne-Twister",normal.kind = "Inversion")

source("applications/narwhal/getCovs.R")
origextent <- ext(c(469606.176308853, 628137.927473616, 7968123.93718224, 8212082.12420781))

scaleFactor <- 1000

model <- "underdamped"

# expand raster boundaries by 200 kilometers in every direction
covs_buffered <- getCovs(buffer_meters = 200000, target_res = 821.408)

# flip the mask so water = 1 (allowed) and land = 0 (restricted)
mask <- terra::ifel(covs_buffered$mask == 0, 1, 0)
plotRaster(mask,extent=origextent)
maskBuff <- maskBuffer(mask,bufferCells = 0)
plotRaster(maskBuff,extent=origextent)

barrier <- prepBarrier(maskBuff)

d2c <- barrier

bathy <- covs_buffered$bathy

plot(bathy, main = "bathymetry")
plotRaster(d2c,extent=origextent)
plot(barrier, main = "signed distance field")

narRaw <- read_csv("https://raw.githubusercontent.com/Fanny-Dupont/Langevin_SSM/284de22e4746cc393a9019f0f13c8ee6fef45dcc/Code_Case_Study/Narwhal_Case_Study.csv")
narRaw <- narRaw %>% mutate(X=X*1000,Y=Y*1000,smaj=smaj*1000,smin=smin*1000)

# Set Fastloc to LS with 50m error and Argos to KF
narRaw$x.err <- narRaw$y.err <- NA
floc_idx <- which(narRaw$loc_class=="G")
narRaw$smaj[floc_idx] <- narRaw$smin[floc_idx] <- narRaw$eor[floc_idx] <- NA
narRaw$x.err[floc_idx] <- narRaw$y.err[floc_idx] <- 50

narDat <- narFilt <- formatData(narRaw, id="ID", date="datetime_UTC", coord=c("X","Y"), lc = "loc_class")

covs <- list(bathy=scale(bathy * maskBuff), barrier = barrier)

leaksDat <- maskLeakage(narDat,mask, level=1, tolerance=0,coord=c("x","y"))

plotRaster(mask,extent = origextent)+geom_point(aes(x=x,y=y),data=narFilt,col="#E69F00")

# fit model with no barrier penalty
narFit0 <- fitLangevin(narFilt, model = model, spatialCovs = covs, barrier = "barrier",
                       par = list(psi=1, tau=c(1,1)),
                       scaleFactor=scaleFactor,
                       lambda = 0,
                       polishOptim = TRUE,
                       silent = TRUE, control=list(trace=1))

lambda <- suggestLambda(narFit0,
                        max_dt=median(narFilt$dt))

narFilt <- routeTracks(narFilt,maskBuff)
leaksRoute <- maskLeakage(narFilt,mask, level=1, tolerance=0,coord=c("mu.x_pr","mu.y_pr"))

# random initial values that yield maximum log-likelihood
narPar <- list(beta = -0.588223050464112,
               sigma = 2599.48381240145,
               gamma = 1.99934981183229,
               psi = 21.0399917348219,
               tau = c(1.42885457494551, 1.25090566160653))

narPar$mu <- cbind(narFilt$mu.x_pr,narFilt$mu.y_pr)

narFit <- fitLangevin(narFilt, model = model, spatialCovs = covs, barrier = "barrier",
                      par = narPar,
                      scaleFactor=scaleFactor,
                      lambda = lambda,
                      control=list(trace=1),
                      polishOptim = TRUE,
                      silent = TRUE)
narFit
plot(narFit,spatialCovs=covs,data=narFilt,extent=origextent,normalize=TRUE,maskRast=maskBuff)
narLeaks <- maskLeakage(narFit,mask,level=1,tolerance = res(mask)[1]*2+1)

narRes <- residuals(narFit,data=narFilt,spatialCovs=covs)
narRes
plot(narRes)

narSim <- simLangevin(narFit,data=narFilt,spatialCovs=covs,timeStep="1 min")
leaksSim <- maskLeakage(narSim,mask,level=1,tolerance=res(mask)[1]*2+1)
plotRaster(mask,extent = origextent)+geom_point(aes(x=mu.x,y=mu.y),data=narSim,col="#E69F00")

narBest <- narFit
parBest <- narPar

for(i in 1:100){
  message("Simulation ",i)
  #par <- list(beta=rnorm(1,0,2),
  #            sigma = runif(1,1500,3000),
  #            gamma = runif(1,0.5,4),
  #            psi = runif(1,1,5),
  #            tau = runif(2,0.5,10))

  par <- getPar(narBest)
  par$beta <- par$beta + rnorm(1,0,1)
  par$sigma <- par$sigma + runif(1,-100,100)
  par$gamma <- par$gamma + runif(1,-1,1)
  par$psi <- max(par$psi + runif(1,-1,1),1)
  par$tau <- pmax(par$tau + runif(2,-0.5,2),0.5)
  par$mu <- cbind(narFilt$mu.x_pr,narFilt$mu.y_pr)
  par$vel <- NULL

  fit <- tryCatch(suppressWarnings(suppressMessages(fitLangevin(narFilt, model = model, spatialCovs = covs, barrier = "barrier",
                                                                par = par,
                                                                scaleFactor=scaleFactor,
                                                                lambda = lambda,
                                                                #control=list(trace=1),
                                                                polishOptim=TRUE,
                                                                silent = TRUE))),error=function(e) e)

  if(!inherits(fit,"error")){
    print(fit$objective)
    if(AIC(fit)<=AIC(narBest)) {
      narBest <- fit
      parBest <- par
      print(narBest)
      print(plot(narBest,spatialCovs=covs,data=narFilt,extent=origextent,normalize=TRUE,maskRast=maskBuff))
      leaksBest <- maskLeakage(narBest,mask,level=1,tolerance = res(mask)[1]*2+1)
    }
  }
}

