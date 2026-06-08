#remotes::install_github('bmcclintock/langevinSSM')
library(langevinSSM)
library(terra)
library(readr)
library(dplyr)
library(ggplot2)

set.seed(1,kind="Mersenne-Twister",normal.kind = "Inversion")

source("getCovs.R")
origextent <- ext(c(469606.176308853, 628137.927473616, 7968123.93718224, 8212082.12420781))

scaleFactor <- 1000

model <- "underdamped"

# expand raster boundaries by 200 kilometers in every direction
covs_buffered <- getCovs(buffer_meters = 200000, target_res = 821.408)

# flip the mask so water = 1 (allowed) and land = 0 (restricted)
mask <- terra::ifel(covs_buffered$mask == 0, 1, 0)
plotRaster(mask,extent=origextent)
maskBuff <- maskBuffer(mask,bufferCells = 2)
plotRaster(maskBuff,extent=origextent)

barrier <- prepBarrier(maskBuff)

d2c <- barrier / scaleFactor

bathy <- covs_buffered$bathy / scaleFactor

plot(bathy, main = "bathymetry")
plot(d2c, main = "distance to land")
plot(barrier, main = "signed distance field")

narRaw = read_csv("Narwhal_Case_Study.csv")
narRaw <- narRaw %>% mutate(X=X*1000,Y=Y*1000,smaj=smaj*1000,smin=smin*1000)

# Set Fastloc to LS with 50m error and Argos to KF
narRaw$x.err <- narRaw$y.err <- NA
floc_idx <- which(narRaw$loc_class=="G")
narRaw$smaj[floc_idx] <- narRaw$smin[floc_idx] <- narRaw$eor[floc_idx] <- NA
narRaw$x.err[floc_idx] <- narRaw$y.err[floc_idx] <- 50

narDat <- narFilt <- formatData(narRaw, id="ID", date="datetime_UTC", coord=c("X","Y"), lc = "loc_class")

covs <- list(bathy=bathy * mask, barrier = barrier)

leaksDat <- maskLeakage(narDat,mask, level=1, tolerance=0,coord=c("x","y"))

# remove impossible Argos locations
makeNA <- leaksDat$leaked_data$obs_index[which(leaksDat$leaked_data$type == "Argos")]
narFilt[makeNA, c("x","y","smaj","smin","eor")] <- NA

plotRaster(mask,extent = origextent)+geom_point(aes(x=x,y=y),data=narFilt,col="#E69F00")

# fit model with no barrier penalty
narFit0 <- fitLangevin(narFilt, model = model, spatialCovs = covs, barrier = "barrier",
                       par = list(psi=1, tau=c(1,1)),
                       scaleFactor=scaleFactor,
                       lambda = 0,
                       silent = TRUE, control=list(trace=1))

lambda <- suggestLambda(narFit0,
                        max_dt=median(narFilt$dt))

# random initial values that yield maximum log-likelihood
narPar <- list(beta=-0.1279031,
               sigma=1542.965,
               gamma=4.631041,
               psi=13.89862,
               tau=c(1.163340, 2.325975))

narFit <- fitLangevin(narFilt, model = model, spatialCovs = covs, barrier = "barrier",
                               par = narPar,
                               scaleFactor=scaleFactor,
                               lambda = lambda,
                               control=list(trace=1),
                               silent = TRUE)
narFit
plot(narFit,spatialCovs=covs,data=narFilt,extent=origextent,maskRast=maskBuff)
narLeaks <- maskLeakage(narFit,mask,level=1,tolerance = res(mask)[1]*2+1)

narPred <- predLangevin(narFit, data=narFilt, spatialCovs=covs, max_iter=1000, silent=TRUE)
plot(narPred,spatialCovs=covs,data=narFilt,extent=origextent,maskRast=maskBuff)
leaksPred <- maskLeakage(narPred,mask,level=1,tolerance=res(mask)[1]*2+1) # leakage tolerance of 2 cells

parPred <- getPar(narPred)
parPred$rho_o <- NULL
narFit <- fitLangevin(narFilt, model = model, spatialCovs = covs, barrier = "barrier",
                      par = parPred,
                      scaleFactor=scaleFactor,
                      lambda = lambda,
                      control=list(trace=1),
                      silent = TRUE)
narFit
plot(narFit,spatialCovs=covs,data=narFilt,extent=origextent,maskRast=maskBuff)
narLeaks <- maskLeakage(narFit,mask,level=1,tolerance = res(mask)[1]*2+1)

narRes <- residuals(narFit,data=narFilt,spatialCovs=covs)
narRes
plot(narRes)

narSim <- simLangevin(narFit,data=narFilt,spatialCovs=covs,timeStep="1 min")
leaksSim <- maskLeakage(narSim,mask,level=1,tolerance=res(mask)[1]*2+1)
plotRaster(mask,extent = origextent)+geom_point(aes(x=mu.x,y=mu.y),data=narSim,col="#E69F00")



