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

# number of locations in land cells
leaksDat <- maskLeakage(narDat,mask, level=1.e-6, tolerance=0,coord=c("x","y"))

plotRaster(mask,extent = origextent)+geom_point(aes(x=x,y=y),data=narFilt,col="#E69F00")

narFilt <- routeTracks(narFilt,maskBuff)
leaksRoute <- maskLeakage(narFilt,mask, level=1, tolerance=res(mask)[1]*2+1,coord=c("mu.x_pr","mu.y_pr"))

onset <- 5000

narPar <- list(beta = 0, sigma = 3000, gamma = 1,
             psi = 15, psi_maj = 0.6, b_maj = 0.4, b_min = 0.4,
             tau = c(1, 1), rho_o = 0,
             mu  = cbind(narFilt$mu.x_pr, narFilt$mu.y_pr))

narFit <- fitLangevin(narFilt, model = model, spatialCovs = covs, barrier = "barrier", barrierOnset = onset,
                   obs_dist = "Gaussian", par = narPar, scaleFactor = scaleFactor,
                   control = list(trace = 1))

plot(narFit,spatialCovs=covs,data=narFilt,extent=origextent,normalize=TRUE,maskRast=maskBuff)
# number of estimated locations with no support less than 1 cell inland
narLeaks <- maskLeakage(narFit,mask,level=1,tolerance = res(mask)[1]+1.e-6)
