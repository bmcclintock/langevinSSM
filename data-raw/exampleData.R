
#' Example data simulation
#'
#' Generate data used in other functions' examples and unit tests.
#'
library(ggplot2)
library(usethis)
library(langevinSSM)
library(patchwork)

set.seed(kind="Mersenne-Twister",normal.kind="Inversion",seed=1)

# simulate data
nbAnimals <- 3
obsPerAnimal <- 500
sca <- 200
examplePar=list(beta=c(-4, 6, 5, -0.1),sigma=5,gamma=0.5)
ncov <- length(examplePar$beta)
measurementError <- list(smaj.sd=1.5,smin.sd=0.75,eor=c(0,180))
exampleCovs <- lapply(rep(sca,ncov-1),simCov)
coords <- terra::crds(exampleCovs[[1]])
dist2 <- (coords[, "x"]^2 + coords[, "y"]^2) / sca

exampleCovs[[4]] <- exampleCovs[[1]]
terra::values(exampleCovs[[4]]) <- dist2

names(exampleCovs) <- c(paste0("cov",1:(ncov-1)),"d2c")

UD <- getUD(exampleCovs,examplePar$beta)

set.seed(kind="Mersenne-Twister",normal.kind="Inversion",seed=1)
exampleDat <- simLangevin(par=examplePar,spatialCovs=exampleCovs,nbAnimals=nbAnimals,obsPerAnimal=obsPerAnimal,measurementError = measurementError)

plotRaster(UD)+geom_point(aes(x=x,y=y),data=exampleDat,col=2)+geom_point(aes(x=mu.x,y=mu.y),data=exampleDat)

fit <- fitLangevin(exampleDat,spatialCovs = exampleCovs,silent=TRUE,control=list(trace=1),calcOSA=TRUE)
fit

plot(fit,spatialCovs=exampleCovs,data=exampleDat)

p <- plotResiduals(fit)
p$qq_x + p$qq_y + p$acf_x + p$acf_y + plot_layout(ncol=2)

exampleDat$date <- as.POSIXlt(exampleDat$date*100*60, tz = "UTC")
exampleDat <- exampleDat[,c("id","date","dt","x","y","smaj","smin","eor","x.sd","y.sd","mu.x","mu.y","vel.x","vel.y")]
#exampleDat$dt <- NULL
#exampleDat$eor <- exampleDat$eor * 180 / pi
attr(exampleDat,"time.unit") <- "mins"

lapply(1:ncov,function(x) terra::writeRaster(exampleCovs[[x]], paste0("inst/extdata/exampleCov",x,".tif"),overwrite=TRUE))
usethis::use_data(exampleDat,compress="xz",overwrite=TRUE)


