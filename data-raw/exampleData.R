
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
measurementError <- list(smaj.sd=1.5,smin.sd=0.75,eor.lim=c(0,180))
exampleCovs <- lapply(rep(sca,ncov-1),simCov)
exampleCovs <- lapply(1:(ncov-1),function(x) {
  names(exampleCovs[[x]]) <- paste0("cov",x);
  return(exampleCovs[[x]])
})
coords <- terra::crds(exampleCovs[[1]])
dist2 <- (coords[, "x"]^2 + coords[, "y"]^2) / sca

exampleCovs[[4]] <- exampleCovs[[1]]
terra::values(exampleCovs[[4]]) <- dist2
names(exampleCovs[[4]]) <- "d2c"

names(exampleCovs) <- c(paste0("cov",1:(ncov-1)),"d2c")

# shift covariates so coordinates don't look like they could be lat/long
exampleCovs <- lapply(exampleCovs,terra::shift, dx=1000, dy=1000)

UD <- getUD(exampleCovs, beta = examplePar$beta)

set.seed(kind="Mersenne-Twister",normal.kind="Inversion",seed=1)
exampleDat <- simLangevin(par=examplePar,spatialCovs=exampleCovs,nbAnimals=nbAnimals,obsPerAnimal=obsPerAnimal,measurementError = measurementError)

plotUD(UD, extent=c(900,1100,900,1100))+geom_point(aes(x=x,y=y),data=exampleDat,col=2)+geom_point(aes(x=mu.x,y=mu.y),data=exampleDat)

fit <- fitLangevin(exampleDat,spatialCovs = exampleCovs,silent=TRUE,control=list(trace=1))
fit
res <- residuals(fit,exampleDat, exampleCovs, run_tests = TRUE, ncores=nbAnimals)
res

sim <- simLangevin(model=fit,data=exampleDat,spatialCovs=exampleCovs,jointPrecision = TRUE)

estUD <- getUD(fit,spatialCovs=exampleCovs,nSims=1000)
rasterOverlap(exp(estUD),exp(UD))

plot(fit,spatialCovs=exampleCovs,data=exampleDat)

p <- plot(res)
p$qq_x + p$qq_y + p$acf_x + p$acf_y + plot_layout(ncol=2)

d2c <- exampleCovs$d2c < 2.5
reg_prob <- regionProb(fit,
                       spatialCovs = exampleCovs,
                       mask = d2c, # region of interest
                       nSims = 1000)
reg_prob

# convert mask to SDF and add to the spatial covariates list
coast_barrier <- exampleCovs[[1]]
terra::values(coast_barrier) <- ifelse(terra::crds(coast_barrier)[, "x"]
                                       >= mean(terra::crds(coast_barrier)[, "x"]), 1, 0)
names(coast_barrier) <- "coast_barrier"
exampleCovs_barrier <- exampleCovs
exampleCovs_barrier$coast_barrier <- prepBarrier(coast_barrier)

# add a beta coefficient for the barrier to the parameter list
par_barrier <- examplePar
par_barrier$beta <- c(par_barrier$beta, -0.2)

# simulate the data
set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
simDat_barrier <- simLangevin(par = par_barrier,
                              nbAnimals = 3,
                              spatialCovs = exampleCovs_barrier,
                              measurementError = list(smaj.sd = 1.5,
                                                      smin.sd = 0.75,
                                                      eor.lim = c(0,180)))

# Because simDat_barrier is a simLangevin object, fitLangevin will automatically
# detect and use the exact barrier penalty (lambda) that generated the data
fit_barrier <- fitLangevin(data = simDat_barrier,
                           spatialCovs = exampleCovs_barrier,
                           silent = TRUE)
fit_barrier

plot(fit_barrier,data=simDat_barrier,spatialCovs = exampleCovs_barrier, maskRast = coast_barrier)


res_barrier <- residuals(fit_barrier,simDat_barrier, exampleCovs_barrier, run_tests = TRUE, ncores=nbAnimals)



fit_barrier_ks <- tuneBarrier(data = simDat_barrier,
                              spatialCovs = exampleCovs_barrier,
                              silent = TRUE)

plot(fit_barrier_ks,data=simDat_barrier,spatialCovs = exampleCovs_barrier, maskRast = coast_barrier)


start_time <- as.POSIXct(paste(Sys.Date(), "00:00:00"), tz = "UTC")
exampleDat$date <- start_time + (exampleDat$date * 3600)
exampleDat$lc <- NA_character_
exampleDat$mu.x <- exampleDat$mu.y <- exampleDat$vel.x <- exampleDat$vel.y <- NULL
exampleDat <- exampleDat[,c("id","date","dt","x","y","lc","smaj","smin","eor","x.err","y.err")]
exampleDat$lc <- as.factor(exampleDat$lc)
exampleDat$id <- as.factor(exampleDat$id)
attr(exampleDat,"time.unit") <- "hours"
class(exampleDat) <- c("dataLangevin","data.frame")

unformatDat <- exampleDat
unformatDat$dt <- unformatDat$lc <- NULL
unformatDat$mu.x <- unformatDat$mu.y <- unformatDat$vel.x <- unformatDat$vel.y <- NULL
unformatDat$eor <- unformatDat$eor * 180 / pi
class(unformatDat) <- "data.frame"
attr(unformatDat,"time.unit") <- NULL

reformatDat <- formatData(unformatDat, time.unit="hours")

if(!all.equal(exampleDat,reformatDat)) stop("example data sets don't match")


lapply(1:ncov,function(x) terra::writeRaster(exampleCovs[[x]], paste0("inst/extdata/exampleCov",x,".tif"),overwrite=TRUE))
usethis::use_data(exampleDat,unformatDat,compress="xz",overwrite=TRUE)


