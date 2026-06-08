library(langevinSSM)
library(terra)
library(ggplot2)
library(dplyr)

model <- "underdamped"

sca <- 600
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")
timeStep <- 0.01
samplingRate <- 1
propMissing <- 0.4
nbAnimals <- 5
obsPerAnimal <- 5000 / (1-propMissing)
covRange <- c(0.1, 0.5)
includeBarrier <- FALSE

measurementError = list(smaj.sd = 1.5, smin.sd = 0.75)

sim_pars <- list(
  beta = c(4, -1, -0.2 , -0.1),
  sigma = 5,
  psi = 1
)
if(model=="underdamped") sim_pars$gamma <- 0.5

nSims <- 100
if(includeBarrier){
  parMat <- matrix(NA, nSims, 11)
  colnames(parMat) <- c(paste0("beta", 1:length(sim_pars$beta)), "sigma", "gamma", "lambda","muBias", "muSE", "muCov", "BA")
} else {
  parMat <- matrix(NA, nSims, 10)
  colnames(parMat) <- c(paste0("beta", 1:length(sim_pars$beta)), "sigma", "gamma","muBias", "muSE", "muCov", "BA")
}
parMat_true <- parMat

for(isim in 1:nSims){
  cat("\n============================================\n")
  cat("Simulation", isim, "\n")
  cat("============================================\n")

  sim_data <- tryCatch(stop(), error=function(e) e)
  while(inherits(sim_data, "error")){

    cov1 <- simCov(sca = sca, irange = runif(1, covRange[1], covRange[2]), sigma2 = 0.1, kappa = 0.5)
    cov2 <- simCov(sca = sca, irange = runif(1, covRange[1], covRange[2]), sigma2 = 0.1, kappa = 0.5)
    coords <- terra::crds(cov1)

    # simulate complex coastline and islands
    r <- terra::rast(nrows = 1200, ncols = 1200, ext = ext(cov1))
    terra::values(r) <- runif(terra::ncell(r))
    w_size <- sample(seq(11, 21, by = 2), 1)
    w <- matrix(1, nrow = w_size, ncol = w_size)
    noise <- r
    for (i in 1:4) {
      noise <- terra::focal(noise, w = w, fun = mean, na.policy = "omit", expand=TRUE)
    }
    n_min <- terra::global(noise, "min", na.rm = TRUE)[[1]]
    n_max <- terra::global(noise, "max", na.rm = TRUE)[[1]]
    noise <- (noise - n_min) / (n_max - n_min)
    x_coords <- terra::init(noise, "x")
    base_width <- runif(1, 0.10, 0.25)
    wiggle_amp <- runif(1, 0.15, 0.35)
    x_min <- terra::xmin(r)
    x_range <- terra::xmax(r) - x_min
    mainland <- x_coords < (x_min + x_range * base_width + noise * x_range * wiggle_amp)
    isl_thresh <- runif(1, 0.65, 0.85)
    islands <- noise > isl_thresh
    land_mask <- mainland | islands
    water_mask <- terra::ifel(land_mask, 0, 1)
    names(water_mask) <- "coast_barrier"

    dist2 <- ((coords[, "x"])^2 + coords[, "y"]^2) / sca
    d2c <- terra::setValues(cov1, dist2)
    names(d2c) <- "d2c"

    barrier <- suppressMessages(prepBarrier(water_mask))

    covs <- list(cov1 = cov1, cov2 = cov2, d2c = d2c, d2coast = barrier/100)
    if(includeBarrier) covs$coast_barrier <- barrier

    if(!includeBarrier){
      barrier <- NULL
    } else {
      barrier <- "coast_barrier"
    }

    sim_data <- tryCatch({
      out_data <- suppressMessages(simLangevin(
        model = model,
        par = sim_pars,
        spatialCovs = covs,
        nbAnimals = nbAnimals,
        obsPerAnimal = obsPerAnimal,
        barrier=barrier,
        timeStep = timeStep,
        measurementError = measurementError,
        subSample = list(samplingRate=samplingRate, propMissing=propMissing)
      ))

      if(includeBarrier){
        pts <- cbind(out_data$x, out_data$y)
        dist_vals <- terra::extract(covs$coast_barrier, pts)[, 1]
        land_idx <- which(dist_vals <= 0)

        if (length(land_idx) == 0) stop("No observed locations on land. Resimulating...")
      }

      out_data
    }, error=function(e) e)

  }



  trueUD <- getUD(covs, beta=sim_pars$beta, barrier=barrier, lambda=attr(sim_data,"lambda"), log=TRUE, plot=FALSE)

  missDat <- sim_data %>% filter(!is.na(x))
  missDat$eor <- missDat$eor * 180 / pi
  class(missDat) <- "data.frame"
  missDat <- suppressWarnings(formatData(missDat))

  lambda <- NULL
  if(includeBarrier) {
    fit0 <- tryCatch(suppressMessages(fitLangevin(
      data = missDat,
      model = model,
      spatialCovs = covs,
      barrier = barrier,
      lambda = 0,
      silent = TRUE
    )),error=function(e) e)

    if(!inherits(fit0,"error")) lambda <- suggestLambda(fit0,max(sim_data$dt))
  } else fit0 <- NULL

  if(!inherits(fit0,"error")){
    fit_miss <- tryCatch(suppressMessages(fitLangevin(
      data = missDat,
      model = model,
      spatialCovs = covs,
      barrier=barrier,
      lambda=lambda,
      silent = TRUE
    )),error=function(e) e)

    if(!inherits(fit_miss,"error")){
      fit_pred <- tryCatch(predLangevin(fit_miss,data=sim_data,spatialCovs=covs,model=model,silent=TRUE,max_iter=25),error=function(e) e)
      if(!inherits(fit_pred,"error")){
        predMu <- fitted(fit_pred,parm="mu")
        predCI <- confint(fit_pred,parm="mu")

        predmuBias <- mean(abs(as.matrix(predMu[,c("mu.x","mu.y")]-sim_data[,c("mu.x","mu.y")])))
        predmuSE <- mean(as.matrix(fit_pred$estimates$random$mu$se[,c("mu.x","mu.y")]))
        predmuCov <- mean(predCI$`mu.x_2.5%` <= sim_data$mu.x & sim_data$mu.x <= predCI$`mu.x_97.5%` &
                            predCI$`mu.y_2.5%` <= sim_data$mu.y & sim_data$mu.y <= predCI$`mu.y_97.5%` )
      }
    }
  } else fit_miss <- fit_pred <- tryCatch(stop(),error=function(e) e)

  fit_true <- tryCatch(suppressMessages(fitLangevin(
    data = sim_data,
    model = model,
    spatialCovs = covs,
    barrier=barrier,
    silent = TRUE
  )),error=function(e) e)

  if(!inherits(fit_true,"error")){
    trueMu <- fitted(fit_true,parm="mu")
    trueCI <- confint(fit_true,parm="mu")

    truemuBias <- mean(abs(as.matrix(trueMu[,c("mu.x","mu.y")]-sim_data[,c("mu.x","mu.y")])))
    truemuSE <- mean(as.matrix(fit_true$estimates$random$mu$se[,c("mu.x","mu.y")]))
    truemuCov <- mean(trueCI$`mu.x_2.5%` <= sim_data$mu.x & sim_data$mu.x <= trueCI$`mu.x_97.5%` &
                        trueCI$`mu.y_2.5%` <= sim_data$mu.y & sim_data$mu.y <= trueCI$`mu.y_97.5%` )
  }

  if(!inherits(fit_miss, "error")){
    estUD <- suppressMessages(getUD(covs, fit=fit_miss, log=TRUE, plot=FALSE))
    parMat[isim, 1:6] <- fit_miss$estimates$natural$Estimate[1:6]
    parMat[isim, "BA"] <- rasterOverlap(exp(estUD), exp(trueUD))#,local=fit_miss)
    parMat[isim, "muBias"] <- predmuBias
    parMat[isim, "muCov"] <- predmuCov
    parMat[isim, "muSE"] <- predmuSE
    if(includeBarrier) {
      parMat[isim,"lambda"] <- fit_miss$conditions$lambda
      maskRast <- water_mask
    } else {
      maskRast <- NULL
    }
    if(!inherits(fit_pred,"error")) print(plot(fit_pred,data=sim_data,spatialCovs=covs,maskRast=maskRast,log=TRUE)+labs(title=paste0("fit_pred Sim ",isim)))
  }

  if(!inherits(fit_true, "error")){
    estUD_true <- suppressMessages(getUD(covs, fit=fit_true, log=TRUE, plot=FALSE))
    parMat_true[isim, 1:6] <- fit_true$estimates$natural$Estimate[1:6]
    parMat_true[isim, "BA"] <- rasterOverlap(exp(estUD_true), exp(trueUD))#,local=fit_true)
    parMat_true[isim, "muBias"] <- truemuBias
    parMat_true[isim, "muCov"] <- truemuCov
    parMat_true[isim, "muSE"] <- truemuSE
    if(includeBarrier) {
      parMat_true[isim,"lambda"] <- fit_true$conditions$lambda
      maskRast <- water_mask
    } else {
      maskRast <- NULL
    }
    print(plot(fit_true,data=sim_data,spatialCovs=covs,maskRast=maskRast,log=TRUE)+labs(title=paste0("fit_true Sim ",isim)))
  }

  message(" predLangevin model ")
  print(paste0("            ",paste0(colnames(parMat),collapse="    ")))
  print(paste0("current ",paste0(round(parMat[isim,],6),collapse=" ")))
  print(paste0("overall ", paste0(round(apply(parMat[1:isim,,drop=FALSE], 2, mean, na.rm=TRUE), 6), collapse=" ")))

  message(" Full model ")
  print(paste0("            ",paste0(colnames(parMat_true),collapse="    ")))
  print(paste0("current ",paste0(round(parMat_true[isim,],6),collapse=" ")))
  print(paste0("overall ", paste0(round(apply(parMat_true[1:isim,,drop=FALSE], 2, mean, na.rm=TRUE), 6), collapse=" ")))

}
