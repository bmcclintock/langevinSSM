library(langevinSSM)
library(terra)
library(ggplot2)

sca <- 600
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")
timeStep <- 0.01
samplingRate <- 1
propMissing <- 0
nbAnimals <- 5
obsPerAnimal <- 5000
covRange <- c(0.1, 0.5)
lambda_max <- NULL

measurementError = list(smaj.sd = 1.5, smin.sd = 0.75)

sim_pars <- list(
  beta = c(4, -1, -0.2 , -0.1),
  sigma = 5,
  gamma = 0.5,
  psi = 1
)

nSims <- 100
parMat <- matrix(NA, nSims, 7)
colnames(parMat) <- c(paste0("beta", 1:length(sim_pars$beta)), "sigma", "gamma", "BA")
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

    covs <- list(cov1 = cov1, cov2 = cov2, d2c = d2c, coast_barrier = prepBarrier(water_mask))

    sim_data <- tryCatch({
      out_data <- simLangevin(
        model = "underdamped",
        par = sim_pars,
        spatialCovs = covs,
        nbAnimals = nbAnimals,
        obsPerAnimal = obsPerAnimal,
        timeStep = timeStep,
        measurementError = measurementError,
        subSample = list(samplingRate=samplingRate, propMissing=propMissing)
      )

      pts <- cbind(out_data$x, out_data$y)
      dist_vals <- terra::extract(covs$coast_barrier, pts)[, 1]
      land_idx <- which(dist_vals <= 0)

      if (length(land_idx) == 0) stop("No observed locations on land. Resimulating...")

      out_data
    }, error=function(e) e)

  }

  trueUD <- getUD(covs, beta=sim_pars$beta, lambda=attr(sim_data,"lambda"), log=TRUE, plot=FALSE)

  fit <- tryCatch(fitLangevin_barrier(sim_data,spatialCovs=covs, lambda_max = lambda_max, timeStep=timeStep, n_sims = 10, n_coarse=5, n_fine=5, ncores = 5, silent=TRUE),error=function(e) e)

  fit_true <- tryCatch(fitLangevin(
    data = sim_data,
    model = "underdamped",
    spatialCovs = covs,
    silent = TRUE
  ),error=function(e) e)

  if(!inherits(fit, "error")){
    estUD <- suppressMessages(getUD(covs, fit=fit, log=TRUE, plot=FALSE))
    parMat[isim, 1:6] <- fit$estimates$natural$Estimate[1:6]
    parMat[isim, "BA"] <- rasterOverlap(exp(estUD), exp(trueUD))#,local=fit)
    print(plot(fit,data=sim_data,spatialCovs=covs,log=TRUE,maskBarrier=TRUE)+labs(title=paste0("Sim ",isim)))
  }

  if(!inherits(fit_true, "error")){
    estUD_true <- suppressMessages(getUD(covs, fit=fit_true, log=TRUE, plot=FALSE))
    parMat_true[isim, 1:6] <- fit_true$estimates$natural$Estimate[1:6]
    parMat_true[isim, "BA"] <- rasterOverlap(exp(estUD_true), exp(trueUD))#,local=fit_true)
  }

  message(" Iterative lambda ")
  print(paste0("            ",paste0(colnames(parMat),collapse="    ")))
  print(paste0("current ",paste0(round(parMat[isim,],6),collapse=" ")))
  print(paste0("overall ", paste0(round(apply(parMat[1:isim,,drop=FALSE], 2, mean, na.rm=TRUE), 6), collapse=" ")))

  message(" True lambda ")
  print(paste0("            ",paste0(colnames(parMat_true),collapse="    ")))
  print(paste0("current ",paste0(round(parMat_true[isim,],6),collapse=" ")))
  print(paste0("overall ", paste0(round(apply(parMat_true[1:isim,,drop=FALSE], 2, mean, na.rm=TRUE), 6), collapse=" ")))

}
