library(langevinSSM)
library(terra)
library(ggplot2)

sca <- 200
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")
timeStep <- 0.01
samplingRate <- 1
propMissing <- 0
nbAnimals <- 5
obsPerAnimal <- 5000
covRange <- c(0.1, 0.5)

measurementError = list(smaj.sd = 7.5, smin.sd = 2.5, eor.lim=c(45,90))

sim_pars <- list(
  beta = c(6, -1, -0.2 , -0.1),
  sigma = 10,
  gamma = 2,
  psi = 1
)

nSims <- 100
parMat <- matrix(NA, nSims, 7)
colnames(parMat) <- c(paste0("beta", 1:length(sim_pars$beta)), "sigma", "gamma", "BA")
parMat_filter <- parMat_suggest <- parMat
parMat_psi <- matrix(NA, nSims, 8)
colnames(parMat_psi) <- c(paste0("beta", 1:length(sim_pars$beta)), "sigma", "gamma","psi", "BA")

for(isim in 1:nSims){
  cat("\n============================================\n")
  cat("Simulation", isim, "\n")
  cat("============================================\n")

  sim_data <- barrier_spec <- tryCatch(stop(), error=function(e) e)
  while(inherits(sim_data, "error") || inherits(barrier_spec,"error")){

    cov1 <- simCov(sca = sca, irange = runif(1, covRange[1], covRange[2]), sigma2 = 0.1, kappa = 0.5)
    cov2 <- simCov(sca = sca, irange = runif(1, covRange[1], covRange[2]), sigma2 = 0.1, kappa = 0.5)
    coords <- terra::crds(cov1)

    r <- terra::rast(nrows = 400, ncols = 400, ext = ext(cov1))
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

    covs <- list(cov1 = cov1, cov2 = cov2, d2c = d2c, coast_barrier = water_mask)

    sim_data <- tryCatch({
      out_data <- suppressMessages(simLangevin(
        model = "underdamped",
        par = sim_pars,
        spatialCovs = covs,
        nbAnimals = nbAnimals,
        obsPerAnimal = obsPerAnimal,
        timeStep = timeStep,
        barrier = "coast_barrier",
        measurementError = measurementError,
        subSample = list(samplingRate=samplingRate, propMissing=propMissing)
      ))

      sdf_rast <- langevinSSM:::.get_barrier_sdf("coast_barrier", covs)
      pts <- cbind(out_data$x, out_data$y)
      dist_vals <- terra::extract(sdf_rast, pts)[, 1]
      land_idx <- which(dist_vals <= 0)

      if (length(land_idx) == 0) stop("No observed locations on land. Resimulating...")

      attr(out_data, "land_idx") <- land_idx
      attr(out_data, "dist_vals") <- dist_vals

      out_data
    }, error=function(e) e)

    if(!inherits(sim_data,"error")){
      barrier_spec <- tryCatch(checkBarrier(
        data = sim_data,
        spatialCovs=covs,
        barrier="coast_barrier",
        silent=TRUE
      ),error=function(e) e)
    }
  }

  trueUD <- getUD(covs, beta=sim_pars$beta, barrier="coast_barrier", lambda=attr(sim_data,"lambda"), log=TRUE, plot=FALSE)

  land_idx <- attr(sim_data, "land_idx")
  keep_idx <- setdiff(1:nrow(sim_data), land_idx)

  data_filter <- sim_data[keep_idx, ]
  class(data_filter) <- class(sim_data)
  attr(data_filter, "time.unit") <- attr(sim_data, "time.unit")
  attr(data_filter, "lambda") <- attr(sim_data, "lambda")
  attr(data_filter, "proj") <- attr(sim_data, "proj")

  data_suggest <- barrier_spec$filtered_data

  message("Fitting checkBarrier model...")
  fit_suggest <- tryCatch(suppressMessages(fitLangevin(
    data = data_suggest,
    model = "underdamped",
    spatialCovs = covs,
    barrier = "coast_barrier",
    lambda = barrier_spec$recommended,
    silent = TRUE
  )), error=function(e) e)

  message("Fitting filter model...")
  fit_filter <- tryCatch(suppressMessages(fitLangevin(
    data = data_filter, model = "underdamped", spatialCovs = covs,
    barrier = "coast_barrier", lambda = barrier_spec$recommended, silent = TRUE
  )), error=function(e) e)

  message("Fitting checkBarrier model with psi...")
  par_psi <- list(psi = sim_pars$psi)
  fit_psi <- tryCatch(suppressMessages(fitLangevin(
    data = data_suggest,
    model = "underdamped",
    par = par_psi,
    spatialCovs = covs,
    barrier = "coast_barrier",
    lambda = barrier_spec$recommended,
    silent = TRUE
  )), error=function(e) e)

  if(!inherits(fit_filter, "error")){
    estUD_filter <- getUD(covs, fit=fit_filter, log=TRUE, plot=FALSE)
    parMat_filter[isim, 1:6] <- fit_filter$estimates$natural$Estimate[1:6]
    parMat_filter[isim, "BA"] <- rasterOverlap(exp(estUD_filter), exp(trueUD))#,local=fit_filter)
  }
  if(!inherits(fit_psi, "error")){
    estUD_psi <- getUD(covs, fit=fit_psi, log=TRUE, plot=FALSE)
    parMat_psi[isim, 1:6] <- fit_psi$estimates$natural$Estimate[1:6]
    parMat_psi[isim, 7] <- fit_psi$estimates$natural["psi",1]
    parMat_psi[isim, "BA"] <- rasterOverlap(exp(estUD_psi), exp(trueUD))#,local=fit_psi)
  }
  if(!inherits(fit_suggest, "error")){
    print(plot(fit_suggest,data=data_suggest,spatialCovs=covs,log=TRUE,maskBarrier=TRUE)+labs(title=paste0("Sim ",isim)))
    estUD_suggest <- getUD(covs, fit=fit_suggest, log=TRUE, plot=FALSE)
    parMat_suggest[isim, 1:6] <- fit_suggest$estimates$natural$Estimate[1:6]
    parMat_suggest[isim, "BA"] <- rasterOverlap(exp(estUD_suggest), exp(trueUD))#,local=fit_suggest)
  }

  message("\n--- Approach 1: filter all points on land ---")
  print(paste0("overall ", paste0(round(apply(parMat_filter[1:isim,,drop=FALSE], 2, mean, na.rm=TRUE), 4), collapse=" ")))

  message("--- Approach 2: suggested filter ---")
  print(paste0("overall ", paste0(round(apply(parMat_suggest[1:isim,,drop=FALSE], 2, mean, na.rm=TRUE), 4), collapse=" ")))

  message("--- Approach 3: suggested filter with psi ---")
  print(paste0("overall ", paste0(round(apply(parMat_psi[1:isim,,drop=FALSE], 2, mean, na.rm=TRUE), 4), collapse=" ")))

}
