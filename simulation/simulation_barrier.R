library(langevinSSM)
library(terra)
library(ggplot2)

model <- "underdamped"

sca <- 600
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")
timeStep <- 0.01
samplingRate <- 1
propMissing <- 0
nbAnimals <- 5
obsPerAnimal <- 5000
covRange <- c(0.1, 0.5)
lambda_max <- NULL

obs_dist <- "t"
df <- 3 # degrees of freedom (applies to obs_dist <- "t" only)
k <- 2 # threshold (applied to obs_dist <- "Huber" only)

measurementError = list(smaj.sd = 1.5, smin.sd = 0.75)

sim_pars <- list(
  beta = c(4, -1, -0.2 , -0.1),
  sigma = 5,
  psi = 1
)
if(obs_dist=="t") sim_pars$df <- df
if(obs_dist=="Huber") sim_pars$k <- k
if(model=="underdamped") sim_pars$gamma <- 0.5

nSims <- 100
parMat <- matrix(NA, nSims, 8+ifelse(obs_dist=="Gaussian",0,1))
colnames(parMat) <- c(paste0("beta", 1:length(sim_pars$beta)), "sigma", "gamma", paste0(ifelse(obs_dist=="Gaussian","",ifelse(obs_dist=="t","df","k"))),"lambda", "BA")
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

    barrier <- prepBarrier(water_mask)

    covs <- list(cov1 = cov1, cov2 = cov2, d2c = d2c, d2coast = barrier/100, coast_barrier = barrier)

    sim_data <- tryCatch({
      out_data <- simLangevin(
        model = model,
        obs_dist = obs_dist,
        par = sim_pars,
        spatialCovs = covs,
        barrier = "coast_barrier",
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

  trueUD <- getUD(covs, beta=sim_pars$beta, barrier="coast_barrier", lambda=attr(sim_data,"lambda"), log=TRUE, plot=FALSE)


  par <- initialValues(sim_data,model=model,obs_dist = obs_dist,spatialCovs=covs[1:length(sim_pars$beta)],lambda="estimate")
  if(obs_dist=="t") par$df <- 5
  if(obs_dist=="Huber") par$k <- 5

  fit <- tryCatch(suppressMessages(fitLangevin(
    data = sim_data,
    model = model,
    obs_dist = obs_dist,
    par = par,
    spatialCovs = covs,
    lambda = "estimate",
    silent = TRUE
  )),error=function(e) e)

  par <- initialValues(sim_data,model=model,obs_dist = obs_dist,spatialCovs=covs[1:length(sim_pars$beta)],lambda="tied")
  if(obs_dist=="t") par$df <- df
  if(obs_dist=="Huber") par$k <- k

  fit_true <- tryCatch(suppressMessages(fitLangevin(
    data = sim_data,
    model = model,
    obs_dist = obs_dist,
    par = par,
    spatialCovs = covs,
    silent = TRUE
  )),error=function(e) e)

  if(!inherits(fit, "error")){
    estUD <- suppressMessages(getUD(covs, fit=fit, log=TRUE, plot=FALSE))
    parMat[isim, 1:6] <- fit$estimates$natural$Estimate[1:6]
    if(obs_dist=="t") parMat[isim,"df"] <- fit$estimates$natural["df",1]
    if(obs_dist=="Huber") parMat[isim,"k"] <- fit$estimates$natural["k",1]
    parMat[isim, "BA"] <- rasterOverlap(exp(estUD), exp(trueUD))#,local=fit)
    parMat[isim, "lambda"] <- fit$conditions$lambda
    print(plot(fit,data=sim_data,spatialCovs=covs,log=TRUE,maskRast=water_mask)+labs(title=paste0("Sim ",isim)))
  }

  if(!inherits(fit_true, "error")){
    estUD_true <- suppressMessages(getUD(covs, fit=fit_true, log=TRUE, plot=FALSE))
    parMat_true[isim, 1:6] <- fit_true$estimates$natural$Estimate[1:6]
    if(obs_dist=="t") parMat_true[isim,"df"] <- fit_true$estimates$natural["df",1]
    if(obs_dist=="Huber") parMat_true[isim,"k"] <- fit_true$estimates$natural["k",1]
    parMat_true[isim, "BA"] <- rasterOverlap(exp(estUD_true), exp(trueUD))#,local=fit_true)
    parMat_true[isim, "lambda"] <- fit_true$conditions$lambda
  }

  message("\n Estimated lambda ")
  print(paste0("            ",paste0(colnames(parMat),collapse="    ")))
  print(paste0("current ",paste0(round(parMat[isim,],6),collapse=" ")))
  print(paste0("overall ", paste0(round(apply(parMat[1:isim,,drop=FALSE], 2, mean, na.rm=TRUE), 6), collapse=" ")))

  message(" Tied lambda ")
  print(paste0("            ",paste0(colnames(parMat_true),collapse="    ")))
  print(paste0("current ",paste0(round(parMat_true[isim,],6),collapse=" ")))
  print(paste0("overall ", paste0(round(apply(parMat_true[1:isim,,drop=FALSE], 2, mean, na.rm=TRUE), 6), collapse=" ")))

}
