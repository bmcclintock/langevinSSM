library(langevinSSM)
library(terra)
library(ggplot2)

sca <- 200
set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
timeStep <- 0.01
samplingRate <- 1
propMissing <- 0
obsPerAnimal <- 5000
covRange <- c(0.1,0.5)

measurementError = list(smaj.sd = 1.5, smin.sd = 0.75)

sim_pars <- list(
  beta = c(4, -1, -0.2 , -0.1), # Attraction to cov1, repulsion from cov2, attraction to center, attraction to shore
  sigma = 5,
  gamma = 2
)

nSims <- 100
parMat <- matrix(NA,nSims,7)
colnames(parMat) <- c(paste0("beta",1:length(sim_pars$beta)),"sigma","gamma","BA")
parMat_true <- parMat_min <- parMat_max <- parMat_nopen <- parMat
for(isim in 1:nSims){
  cat("Simulation",isim,"\n")

  sim_data <- tryCatch(stop(),error=function(e) e)
  while(inherits(sim_data,"error")){

    cov1 <- simCov(sca = sca, irange = runif(1,covRange[1],covRange[2]), sigma2 = 0.1, kappa = 0.5)
    cov2 <- simCov(sca = sca, irange = runif(1,covRange[1],covRange[2]), sigma2 = 0.1, kappa = 0.5)
    coords <- terra::crds(cov1)

    # create coast and islands
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

    dist2 <- ((coords[, "x"])^2 + coords[, "y"]^2) / sca # warped to contain y-axis
    d2c <- terra::setValues(cov1, dist2)

    names(d2c) <- "d2c"

    covs <- list(cov1 = cov1, cov2 = cov2, d2c = d2c, coast_barrier = water_mask)

    # ==========================================
    # 2. Simulation Parameters
    # ==========================================

    sim_data <- tryCatch({
      out_data <- simLangevin(
        model = "underdamped",
        par = sim_pars,
        spatialCovs = covs,
        nbAnimals = 3,
        obsPerAnimal = obsPerAnimal,
        timeStep = timeStep,
        barrier = "coast_barrier",
        measurementError = measurementError,
        subSample = list(samplingRate=samplingRate,propMissing=propMissing)
      )

      pts <- cbind(out_data$x, out_data$y)
      extracted_vals <- terra::extract(covs$coast_barrier, pts)

      barrier_vals <- extracted_vals[[ncol(extracted_vals)]]

      if (!any(barrier_vals <= 0, na.rm = TRUE)) {
        stop("No observed locations on land. Resimulating...")
      }

      out_data
    }, error=function(e) e)
  }

  trueUD <- getUD(covs,beta=sim_pars$beta,barrier="coast_barrier",lambda=attr(sim_data,"lambda"),log=TRUE,plot=FALSE)


  #plot(sim_data,beta=sim_pars$beta,spatialCovs=covs,log=TRUE,maskBarrier = TRUE)


  # ==========================================
  # 3. Model Fitting
  # ==========================================
  fit_uncon <- tryCatch(suppressMessages(fitLangevin(
    data = sim_data,
    model = "underdamped",
    spatialCovs = covs,
    map=list(beta=factor(c(1:(length(sim_pars$beta)-1),NA))),
    silent = TRUE
  )),error=function(e) e)
  #fit_uncon

  if(!inherits(fit_uncon,"error")){

    lambda_est <- suggestLambda(sim_data, fit_uncon, covs, "coast_barrier")

    fit_suggest <- tryCatch(fitLangevin(
      data = sim_data,
      model = "underdamped",
      spatialCovs = covs,
      barrier = "coast_barrier",
      lambda = lambda_est$recommended,
      silent = TRUE
    ),error=function(e) e)

    fit_nopen <- tryCatch(fitLangevin(
      data = sim_data,
      model = "underdamped",
      spatialCovs = covs,
      barrier = "coast_barrier",
      lambda = 0,
      silent = TRUE
    ),error=function(e) e)

    fit_min <- tryCatch(fitLangevin(
      data = sim_data,
      model = "underdamped",
      spatialCovs = covs,
      barrier = "coast_barrier",
      lambda = lambda_est$lambda_min,
      silent = TRUE
    ),error=function(e) e)

    fit_max <- tryCatch(fitLangevin(
      data = sim_data,
      model = "underdamped",
      spatialCovs = covs,
      barrier = "coast_barrier",
      lambda = lambda_est$lambda_max,
      silent = TRUE
    ),error=function(e) e)
  }

  fit_true <- tryCatch(fitLangevin(
    data = sim_data,
    model = "underdamped",
    spatialCovs = covs,
    barrier = "coast_barrier",
    silent = TRUE
  ),error=function(e) e)

  if(!inherits(fit_min,"error")){
    estUD_min <- getUD(covs,fit=fit_min,log=TRUE,plot=FALSE)
    parMat_min[isim,1:6] <- fit_min$estimates$natural$Estimate[1:6]
    parMat_min[isim,"BA"] <- rasterOverlap(exp(estUD_min),exp(trueUD))
  }
  if(!inherits(fit_max,"error")){
    estUD_max <- getUD(covs,fit=fit_max,log=TRUE,plot=FALSE)
    parMat_max[isim,1:6] <- fit_max$estimates$natural$Estimate[1:6]
    parMat_max[isim,"BA"] <- rasterOverlap(exp(estUD_max),exp(trueUD))
  }
  if(!inherits(fit_suggest,"error")){
    print(plot(fit_suggest,data=sim_data,spatialCovs=covs,log=TRUE,maskBarrier=TRUE)+labs(title=paste0("Sim ",isim)))
    estUD_suggest <- getUD(covs,fit=fit_suggest,log=TRUE,plot=FALSE)
    parMat[isim,1:6] <- fit_suggest$estimates$natural$Estimate[1:6]
    parMat[isim,"BA"] <- rasterOverlap(exp(estUD_suggest),exp(trueUD))
  }
  if(!inherits(fit_nopen,"error")){
    print(plot(fit_nopen,data=sim_data,spatialCovs=covs,log=TRUE,maskBarrier=TRUE)+labs(title=paste0("Sim ",isim)))
    estUD_nopen <- getUD(covs,fit=fit_nopen,log=TRUE,plot=FALSE)
    parMat_nopen[isim,1:6] <- fit_nopen$estimates$natural$Estimate[1:6]
    parMat_nopen[isim,"BA"] <- rasterOverlap(exp(estUD_nopen),exp(trueUD))
  }
  if(!inherits(fit_true,"error")){
    estUD_true <- getUD(covs,fit=fit_true,log=TRUE,plot=FALSE)
    parMat_true[isim,1:6] <- fit_true$estimates$natural$Estimate[1:6]
    parMat_true[isim,"BA"] <- rasterOverlap(exp(estUD_true),exp(trueUD))
  }
  message("suggest lamda")
  print(paste0("            ",paste0(colnames(parMat),collapse="    ")))
  print(paste0("current ",paste0(round(parMat[isim,],6),collapse=" ")))
  print(paste0("overall ",paste0(round(apply(parMat[1:isim,,drop=FALSE],2,mean,na.rm=TRUE),6),collapse=" ")))
  message("suggest lamda no penalty")
  print(paste0("            ",paste0(colnames(parMat_nopen),collapse="    ")))
  print(paste0("current ",paste0(round(parMat_nopen[isim,],6),collapse=" ")))
  print(paste0("overall ",paste0(round(apply(parMat_nopen[1:isim,,drop=FALSE],2,mean,na.rm=TRUE),6),collapse=" ")))
  message("min lambda")
  print(paste0("            ",paste0(colnames(parMat_min),collapse="    ")))
  print(paste0("current ",paste0(round(parMat_min[isim,],6),collapse=" ")))
  print(paste0("overall ",paste0(round(apply(parMat_min[1:isim,,drop=FALSE],2,mean,na.rm=TRUE),6),collapse=" ")))
  message("max lambda")
  print(paste0("            ",paste0(colnames(parMat_max),collapse="    ")))
  print(paste0("current ",paste0(round(parMat_max[isim,],6),collapse=" ")))
  print(paste0("overall ",paste0(round(apply(parMat_max[1:isim,,drop=FALSE],2,mean,na.rm=TRUE),6),collapse=" ")))
  message("true lambda")
  print(paste0("            ",paste0(colnames(parMat_true),collapse="    ")))
  print(paste0("current ",paste0(round(parMat_true[isim,],6),collapse=" ")))
  print(paste0("overall ",paste0(round(apply(parMat_true[1:isim,,drop=FALSE],2,mean,na.rm=TRUE),6),collapse=" ")))
}

# ==========================================
# 4. Diagnostics & Prediction
# ==========================================

# A. Calculate and Plot OSA Residuals
message("Calculating residuals...")
res <- residuals(fit, data = sim_data, spatialCovs = covs)
res_plots <- plot(res)
# View specific diagnostic (e.g., Mahalanobis distance Q-Q)
print(res_plots$qq_mah)

# B. Estimate Utilization Distribution from Fit
# This automatically retrieves the barrier/lambda from fit$conditions
message("Generating UD from fitted model...")
fitted_ud <- getUD(spatialCovs = covs, fit = fit)

# C. Simulate from the Fitted Model (Posterior Predictive Check)
message("Simulating from fitted parameters...")
sim_from_fit <- simLangevin(
  model = fit,
  data = sim_data,
  spatialCovs = covs,
  conditional = FALSE # Generate a new forward track from the start point
)

# ==========================================
# 5. Comparative Visualization
# ==========================================
message("Plotting results...")

# Overlay original track and simulation from fit on the UD
plot(fit, spatialCovs = covs, data = sim_data) +
  labs(title = "Fitted UD with Original and Predicted Tracks")

# Plot standard error of the UD (Delta Method)
# getUD returns a stack; index 2 and 3 are typically SE and CV
plotUD(fitted_ud)
