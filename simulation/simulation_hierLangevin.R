library(langevinSSM)
library(terra)

# ==============================================================================
# 1. SETUP & POPULATION PARAMETERS
# ==============================================================================
set.seed(2026)
ncores <- 5 # Number of cores for parallel processing (fitLangevin uses parallel::mclapply)
mpl <- 1.5  # Maximum Penalized Likelihood (MPL) factor for hierLangevin

# Simulation settings
nSims <- 100          # Number of simulations to run
n_ind <- 30          # Number of individuals to simulate per simulation
n_obs <- 60000        # Number of observations per individual

sigma <- 5 # true population scale parameter
gamma <- 0.5 # true population friction parameter

# Define true population means (working scale)
true_mu_beta  <- c(-4, 6, 5, -0.1)
ncov <- length(true_mu_beta) - 1 # number of spatial covariates to be generated using simCov
true_mu_log_sigma <- log(sigma)
true_mu_log_gamma <- log(gamma)

# Define true between-individual standard deviations
true_sd_beta  <- c(0.4, 0.6, 0.5, 0.01)
true_sd_log_sigma <- 0.2
true_sd_log_gamma <- 0.15

dt    <- 0.01        # Time step
samplingRate <- 2
scaleFactor <- 10

Dt <- dt * samplingRate # Time step after subsampling (i.e. the time between observations in the final dataset)
SNR <- 2 # signal-to-noise ratio for measurement error
sigma_err <- sqrt((sigma^2/gamma^3 * (gamma*Dt - 2*(1-exp(-gamma*Dt))+1/2*(1-exp(-2*gamma*Dt))))/SNR)
measurementError <- list(x.sd = sigma_err, y.sd = sigma_err) # location measurement error

## specify scale and spatial autocorrelation for covariates
sca <- 200 # bounding box scale
covRange <- c(0.1,0.5) # lower and upper bounds for covariate spatial range parameter (lower has less spatial autocorrelation)
covTimes <- 1 # number of time-varying covariate layers (if >1, then each covariate is a SpatRaster with covTimes layers; if 1, then each covariate is a single-layer SpatRaster)
ztimes <- seq(0,n_obs*dt,length=covTimes) # z (i.e. time) values for raster stacks

# Combine true generative parameters for easier tracking and alignment
true_mu <- c(true_mu_beta, true_mu_log_sigma, true_mu_log_gamma)
true_sd <- c(true_sd_beta, true_sd_log_sigma, true_sd_log_gamma)
p <- length(true_mu)

# Storage arrays for estimates across all nSims iterations
est_mu_matrix <- matrix(NA, nrow = nSims, ncol = p)
est_sd_matrix <- matrix(NA, nrow = nSims, ncol = p)
param_names <- NULL

# ==============================================================================
# 2. RUN SIMULATION LOOP
# ==============================================================================
for (sim in 1:nSims) {
  cat("\n\n======================================================\n")
  cat("STARTING SIMULATION", sim, "OF", nSims, "\n")
  cat("======================================================\n")

  # ----------------------------------------------------------------------------
  # 2a. GENERATE COVARIATES
  # ----------------------------------------------------------------------------
  message("Generating spatial covariates...")
  covNames <- c("cov1","cov2","cov3","d2c")
  covs <- list()

  for(i in 1:ncov) {
    covs[[i]] <- list()
    if(covTimes>1){
      for(ztime in 1:length(ztimes)){
        irange <- runif(1,covRange[1],covRange[2])
        covs[[i]][[ztime]] <- simCov(sca = sca, irange=irange, sigma2 = 0.1, kappa = 0.5)
        terra::time(covs[[i]][[ztime]]) <- ztimes[ztime]
        names(covs[[i]][[ztime]]) <- covNames[i]
      }
      covs[[i]] <- terra::rast(covs[[i]])
    } else {
      irange <- runif(1,covRange[1],covRange[2])
      covs[[i]] <- simCov(sca = sca, irange=irange, sigma2 = 0.1, kappa = 0.5)
      names(covs[[i]]) <- paste0("cov",i)
    }
  }

  coords <- terra::crds(covs[[1]])
  dist2 <- (coords[, "x"]^2 + coords[, "y"]^2) / sca

  covs[[4]] <- terra::setValues(covs[[1]][[1]], dist2)
  names(covs[[4]]) <- covNames[4]
  names(covs) <- covNames

  # ----------------------------------------------------------------------------
  # 2b. SIMULATE INDIVIDUAL TRACKS
  # ----------------------------------------------------------------------------
  message("Simulating tracks for ", n_ind, " individuals...")
  sim_data_list <- list()

  for (i in seq_len(n_ind)) {
    # Draw individual working-scale parameters from the population distribution
    ind_beta      <- rnorm(length(true_mu_beta), true_mu_beta, true_sd_beta)
    ind_log_sigma <- rnorm(1, true_mu_log_sigma, true_sd_log_sigma)
    ind_log_gamma <- rnorm(1, true_mu_log_gamma, true_sd_log_gamma)

    # Format parameters for simLangevin (requires natural scale for sigma/gamma)
    par_i <- list(
      beta  = ind_beta,
      sigma = exp(ind_log_sigma),
      gamma = exp(ind_log_gamma)
    )

    # We add a small amount of GPS-style measurement error to make the fits realistic
    sim_track <- suppressMessages(simLangevin(
      model = "underdamped",
      spatialCovs = covs,
      par = par_i,
      nbAnimals = 1,
      obsPerAnimal = n_obs,
      timeStep = dt,
      measurementError = measurementError,
      subSample = list(samplingRate=samplingRate)
    ))

    # Assign ID as a factor to prevent TMB errors downstream
    sim_track$id <- as.factor(paste0("ind_", i))
    sim_data_list[[i]] <- sim_track
  }

  # ----------------------------------------------------------------------------
  # 2c. FIT INDIVIDUAL MODELS (STAGE I)
  # ----------------------------------------------------------------------------
  message("Fitting Stage I independent models...")

  # NOTE: In a real-world scenario with heavy data, this lapply loop is easily
  # replaced with parallel::mclapply (Mac/Linux) or future.apply::future_lapply
  # to run all individuals concurrently on multiple cores.
  fit_list <- parallel::mclapply(seq_along(sim_data_list), function(i) {
    # Wrap in try() to ensure one failed individual doesn't break the pipeline
    try({
      suppressMessages(fitLangevin(
        data = sim_data_list[[i]],
        model = "underdamped",
        spatialCovs = covs,
        scaleFactor = scaleFactor,
        silent = TRUE
      ))
    }, silent = TRUE)
  }, mc.cores = ncores)

  # Name the list elements by their individual IDs (hierLangevin uses these)
  names(fit_list) <- paste0("ind_", seq_len(n_ind))

  # ----------------------------------------------------------------------------
  # 2d. HIERARCHICAL POOLING (STAGE II)
  # ----------------------------------------------------------------------------
  message("Fitting Stage II hierarchical model...")

  # Pass the list of fitLangevin objects directly to hierLangevin
  # We use mpl = 1.5 (Maximum Penalized Likelihood) to guard against
  # zero-variance boundary collapse, as recommended for this model
  stage2_fit <- try({
    hierLangevin(fit_list, mpl = mpl, silent = TRUE)
  }, silent = TRUE)

  # ----------------------------------------------------------------------------
  # 2e. STORE RESULTS & CALCULATE RUNNING AVERAGE
  # ----------------------------------------------------------------------------
  if (inherits(stage2_fit, "hierLangevin")) {
    # Dynamically extract parameter names from the first successful fit
    if (is.null(param_names)) {
      param_names <- rownames(stage2_fit$estimates$working)[1:p]
    }

    # Store parameter estimates for this simulation
    est_mu_matrix[sim, ] <- stage2_fit$estimates$working$Estimate[1:p]
    est_sd_matrix[sim, ] <- stage2_fit$estimates$working$Estimate[(p+1):(2*p)]
  } else {
    message("  -> WARNING: Stage II hierarchical fit failed for simulation ", sim)
  }

  # Calculate running averages of population means and standard deviations
  run_avg_mu <- colMeans(est_mu_matrix[1:sim, , drop = FALSE], na.rm = TRUE)
  run_avg_sd <- colMeans(est_sd_matrix[1:sim, , drop = FALSE], na.rm = TRUE)

  # Print running average to the console
  cat("\n--- Running Average (Working Scale) After", sim, "Simulations ---\n")
  run_avg_df <- data.frame(
    Parameter = if (!is.null(param_names)) param_names else paste0("Param_", 1:p),
    True_Mean = true_mu,
    RunAvg_Mean = run_avg_mu,
    True_SD   = true_sd,
    RunAvg_SD = run_avg_sd
  )
  print(run_avg_df, row.names = FALSE)
}

# ==============================================================================
# 3. FINAL SUMMARY
# ==============================================================================
cat("\n\n======================================================\n")
cat("FINAL SIMULATION SUMMARY ACROSS", nSims, "ITERATIONS\n")
cat("======================================================\n")
print(run_avg_df, row.names = FALSE)
