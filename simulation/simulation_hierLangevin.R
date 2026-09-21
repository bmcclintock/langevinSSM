library(langevinSSM)
library(terra)
library(MASS)

# ==============================================================================
# 1. SETUP & POPULATION PARAMETERS
# ==============================================================================
set.seed(2026)
ncores <- 5 # Number of cores for parallel processing
mpl <- 1.5  # Maximum Penalized Likelihood (MPL) factor for hierLangevin

# Simulation settings
nSims <- 100          # Number of simulations to run
n_ind <- 30           # Number of individuals to simulate per simulation
n_obs <- 60000       # Number of observations per individual

sigma <- 5   # true population scale parameter
gamma <- 0.5 # true population friction parameter

# Define true population means (working scale)
true_mu_beta  <- c(0, 0, 0, -0.1)
ncov <- length(true_mu_beta) - 1
true_mu_log_gamma <- log(gamma)
true_mu_log_sigma <- log(sigma)

# Define true between-individual standard deviations
true_sd_beta  <- c(4,4,4, 0.01)
true_sd_log_sigma <- 0.2
true_sd_log_gamma <- 0.15

dt    <- 0.01        # Time step
samplingRate <- 2
scaleFactor <- 10

Dt <- dt * samplingRate
SNR <- 2
sigma_err <- sqrt((sigma^2/gamma^3 * (gamma*Dt - 2*(1-exp(-gamma*Dt))+1/2*(1-exp(-2*gamma*Dt))))/SNR)
measurementError <- list(x.sd = sigma_err, y.sd = sigma_err)

## specify scale and spatial autocorrelation for covariates
sca <- 200
covRange <- c(0.1,0.5)

# Combine true generative parameters for easier tracking and alignment
true_mu <- c(true_mu_beta, true_mu_log_sigma, true_mu_log_gamma)
true_sd <- c(true_sd_beta, true_sd_log_sigma, true_sd_log_gamma)
p <- length(true_mu)

# Storage arrays
overlap_emp <- numeric(nSims)
overlap_true <- numeric(nSims)
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
    irange <- runif(1,covRange[1],covRange[2])
    covs[[i]] <- simCov(sca = sca, irange=irange, sigma2 = 0.1, kappa = 0.5)
    names(covs[[i]]) <- paste0("cov",i)
  }

  coords <- terra::crds(covs[[1]])
  dist2 <- (coords[, "x"]^2 + coords[, "y"]^2) / sca

  covs[[4]] <- terra::setValues(covs[[1]], dist2)
  names(covs[[4]]) <- covNames[4]
  names(covs) <- covNames

  # ----------------------------------------------------------------------------
  # 2b. SIMULATE INDIVIDUAL TRACKS & BUILD EMPIRICAL POPULATION UD
  # ----------------------------------------------------------------------------
  message("Simulating tracks for ", n_ind, " individuals...")
  sim_data_list <- list()

  # Initialize an empty raster to tally all true simulated locations
  emp_UD <- terra::setValues(covs[[1]], 0)
  names(emp_UD) <- "UD"
  all_true_x <- numeric()
  all_true_y <- numeric()

  for (i in seq_len(n_ind)) {
    ind_beta      <- rnorm(length(true_mu_beta), true_mu_beta, true_sd_beta)
    # Generate original scale sigma for simLangevin
    ind_log_sigma <- rnorm(1, true_mu_log_sigma, true_sd_log_sigma)
    ind_log_gamma <- rnorm(1, true_mu_log_gamma, true_sd_log_gamma)

    par_i <- list(
      beta  = ind_beta,
      sigma = exp(ind_log_sigma),
      gamma = exp(ind_log_gamma)
    )

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

    # Store true continuous-time generative locations (before measurement error)
    all_true_x <- c(all_true_x, sim_track$mu.x)
    all_true_y <- c(all_true_y, sim_track$mu.y)

    sim_track$id <- as.factor(paste0("ind_", i))
    sim_data_list[[i]] <- sim_track
  }

  # Calculate the true empirical population UD by tallying cell hits
  message("Calculating True Empirical Population UD...")
  cells <- terra::cellFromXY(emp_UD, cbind(all_true_x, all_true_y))
  cells <- cells[!is.na(cells)]
  cell_counts <- table(cells)
  emp_UD[as.numeric(names(cell_counts))] <- as.numeric(cell_counts)

  # Normalize to create a valid probability distribution
  emp_UD <- emp_UD / terra::global(emp_UD, "sum", na.rm=TRUE)[1,1]

  # ----------------------------------------------------------------------------
  # 2c. BUILD TRUE THEORETICAL POPULATION UD
  # ----------------------------------------------------------------------------
  message("Calculating True Theoretical Population UD via Monte Carlo integration...")

  n_mc_true <- 1000
  true_cov_mat <- diag(true_sd_beta^2)
  true_beta_draws <- MASS::mvrnorm(n_mc_true, true_mu_beta, true_cov_mat)

  cov_mats <- lapply(covs, terra::as.matrix, wide = FALSE)
  n_cells <- terra::ncell(covs[[1]])
  true_pop_UD_vals <- numeric(n_cells)

  for (m in 1:n_mc_true) {
    b <- true_beta_draws[m, ]
    log_pi <- cov_mats[[1]] * b[1]
    for (j in 2:length(covs)) {
      log_pi <- log_pi + (cov_mats[[j]] * b[j])
    }
    pi_vals <- exp(log_pi - max(log_pi, na.rm = TRUE))
    pi_vals <- pi_vals / sum(pi_vals, na.rm = TRUE)
    true_pop_UD_vals <- true_pop_UD_vals + (pi_vals / n_mc_true)
  }

  true_UD <- terra::setValues(covs[[1]], true_pop_UD_vals)
  names(true_UD) <- "UD"

  # ----------------------------------------------------------------------------
  # 2d. FIT INDIVIDUAL MODELS (STAGE I)
  # ----------------------------------------------------------------------------
  message("Fitting Stage I independent models...")

  fit_list <- parallel::mclapply(seq_along(sim_data_list), function(i) {
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

  names(fit_list) <- paste0("ind_", seq_len(n_ind))

  # ----------------------------------------------------------------------------
  # 2e. HIERARCHICAL POOLING (STAGE II)
  # ----------------------------------------------------------------------------
  message("Fitting Stage II hierarchical model...")

  stage2_fit <- try({
    hierLangevin(fit_list, mpl = mpl, silent = TRUE)
  }, silent = TRUE)

  # ----------------------------------------------------------------------------
  # 2f. STORE ESTIMATES & CALCULATE OVERLAPS
  # ----------------------------------------------------------------------------
  if (inherits(stage2_fit, "hierLangevin")) {

    # Store parameter estimates
    if (is.null(param_names)) {
      param_names <- rownames(stage2_fit$estimates$working)[1:p]
    }
    est_mu_matrix[sim, ] <- stage2_fit$estimates$working$Estimate[1:p]
    est_sd_matrix[sim, ] <- stage2_fit$estimates$working$Estimate[(p+1):(2*p)]

    # Calculate UD Overlaps
    message("Estimating Population UD from hierLangevin model...")
    est_UD_stack <- getUD(spatialCovs = covs,
                          fit = stage2_fit,
                          nSims = 1000,
                          log = FALSE,
                          plot = FALSE)

    overlap_emp[sim] <- rasterOverlap(r1 = est_UD_stack$UD_Population, r2 = emp_UD)
    overlap_true[sim] <- rasterOverlap(r1 = est_UD_stack$UD_Population, r2 = true_UD)

  } else {
    message("  -> WARNING: Stage II hierarchical fit failed for simulation ", sim)
  }

  # Calculate running averages of population parameters
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

  # Print UD Overlaps
  cat("\n--- Pop-level Utilization Distribution Overlap (Bhattacharyya's Affinity) ---\n")
  cat(sprintf("Estimated vs Empirical Tally: %.4f (Running Avg: %.4f)\n",
              overlap_emp[sim], mean(overlap_emp[1:sim], na.rm = TRUE)))
  cat(sprintf("Estimated vs True Theoretical: %.4f (Running Avg: %.4f)\n",
              overlap_true[sim], mean(overlap_true[1:sim], na.rm = TRUE)))
}

# ==============================================================================
# 3. FINAL SUMMARY
# ==============================================================================
cat("\n\n======================================================\n")
cat("FINAL SIMULATION SUMMARY ACROSS", nSims, "ITERATIONS\n")
cat("======================================================\n")
print(run_avg_df, row.names = FALSE)

cat("\n--- Utilization Distribution Overlap Summary ---\n")
cat(sprintf("Avg Overlap w/ Empirical Tally: %.4f\n", mean(overlap_emp, na.rm = TRUE)))
cat(sprintf("Avg Overlap w/ True Theoretical: %.4f\n", mean(overlap_true, na.rm = TRUE)))
