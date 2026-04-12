library(langevinSSM)


n_sims <- 100

beta <- c(-4, 6, 5, -0.1) # resource selection coefficients for the spatial covariates (cov_1, cov_2, ... cov_ncov, d2c)
ncov <- length(beta) - 1 # number of spatial covariates to be generated using simCov
sigma <- 5 # diffusion (or speed) parameter
gamma <- 0.5 # friction parameter (smaller value -> more directional persistence); ignored unless model=="underdamped"

smaj.sd <- 1.5 # SD for semi-major error ellipse axis; smaj ~ abs(Normal(0,smaj.sd))
smin.sd <- smaj.sd/2 # SD for semi-minor error ellipse axis; smin ~ abs(Normal(0,smin.sd))
eor.lim <- c(0,180) # range for error ellipse orientation (in degrees from north); eor ~ Uniform(eor.lim[1],eor.lim[2])
measurementError <- list(smaj.sd=smaj.sd,smin.sd=smin.sd,eor.lim=eor.lim) # setting measurementError <- NULL adds no measurement error to observations

sca <- 200 # bounding box scale
covRange <- c(0.1,0.5) # lower and upper bounds for covariate spatial range parameter (lower has less spatial autocorrelation)

thresh_percent = 0.99

captured <- matrix(NA, nrow = n_sims, ncol = 4, dimnames = list(NULL, c("delta_true", "sim_true", "delta_est", "sim_est")))
spatialCovs <- vector("list", n_sims)

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")
for(isim in 1:n_sims) {

  covNames <- c("cov1","cov2","cov3","d2c")
  spatialCovs[[isim]] <- list()

  for(i in 1:ncov) {

    spatialCovs[[isim]][[i]] <- list()

    irange <- runif(1,covRange[1],covRange[2])
    spatialCovs[[isim]][[i]] <- simCov(sca = sca, irange=irange, sigma2 = 0.1, kappa = 0.5)
    names(spatialCovs[[isim]][[i]]) <- paste0("cov",i)

  }

  coords <- terra::crds(spatialCovs[[isim]][[1]])
  dist2 <- (coords[, "x"]^2 + coords[, "y"]^2) / sca

  spatialCovs[[isim]][[4]] <- terra::setValues(spatialCovs[[isim]][[1]][[1]], dist2)
  names(spatialCovs[[isim]][[4]]) <- covNames[4]
  names(spatialCovs[[isim]]) <- covNames

  true_ud_raster <- getUD(spatialCovs[[isim]], beta = beta, log=FALSE)

  true_vals <- terra::values(true_ud_raster)

  sim_data <- simLangevin(par = list(beta=beta, sigma=sigma, gamma=gamma), obsPerAnimal = 5000, nbAnimals = 3, subSample = list(samplingRate = 10),
                          spatialCovs = spatialCovs[[isim]], measurementError = measurementError)

  fit <- fitLangevin(data = sim_data,
                           spatialCovs = spatialCovs[[isim]],
                           silent = TRUE)

  ud_stack <- getUD(spatialCovs[[isim]], fit, nSims = 4000, plot = FALSE, log = FALSE)

  est_vals <- terra::values(ud_stack[["UD"]])

  delta_cv <- terra::values(ud_stack[["UD_CV_delta"]])
  sim_cv   <- terra::values(ud_stack[["UD_CV_sim"]])

  delta_lwr <- est_vals * exp(-1.96 * delta_cv)
  delta_upr <- est_vals * exp( 1.96 * delta_cv)

  sim_lwr <- est_vals * exp(-1.96 * sim_cv)
  sim_upr <- est_vals * exp( 1.96 * sim_cv)

  delta_upr <- pmin(delta_upr, 1)
  sim_upr   <- pmin(sim_upr, 1)

  delta_captured <- (true_vals >= delta_lwr) & (true_vals <= delta_upr)
  sim_captured   <- (true_vals >= sim_lwr)   & (true_vals <= sim_upr)

  est_vals_vec <- as.numeric(est_vals)

  sorted_est_vals <- sort(est_vals_vec, decreasing = TRUE)

  cum_probs_est <- cumsum(sorted_est_vals)

  threshold_est <- sorted_est_vals[min(which(cum_probs_est >= thresh_percent))]

  meaningful_pixels <- est_vals >= threshold_est

  # calculate the mean coverage only within the estimated thresh_percent% UD contour
  captured[isim, "delta_est"] <- mean(delta_captured[meaningful_pixels], na.rm = TRUE)
  captured[isim, "sim_est"]   <- mean(sim_captured[meaningful_pixels], na.rm = TRUE)

  true_vals_vec <- as.numeric(true_vals)

  sorted_vals <- sort(true_vals_vec, decreasing = TRUE)

  cum_probs <- cumsum(sorted_vals)

  threshold <- sorted_vals[min(which(cum_probs >= thresh_percent))]

  meaningful_pixels_true <- true_vals >= threshold

  # calculate the mean coverage only within the true thresh_percent% UD contour
  captured[isim, "delta_true"] <- mean(delta_captured[meaningful_pixels_true], na.rm = TRUE)
  captured[isim, "sim_true"]   <- mean(sim_captured[meaningful_pixels_true], na.rm = TRUE)

  print(paste0("Simulation ",isim," est  (",sum(meaningful_pixels)," cells): ",paste0(round(apply(captured[1:isim,c("delta_est","sim_est"),drop=FALSE ], 2, mean),3), collapse = " | ")))
  print(paste0("Simulation ",isim," true (",sum(meaningful_pixels_true)," cells): ",paste0(round(apply(captured[1:isim,c("delta_true","sim_true"),drop=FALSE ], 2, mean),3), collapse = " | ")))

}
