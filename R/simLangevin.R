#' Simulate trajectories from the habitat-driven Langevin diffusion
#'
#' This function acts as a generic interface for simulating from the Langevin diffusion model. It handles two distinct workflows:
#' 1) Simulating an entirely new track "from scratch" by providing specific parameters (the default method).
#' 2) Simulating tracks from a previously fitted model (\code{fitLangevin} object) for further analysis (e.g. posterior predictive checks or imputation).
#'
#' @param model Can be either a character string specifying the model to simulate from scratch (\code{"underdamped"} or \code{"overdamped"}), or a fitted Langevin model object of class \code{fitLangevin} (as returned by \code{\link{fitLangevin}}).
#' @param ... Additional arguments passed to the specific methods.
#'
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates. The covariates must be on the same spatial grid and have the same spatial extent.
#' @param barrier Character string. The name of the barrier in \code{spatialCovs} that is represented as a signed distance field (see \code{\link{prepBarrier}}). If provided, this raster is exclusively used for the barrier penalty and is not included in the habitat selection covariates. Default: \code{NULL} (no barrier).
#'
#' @param par List of parameters. For the "underdamped" model, this must include \code{beta} (a vector of length equal to the number of covariates), \code{sigma} (speed parameter), and \code{gamma} (friction parameter). For the "overdamped" model, this must include \code{beta} and \code{sigma}. If \code{measurementError} is specified, optional observation process parameters include \code{psi}, \code{tau}, and \code{rho_o}. See \code{\link{fitLangevin}}.
#' @param lambda Numeric. The penalty weight for the barrier constraint. Larger values create a steeper "wall" preventing locations from crossing into restricted areas. Default: \code{NULL}. Because the true generative parameters are known during simulation, leaving this as \code{NULL} allows the function to perfectly auto-calculate the optimal theoretical stability limit based on the \code{model} type, true speed parameter (\code{sigma}), and \code{timeStep}. See Details.
#' @param timeStep Time step to use for the simulation. Determines the resolution of the discrete-time approximation of the continuous-time process. The smaller the \code{timeStep}, the more accurate the approximation. Ignored if \code{model} is a \code{fitLangevin} object and \code{conditional=TRUE}. Default: 0.01.
#' @param nbAnimals Number of animals to simulate. Default: 1.
#' @param obsPerAnimal Number of observations to simulate per animal. Default: 500.
#' @param initialPosition Initial position(s) for the simulation. A 2-vector, or a list of length \code{nbAnimals} of 2-vectors. If missing, initial positions are randomly generated based on the utilization distribution.
#' @param measurementError A list or data frame of parameters used to simulate observation error. See \code{\link{addMeasurementError}}. Default: \code{NULL}.
#' @param subSample List of specifications for subsampling data from the continuous-time process model, which can include \code{samplingRate} and \code{propMissing}. See \code{\link{subSampleData}}. Default: \code{NULL} (no subsampling or missing data).
#'
#' @param data The \code{dataLangevin} object (as returned by \code{\link{formatData}} or \code{\link{simLangevin}}) used to fit the model.
#' @param jointPrecision Logical. If \code{TRUE}, draws parameters and latent variables from the full joint precision matrix of both fixed parameters and random effects. If \code{FALSE}, fixes the movement parameters at their point estimates and draws only random effects. Default: \code{FALSE}.
#' @param conditional Logical. If \code{TRUE}, simulates tracks conditional on the observed data (imputation), adding observation error to the drawn latent states based on the measurement error information in \code{data}. If \code{FALSE}, simulates an entirely new track forwards in time (e.g. for posterior predictive check). Default: \code{FALSE}.
#'
#' @details
#' \strong{Simulating from Scratch (Default Method):}
#' When \code{model} is a character string (\code{"underdamped"} or \code{"overdamped"}), the function generates a completely new dataset. This requires specifying the movement parameters (\code{par}), the number of animals, and the desired measurement error structure.
#'
#' \strong{Simulating from a Fitted Model (fitLangevin Method):}
#' When \code{model} is a \code{fitLangevin} object, the function behaves as a diagnostic and simulation tool.
#' The \code{conditional} and \code{jointPrecision} arguments define four possible ways to simulate from the fitted model:
#' \itemize{
#'   \item \strong{\code{conditional = TRUE, jointPrecision = TRUE}:} Imputes tracks tied to the \code{data}, drawn from the full joint covariance matrix to account for uncertainty. \code{timeStep} is ignored.
#'   \item \strong{\code{conditional = TRUE, jointPrecision = FALSE}:} Imputes tracks tied to the \code{data}, drawn from the random effects covariance matrix with the movement parameters fixed at their point estimates. \code{timeStep} is ignored.
#'   \item \strong{\code{conditional = FALSE, jointPrecision = TRUE}:} Starting at the initial location for each track, generates unconstrained tracks forward in time using parameters drawn from the full joint covariance matrix
#'   \item \strong{\code{conditional = FALSE, jointPrecision = FALSE}:} Starting at the initial location for each track, generates unconstrained tracks forward in time using parameters drawn from the random effects covariance matrix with the movement parameters fixed at their point estimates.
#' }
#'
#' @template barrier_details
#'
#' @details
#' \strong{Barrier Penalty Auto-Scaling:}
#' When \code{lambda = NULL} and a \code{barrier} is provided, \code{simLangevin} automatically calculates the maximum theoretical stability limit for the barrier penalty based on the SDE numerical integration limits. Because the true generative speed parameter (\eqn{\sigma}) and the maximum simulation time step (\eqn{\max(\Delta t)}) are known, the optimal spring constant can be deterministically calculated to create the "hardest" possible boundary that will not cause the numerical solver to explode.
#'
#' For the \strong{overdamped} model, the restoring force acts directly on the animal's position. The stability ceiling scales linearly with the inverse of the time step:
#' \deqn{\lambda = \frac{2}{\sigma^2 \max(\Delta t)}}
#'
#' For the \strong{underdamped} model, the restoring force acts on the animal's velocity, creating a true harmonic oscillator. The stability ceiling scales with the inverse square of the time step, allowing for much stiffer penalties:
#' \deqn{\lambda = \frac{1}{\sigma^2 \max(\Delta t)^2}}
#'
#' This auto-calculated value is applied during simulation and stored as a \code{lambda} attribute in the returned \code{dataLangevin} object.
#'
#' @return A data frame of class \code{dataLangevin} containing the simulated trajectories.
#'
#' @examples
#' # underdamped model with measurement error
#'
#' par <- list(beta = c(-4, 6, 5, -0.1),
#'             sigma = 5,
#'             gamma = 0.5)
#'
#' # error ellipse
#' # exampleCovs included in package; see ?exampleCovs for details
#' simDat_ee <- simLangevin(par = par,
#'                       spatialCovs = exampleCovs,
#'                       measurementError = list(smaj.sd = 1.5,
#'                                               smin.sd = 0.75,
#'                                               eor.lim = c(0,180)))
#' # x- and y-axis errors
#' simDat_xy <- simLangevin(par = par,
#'                       spatialCovs = exampleCovs,
#'                       measurementError = list(x.sd = 1.5,
#'                                               y.sd = 1.5))
#'
#' # location quality classes based on provided probabilities
#' emf_df <- getEMF()
#'
#' simDat_lc <- simLangevin(par = par,
#'                          spatialCovs = exampleCovs,
#'                          measurementError = emf_df)
#'
#' \dontrun{
#' # simulate from fitted model
#' fit <- fitLangevin(data = exampleDat,
#'                    spatialCovs = exampleCovs,
#'                    silent = TRUE)
#'
#' # unconditional
#' simFit <- simLangevin(fit,
#'                       data = exampleDat,
#'                       spatialCovs = exampleCovs,
#'                       timeStep = "10 secs")
#'
#' # conditional on observed tracks
#' simFit_cond <- simLangevin(fit,
#'                            data = exampleDat,
#'                            spatialCovs = exampleCovs,
#'                            timeStep = "10 secs",
#'                            conditional = TRUE)
#'
#' # simulating with a barrier
#' # create a dummy barrier mask (left half restricted = 0, right half allowed = 1)
#' coast_barrier <- exampleCovs[[1]]
#' terra::values(coast_barrier) <- ifelse(terra::crds(coast_barrier)[, "x"]
#'                                  >= mean(terra::crds(coast_barrier)[, "x"]), 1, 0)
#' names(coast_barrier) <- "coast_barrier"
#'
#' # convert mask to SDF and add to the spatial covariates list
#' exampleCovs_barrier <- exampleCovs
#' exampleCovs_barrier$coast_barrier <- prepBarrier(coast_barrier)
#' exampleCovs_barrier$d2coast <- exampleCovs_barrier$coast_barrier / 100
#'
#' # add a beta coefficient for d2coast to the parameter list
#' par_barrier <- par
#' par_barrier$beta <- c(par_barrier$beta, -0.2)
#'
#' # simulate the data
#' set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
#' simDat_barrier <- simLangevin(par = par_barrier,
#'                               nbAnimals = 3,
#'                               spatialCovs = exampleCovs_barrier,
#'                               barrier = "coast_barrier",
#'                               measurementError = list(smaj.sd = 1.5,
#'                                                       smin.sd = 0.75,
#'                                                       eor.lim = c(0,180)))
#'
#' plot(simDat_barrier,beta=par_barrier$beta,
#'      spatialCovs=exampleCovs_barrier,
#'      maskRast = coast_barrier)
#' }
#' @references
#' Dupont F, McClintock BT, Fischer J-O, Marcoux M, Hussey N, Auger-Methe M. 2025. Inferring resource selection and utilization distributions from irregular and error-prone animal tracking data using the habitat-driven Langevin diffusion.
#'
#' McClintock BT, London JM, Camneron MF, Boveng PL. 2015. Modelling animal movement using the Argos satellite telemetry location error ellipse. Methods Ecol Evol 6:266–277. doi: 10.1111/2041-210X.12286.
#'
#' Michelot T, Gloaguen P, Blackwell PG, Etienne M-P. 2019. The Langevin diffusion as a continuous-time model of animal movement and habitat selection. Methods Ecol Evol 10:1894–1907. doi: 10.1111/2041-210X.13275.
#'
#' Michelot T, Hanks E. 2025. Multiscale modelling of animal movement with persistent dynamics. arXiv doi: 10.48550/arXiv.2406.15195.
#' @useDynLib langevinSSM
#' @export
simLangevin <- function(model, ...) {
  UseMethod("simLangevin")
}


#' @rdname simLangevin
#' @export
simLangevin.default <- function(model = c("underdamped", "overdamped"),
                                spatialCovs,
                                barrier = NULL,
                                par,
                                lambda = NULL,
                                timeStep = 0.01,
                                nbAnimals = 1,
                                obsPerAnimal = 500,
                                initialPosition,
                                measurementError = NULL,
                                subSample = NULL,
                                ...) {

  model <- match.arg(model)
  .validate_lambda(lambda)
  orig_spatialCovs <- spatialCovs

  if (!is.null(barrier)) {
    if (!(barrier %in% names(spatialCovs))) stop(sprintf("Barrier raster '%s' not found in spatialCovs.", barrier))

    barrier_sdf <- spatialCovs[[barrier]]
    if(!isTRUE(attr(barrier_sdf,"barLangevin"))) stop("barrier is not a 'barLangevin' object created by prepBarrier.")
    spatialCovs[[barrier]] <- NULL # Decouple physics constraint from habitat selection
    barrier_dist_mat <- terra::as.matrix(barrier_sdf, wide = TRUE)
  } else {
    barrier_dist_mat <- matrix(0, 1, 1)
  }

  if (length(spatialCovs) == 0) {
    stop("At least one habitat covariate must remain in spatialCovs after isolating the barrier constraint.")
  }

  raster_data <- prepareRaster(spatialCovs)

  if(!is.finite(nbAnimals) || nbAnimals < 1) stop("nbAnimals should be at least 1.")
  if(!is.finite(obsPerAnimal) || obsPerAnimal < 1) stop("obsPerAnimal should be at least 1.")
  if(!is.finite(timeStep) || timeStep <= 0) stop("timeStep should be greater than zero.")

  par <- checkPar(par, model, spatialCovs = spatialCovs)$par
  gamma <- if(!is.null(par$log_gamma)) exp(par$log_gamma) else 1
  sigma <- exp(par$log_sigma)
  beta <- par$beta

  if(!is.null(subSample)){
    if(!is.list(subSample) || !all(names(subSample) %in% c("samplingRate","propMissing"))) stop("subSample must be a list with elements 'samplingRate' and/or 'propMissing'")
    if(!is.null(subSample$samplingRate)){
      if(length(subSample$samplingRate)!=1 || (!is.numeric(subSample$samplingRate) | subSample$samplingRate < 1)) stop("subSample$samplingRate must be a numeric of length 1 that is > 1")
    } else subSample$samplingRate <- 1
    if(!is.null(subSample$propMissing)){
      if(length(subSample$propMissing)!=1 || (!is.numeric(subSample$propMissing) | subSample$propMissing < 0 | subSample$propMissing >= 1)) stop("subSample$propMissing must be a numeric of length 1 that is >= 0 and < 1")
    } else subSample$propMissing <- 0
  }

  dt_vec <- rep(timeStep, obsPerAnimal - 1)

  # --- Temporal Bounding Check for dynamic covariates ---
  sim_min_time <- 0
  sim_max_time <- sum(dt_vec)

  for (j in seq_along(spatialCovs)) {
    if (terra::nlyr(spatialCovs[[j]]) > 1) {
      t_vals <- terra::time(spatialCovs[[j]])

      if (inherits(t_vals, c("POSIXt", "Date"))) {
        stop("When simulating from scratch with dynamic covariates, the 'terra::time' values of 'spatialCovs$", names(spatialCovs)[j], "' must be numeric (not POSIXt or Date) to align with the numeric simulation time (which starts at 0).")
      }

      if (is.numeric(t_vals)) {
        if (sim_min_time < min(t_vals, na.rm = TRUE) || sim_max_time > max(t_vals, na.rm = TRUE)) {
          stop("The simulated tracking times (0 to ", sim_max_time, ") fall outside the temporal boundaries of 'spatialCovs$", names(spatialCovs)[j], "'. Ensure the simulation time span is entirely within the raster's numeric time range.")
        }
      }
    }
  }

  # --- 3. Lambda Auto-Scaling ---
  if (!is.null(barrier) && is.null(lambda)) {
    max_dt <- max(dt_vec, na.rm = TRUE)
    lambda <- .calc_lambda_limit(model, sigma, gamma, max_dt)
    message("   Auto-scaling barrier lambda based on true simulation parameters: ", signif(lambda, 4))
  }

  barrier_pen <- if (!is.null(barrier)) lambda else 0

  # --- 4. Initial Position ---
  init_pos_sim <- if (missing(initialPosition)) {
    getInitialPosition(nbAnimals = nbAnimals, spatialCovs = orig_spatialCovs, barrier = barrier,
                       beta = beta, lambda = lambda)
  } else {
    getInitialPosition(nbAnimals = nbAnimals, initialPosition = initialPosition,
                       spatialCovs = spatialCovs, beta = beta, lambda = lambda)
  }

  out <- simulate_langevin_cpp(
    model = ifelse(model=="underdamped", 1, 0),
    nbAnimals = nbAnimals,
    obsPerAnimal = length(dt_vec),
    dt_vec = dt_vec,
    gamma = gamma,
    sigma = sigma,
    beta = beta,
    raster_data = raster_data,
    initialPosition = init_pos_sim,
    barrier_dist = barrier_dist_mat,
    barrier_penalty = barrier_pen
  )

  class(out) <- c("dataLangevin", class(out))

  if(!is.null(measurementError)) {
    out <- addMeasurementError(out, par, measurementError = measurementError)
  } else {
    out <- out %>% dplyr::mutate(x = mu.x, y = mu.y, smaj = NA, smin = NA, eor = NA, x.err = NA, y.err = NA)
  }

  out$id <- as.factor(out$id)

  if(model == "underdamped") {
    out <- out %>% dplyr::select(id, date, dt, dplyr::any_of("lc"), x, y, smaj, smin, eor, x.err, y.err, mu.x, mu.y, vel.x, vel.y)
  } else if(model == "overdamped") {
    out <- out %>% dplyr::select(id, date, dt, dplyr::any_of("lc"), x, y, smaj, smin, eor, x.err, y.err, mu.x, mu.y)
  }

  out <- class_dataLangevin(out)
  class(out) <- unique(c("simLangevin", class(out)))

  if (!is.null(barrier)) {
    attr(out, "lambda") <- lambda
    attr(out, "barrier") <- barrier
  }

  if(!is.null(subSample)){
    orig_mean_dt <- mean(out$dt[out$dt > 0], na.rm = TRUE)
    out <- subSampleData(out, samplingRate = subSample$samplingRate, propMissing = subSample$propMissing)
    new_mean_dt <- mean(out$dt[out$dt > 0], na.rm = TRUE)

    effective_lambda <- .rescale_lambda(lambda, model, sigma, gamma, orig_mean_dt, new_mean_dt)
    attr(out, "lambda") <- effective_lambda

    if(!is.null(barrier) & isTRUE(subSample$samplingRate>1)) {
      message(sprintf("   Subsampling degraded temporal resolution. Effective lambda: %.4f", effective_lambda))
    }
  }

  return(out)
}


#' @rdname simLangevin
#' @export
simLangevin.fitLangevin <- function(model,
                                    data,
                                    spatialCovs,
                                    timeStep = 0.01,
                                    jointPrecision = FALSE,
                                    conditional = FALSE,
                                    ...) {

  if (is.null(data) || !inherits(data,"dataLangevin")) stop("data must be a dataLangevin object.")
  if (nrow(data) == 0) stop("data contains no observations.")

  time.unit <- attr(data, "time.unit")

  # --- POSIXt Type Enforcement & Parsing ---
  if (inherits(data$date, c("POSIXt", "Date"))) {
    if (!is.character(timeStep)) {
      stop("When the 'date' column is POSIXt or Date, 'timeStep' must be a character string specifying the time interval (e.g., '1 sec', '30 mins', '1 hour', '6 hours').")
    }
    if (is.null(time.unit)) stop("When the 'date' column is POSIXt or Date, the 'data' object must have a 'time.unit' attribute specifying the time unit (e.g., 'secs', 'mins', 'hours').")

    t0 <- as.POSIXct("1970-01-01", tz = "UTC")
    t1 <- tryCatch(seq(t0, by = timeStep, length.out = 2)[2], error = function(e) NA)

    if (is.na(t1)) {
      stop("Invalid 'timeStep' string provided. Use valid base R formats like '1 sec', '30 mins', or '1 hour'.")
    }

    timeStep <- as.numeric(difftime(t1, t0, units = time.unit))

    if (timeStep <= 0) stop("Invalid 'timeStep' string provided. Must result in a positive duration.")
  } else {
    if (!is.numeric(timeStep)) {
      stop("When the 'date' column is numeric, 'timeStep' must also be numeric.")
    }
  }

  if (!conditional && (!is.finite(timeStep) || timeStep <= 0)) {
    stop("timeStep must be a valid positive value.")
  }

  fit <- model
  verify_signatures(fit, data = data, spatialCovs = spatialCovs)

  cond <- fit$conditions
  coord <- cond$coord
  scaleFactor <- cond$scaleFactor
  lambda <- cond$lambda

  if(!is.null(lambda)){

    orig_mean_dt <- mean(data$dt[data$dt > 0], na.rm = TRUE)
    new_mean_dt <- timeStep

    sigma <- exp(fit$tmb_setup$parList$log_sigma)
    gamma <- if (cond$model == "underdamped") exp(fit$tmb_setup$parList$log_gamma) else 1

    lambda <- .rescale_lambda(lambda, cond$model, sigma, gamma, orig_mean_dt, new_mean_dt)
  }

  if (!is.null(cond$barrier)) {
    barrier_sdf <- spatialCovs[[cond$barrier]] / scaleFactor
    spatialCovs[[cond$barrier]] <- NULL
  } else {
    barrier_sdf <- NULL
  }

  dat <- build_tmb_data(data, spatialCovs, cond$model, coord, scaleFactor,
                        cond$smoothGradient, cond$npoints, cond$curweight, cond$zetaScale,
                        barrier_sdf = barrier_sdf, lambda = lambda)

  dat <- c(dat, fit$tmb_setup$priors)

  obj2 <- try({
    TMB::MakeADFun(
      data = c(model = "langevinSSM", dat),
      parameters = fit$tmb_setup$parList,
      map = fit$tmb_setup$map,
      random = fit$tmb_setup$random,
      DLL = "langevinSSM_TMBExports",
      silent = TRUE
    )
  }, silent = TRUE)

  if (inherits(obj2, "try-error")) {
    stop("Failed to reconstruct TMB objective function. Error: ", attr(obj2, "condition")$message)
  }

  obj2$fn(fit$par)

  if ((jointPrecision | conditional) && !requireNamespace("Matrix", quietly = TRUE)) stop("The 'Matrix' package is required for Monte Carlo simulation. Please install it.")

  calc_jp <- jointPrecision && is.null(fit$covariance$random$jointPrecision)
  if(calc_jp) message("   Calculating full covariance matrix using TMB::sdreport(obj, getJointPrecision = TRUE)")
  sdr <- TMB::sdreport(obj2, getJointPrecision = calc_jp)

  full_par <- obj2$env$last.par
  fixed_est <- sdr$par.fixed
  random_est <- sdr$par.random

  samp_fixed <- fixed_est
  samp_random <- random_est

  if (jointPrecision) {
    if (conditional) {
      message("   Imputing tracks tied to data using the full joint covariance matrix...")
    } else {
      message("   Simulating tracks forward using the full joint covariance matrix...")
    }

    Q <- if (!is.null(fit$covariance$random$jointPrecision)) fit$covariance$random$jointPrecision else sdr$jointPrecision
    L <- Matrix::Cholesky(Q, super = TRUE)
    z <- stats::rnorm(ncol(Q))
    step <- as.numeric(Matrix::solve(L, Matrix::solve(L, z, system = "Lt"), system = "Pt"))
    samp_fixed <- fixed_est + step[1:length(fixed_est)]
    samp_random <- random_est + step[(length(fixed_est) + 1):length(step)]

  } else if (conditional) {
    message("   Imputing tracks tied to data using the random effects covariance matrix...")

    Huu <- obj2$env$spHess(random = TRUE)
    L_u <- Matrix::Cholesky(Huu, super = TRUE)
    z_u <- stats::rnorm(ncol(Huu))
    step_random <- as.numeric(Matrix::solve(L_u, Matrix::solve(L_u, z_u, system = "Lt"), system = "Pt"))
    samp_random <- random_est + step_random

  } else {
    message("   Simulating tracks forward using movement parameters fixed at point estimates...")
  }

  full_par[obj2$env$lfixed()] <- samp_fixed
  full_par[obj2$env$lrandom()] <- samp_random

  par_drawn <- obj2$env$parList(par = full_par)

  mu_working <- t(matrix(par_drawn$mu, nrow = 2))
  vel_working <- if (fit$conditions$model == "underdamped") t(matrix(par_drawn$vel, nrow = 2)) else NULL

  # Construct natural parameters using the scaled sigma
  nat_par <- list(
    beta = par_drawn$beta,
    sigma = exp(par_drawn$log_sigma) * scaleFactor, # original scale sigma
    gamma = if(!is.null(par_drawn$log_gamma)) exp(par_drawn$log_gamma) else 1,
    psi = if(!is.null(par_drawn$l_psi)) exp(par_drawn$l_psi) else 1,
    tau = if(!is.null(par_drawn$l_tau)) exp(par_drawn$l_tau) else c(1, 1),
    rho_o = if(!is.null(par_drawn$l_rho_o)) 2/(1 + exp(-par_drawn$l_rho_o)) - 1 else 0
  )

  sigma_working <- exp(par_drawn$log_sigma)

  if (conditional) {
    out <- data
    out$mu.x <- mu_working[, 1] * scaleFactor
    out$mu.y <- mu_working[, 2] * scaleFactor
    if (fit$conditions$model == "underdamped") {
      out$vel.x <- vel_working[, 1] * scaleFactor
      out$vel.y <- vel_working[, 2] * scaleFactor
    }

    out <- addMeasurementError(out, nat_par)
  } else {
    out_list <- vector("list", length(unique(data$id)))
    names(out_list) <- unique(data$id)

    min_dt <- min(data$dt[data$dt > 0], na.rm = TRUE)
    prec <- max(10, ceiling(abs(log10(min_dt))) + 3)

    for (i in unique(data$id)) {
      ind_data <- data[data$id == i, ]

      if (inherits(ind_data$date, c("POSIXt", "Date"))) {
        obs_times <- as.numeric(difftime(ind_data$date, ind_data$date[1], units = time.unit))
      } else {
        obs_times <- as.numeric(ind_data$date) - as.numeric(ind_data$date[1])
      }

      acc_grid <- seq(0, max(obs_times), by = timeStep)

      full_grid <- sort(unique(round(c(obs_times, acc_grid), prec)))
      full_dt_vec <- diff(full_grid)
      obs_idx <- match(round(obs_times, prec), round(full_grid, prec))

      first_idx <- which(data$id == i)[1]

      # Feed the working scale initial positions to the solver
      init_pos <- matrix(c(mu_working[first_idx, 1], mu_working[first_idx, 2]), nrow = 1, ncol = 2)

      raster_names <- c("raster_vals","raster_coords","raster_resolution","raster_extent","n_covs","all_z_values","n_zvals_cov","cov_offset")

      sim_full <- simulate_langevin_cpp(
        model = ifelse(fit$conditions$model=="underdamped", 1, 0),
        nbAnimals = 1,
        obsPerAnimal = length(full_dt_vec),
        dt_vec = full_dt_vec,
        gamma = nat_par$gamma,
        sigma = sigma_working, # Provide working scale sigma
        beta = nat_par$beta,
        raster_data = dat[raster_names],
        initialPosition = init_pos,
        barrier_dist = dat$barrier_dist,
        barrier_penalty = dat$barrier_penalty
      )

      sim_ind <- sim_full[obs_idx, ]

      # rescale the outputs back to the natural scale
      sim_ind$mu.x <- sim_ind$mu.x * scaleFactor
      sim_ind$mu.y <- sim_ind$mu.y * scaleFactor
      if (fit$conditions$model == "underdamped") {
        sim_ind$vel.x <- sim_ind$vel.x * scaleFactor
        sim_ind$vel.y <- sim_ind$vel.y * scaleFactor
      }

      sim_ind$id <- i
      sim_ind$date <- ind_data$date
      sim_ind$dt <- ind_data$dt

      if ("smaj" %in% names(ind_data)) {
        sim_ind$smaj <- ind_data$smaj; sim_ind$smin <- ind_data$smin; sim_ind$eor <- ind_data$eor
      }
      if ("x.err" %in% names(ind_data)) {
        sim_ind$x.err <- ind_data$x.err; sim_ind$y.err <- ind_data$y.err
      }
      if ("lc" %in% names(ind_data)) {
        sim_ind$lc <- ind_data$lc
      }

      out_list[[as.character(i)]] <- sim_ind
    }
    out <- do.call(rbind, out_list)

    class(out) <- c("dataLangevin", class(out))

    out <- addMeasurementError(out, nat_par)
  }

  out$id <- as.factor(out$id)

  if (fit$conditions$model == "underdamped") {
    out <- out %>% dplyr::select(id, date, dt, dplyr::any_of("lc"), x, y, smaj, smin, eor, x.err, y.err, mu.x, mu.y, vel.x, vel.y)
  } else {
    out <- out %>% dplyr::select(id, date, dt, dplyr::any_of("lc"), x, y, smaj, smin, eor, x.err, y.err, mu.x, mu.y)
  }

  attr(out,"time.unit") <- time.unit
  out <- class_dataLangevin(out)
  class(out) <- unique(c("simLangevin", class(out)))

  if (!is.null(cond$barrier)) {
    attr(out, "lambda") <- cond$lambda
    attr(out, "barrier") <- cond$barrier
    attr(out, "scaleFactor") <- scaleFactor
  }

  return(out)
}

# --- Lambda Rescaling Helpers ---
.calc_lambda_limit <- function(model, sigma, gamma, dt) {
  if (model == "overdamped") {
    return(2 / (sigma^2 * dt))
  } else {
    num <- (gamma^2) * (1 - exp(-gamma * dt))
    den <- (sigma^2) * (1 - exp(-gamma * dt) - (gamma * dt * exp(-gamma * dt)))
    return(num / den)
  }
}

.rescale_lambda <- function(lambda, model, sigma, gamma, orig_mean_dt, new_mean_dt) {
  L_orig <- .calc_lambda_limit(model, sigma, gamma, orig_mean_dt)
  L_new <- .calc_lambda_limit(model, sigma, gamma, new_mean_dt)
  return(lambda * (L_new / L_orig))
}
