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
#' @param timeStep Time step to use for the simulation. Determines the resolution of the discrete-time approximation of the continuous-time process. The smaller the \code{timeStep}, the more accurate the approximation. Ignored if \code{model} is a \code{fitLangevin} object and \code{conditional=TRUE}. Default: 0.01.
#'
#' @param par List of parameters. For the "underdamped" model, this must include \code{beta} (a vector of length equal to the number of covariates), \code{sigma} (speed parameter), and \code{gamma} (friction parameter). For the "overdamped" model, this must include \code{beta} and \code{sigma}. If \code{measurementError} is specified, optional observation process parameters include \code{psi}, \code{tau}, and \code{rho_o}. See \code{\link{fitLangevin}}.
#' @param nbAnimals Number of animals to simulate. Default: 1.
#' @param obsPerAnimal Number of observations to simulate per animal. Default: 500.
#' @param initialPosition Initial position(s) for the simulation. A 2-vector, or a list of length \code{nbAnimals} of 2-vectors. If missing, initial positions are randomly generated based on the utilization distribution.
#' @param measurementError A list or data frame of parameters used to simulate observation error. See \code{\link{addMeasurementError}}. Default: \code{NULL}.
#' @param subSample List of specifications for subsampling data from the continuous-time process model, which can include \code{samplingRate} and \code{propMissing}. See \code{\link{subSampleData}}. Default: \code{NULL} (no subsampling or missing data).
#' @param barrier Optional character string specifying the name of the barrier mask within \code{spatialCovs}. This must be a binary raster where 1 indicates allowed movement areas and 0 indicates restricted areas. See Details.
#' @param lambda Numeric. The penalty weight for the barrier constraint. Larger values create a steeper "wall" preventing locations from crossing into restricted areas. Default: \code{NULL}. Because the true generative parameters are known during simulation, leaving this as \code{NULL} allows the function to perfectly auto-calculate the optimal theoretical stability limit based on the \code{model} type, true speed parameter (\code{sigma}), and \code{timeStep}. See Details.
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
#'                       spatialCovs = exampleCovs)
#'
#' # conditional on observed tracks
#' simFit_cond <- simLangevin(fit,
#'                            data = exampleDat,
#'                            spatialCovs = exampleCovs,
#'                            conditional = TRUE)
#'
#' # simulating with a barrier
#' # create a dummy barrier mask (left half restricted = 0, right half allowed = 1)
#' coast_barrier <- exampleCovs[[1]]
#' terra::values(coast_barrier) <- ifelse(terra::crds(coast_barrier)[, "x"] >=
#'                                        mean(terra::crds(coast_barrier)[, "x"]), 1, 0)
#' names(coast_barrier) <- "coast_barrier"
#'
#' # add the mask to the spatial covariates list
#' exampleCovs_barrier <- exampleCovs
#' exampleCovs_barrier$coast_barrier <- coast_barrier
#'
#' # add a beta coefficient for the barrier to the parameter list
#' par_barrier <- par
#' par_barrier$beta <- c(par_barrier$beta, -0.2)
#'
#' # simulate the data
#' set.seed(123,kind="Mersenne-Twister",normal.kind="Inversion")
#' simDat_barrier <- simLangevin(par = par_barrier,
#'                               spatialCovs = exampleCovs_barrier,
#'                               barrier = "coast_barrier",
#'                               measurementError = list(smaj.sd = 1.5,
#'                                                       smin.sd = 0.75,
#'                                                       eor.lim = c(0,180)))
#'
#' plot(simDat_barrier,beta=par_barrier$beta,spatialCovs=exampleCovs_barrier)
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
                                par,
                                spatialCovs,
                                nbAnimals = 1,
                                obsPerAnimal = 500,
                                timeStep = 0.01,
                                initialPosition,
                                measurementError = NULL,
                                subSample = NULL,
                                barrier = NULL,
                                lambda = NULL,
                                ...) {

  model <- match.arg(model)

  .validate_barrier(barrier, spatialCovs)
  .validate_lambda(lambda)

  orig_spatialCovs <- spatialCovs
  if (!is.null(barrier)) {
    barrier_sdf <- .get_barrier_sdf(barrier, spatialCovs)
    spatialCovs[[barrier]] <- barrier_sdf
    barrier_dist_mat <- terra::as.matrix(barrier_sdf, wide = TRUE)
  } else {
    barrier_dist_mat <- matrix(0, 1, 1)
  }

  # spatialCovs now contains the continuous SDF instead of the flat 0/1 mask
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
    # Because we know the TRUE sigma and dt, we can perfectly calculate the limit
    max_dt <- max(dt_vec, na.rm = TRUE)
    if (model == "overdamped") {
      lambda <- 2 / (sigma^2 * max_dt)
    } else {
      lambda <- 1 / (sigma^2 * max_dt^2)
    }
  }

  barrier_pen <- if (!is.null(barrier)) lambda else 0

  # --- 4. Initial Position ---
  init_pos_sim <- if (missing(initialPosition)) {
    getInitialPosition(nbAnimals = nbAnimals, spatialCovs = orig_spatialCovs,
                       beta = beta, barrier = barrier, lambda = lambda)
  } else {
    getInitialPosition(nbAnimals = nbAnimals, initialPosition = initialPosition,
                       spatialCovs = orig_spatialCovs, beta = beta)
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
    out <- subSampleData(out, samplingRate = subSample$samplingRate, propMissing = subSample$propMissing)
    max_dt <- max(out$dt, na.rm = TRUE)
    if (model == "overdamped") {
      lambda <- 2 / (sigma^2 * max_dt)
    } else {
      lambda <- 1 / (sigma^2 * max_dt^2)
    }
    attr(out, "lambda") <- lambda
  }

  if(!is.null(barrier)) message("   Auto-scaling barrier lambda based on true simulation parameters: ", signif(lambda, 4))

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

  if (!conditional && (!is.finite(timeStep) || timeStep <= 0)) {
    stop("timeStep should be greater than zero.")
  }

  fit <- model
  verify_signatures(fit, data = data, spatialCovs = spatialCovs)

  cond <- fit$conditions
  time.unit <- attr(data, "time.unit")
  coord <- cond$coord
  scaleFactor <- cond$scaleFactor

  dat <- build_tmb_data(data, spatialCovs, cond$model, coord, scaleFactor,
                        cond$smoothGradient, cond$npoints, cond$curweight, cond$zetaScale,
                        barrier_names = cond$barrier, lambda = cond$lambda)

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
  mu_mat <- t(matrix(par_drawn$mu, nrow = 2)) * scaleFactor
  vel_mat <- if (fit$conditions$model == "underdamped") t(matrix(par_drawn$vel, nrow = 2)) * scaleFactor else NULL

  nat_par <- list(
    beta = par_drawn$beta,
    sigma = exp(par_drawn$log_sigma) * scaleFactor,
    gamma = if(!is.null(par_drawn$log_gamma)) exp(par_drawn$log_gamma) else 1,
    psi = if(!is.null(par_drawn$l_psi)) exp(par_drawn$l_psi) else 1,
    tau = if(!is.null(par_drawn$l_tau)) exp(par_drawn$l_tau) else c(1, 1),
    rho_o = if(!is.null(par_drawn$l_rho_o)) 2/(1 + exp(-par_drawn$l_rho_o)) - 1 else 0
  )

  if (conditional) {
    out <- data
    out$mu.x <- mu_mat[, 1]
    out$mu.y <- mu_mat[, 2]
    if (fit$conditions$model == "underdamped") {
      out$vel.x <- vel_mat[, 1]
      out$vel.y <- vel_mat[, 2]
    }
    out <- addMeasurementError(out, nat_par)
  } else {
    out_list <- vector("list", length(unique(data$id)))
    names(out_list) <- unique(data$id)

    min_dt <- min(data$dt[data$dt > 0], na.rm = TRUE)
    prec <- max(10, ceiling(abs(log10(min_dt))) + 3)

    for (i in unique(data$id)) {
      ind_data <- data[data$id == i, ]
      obs_times <- as.numeric(ind_data$date) - as.numeric(ind_data$date[1])
      acc_grid <- seq(0, max(obs_times), by = timeStep)

      full_grid <- sort(unique(round(c(obs_times, acc_grid), prec)))
      full_dt_vec <- diff(full_grid)
      obs_idx <- match(round(obs_times, prec), round(full_grid, prec))

      first_idx <- which(data$id == i)[1]
      init_pos <- matrix(c(mu_mat[first_idx, 1], mu_mat[first_idx, 2]), nrow = 1, ncol = 2)

      raster_names <- c("raster_vals","raster_coords","raster_resolution","raster_extent","n_covs","all_z_values","n_zvals_cov","cov_offset")

      sim_full <- simulate_langevin_cpp(
        model = ifelse(fit$conditions$model=="underdamped", 1, 0),
        nbAnimals = 1,
        obsPerAnimal = length(full_dt_vec),
        dt_vec = full_dt_vec,
        gamma = nat_par$gamma,
        sigma = nat_par$sigma / scaleFactor,
        beta = nat_par$beta,
        raster_data = dat[raster_names],
        initialPosition = init_pos,
        barrier_dist = dat$barrier_dist,
        barrier_penalty = dat$barrier_penalty
      )

      sim_ind <- sim_full[obs_idx, ]
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
  }

  return(out)
}
