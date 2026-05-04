# --- Helper: Build TMB Data List ---
build_tmb_data <- function(data, spatialCovs, model, coord, scaleFactor, smoothGradient, npoints, curweight, zetaScale, barrier_names = NULL, lambda = NULL) {

  if (!is.null(barrier_names)) {
    barrier_sdf <- spatialCovs[[barrier_names]]
    barrier_dist_mat <- terra::as.matrix(barrier_sdf, wide = TRUE)
    barrier_pen <- lambda
  } else {
    barrier_dist_mat <- matrix(0, 1, 1)
    barrier_pen <- 0
  }

  time.unit <- attr(data, "time.unit")
  raster_data <- prepareRaster(spatialCovs, scaleFactor = scaleFactor, time.unit = time.unit, data = data, coord = coord)

  if (inherits(data$date, "POSIXt") || inherits(data$date, "Date")) {
    track_times <- as.numeric(difftime(data$date, as.POSIXct("1970-01-01 00:00:00", tz = "UTC"), units = time.unit))
  } else {
    track_times <- as.numeric(data$date)
  }

  dat <- list(
    process_model = ifelse(model == "underdamped", 1, 0),
    Y = t(data[, coord]) / scaleFactor,
    times = track_times,
    dt = data$dt
  )

  dat$skip_step <- as.integer(dat$dt < 1.e-6)
  dat$smaj <- data$smaj / scaleFactor
  dat$smin <- data$smin / scaleFactor
  dat$eor <- data$eor
  dat$K <- as.matrix(data[, c("x.err", "y.err")] / scaleFactor)
  dat$isd <- as.numeric(!is.na(dat$Y[1,]) & ((!is.na(dat$K[,1]) & !is.na(dat$K[,2])) | (!is.na(dat$smaj) & !is.na(dat$smin) & !is.na(dat$eor))))
  dat$obs_mod <- rep(NA, ncol(dat$Y))
  dat$obs_mod[dat$isd == 1 & (!is.na(dat$K[,1]) & !is.na(dat$K[,2]))] <- 0
  dat$obs_mod[dat$isd == 1 & (!is.na(dat$smaj) & !is.na(dat$smin) & !is.na(dat$eor))] <- 1
  dat$ID <- data$id
  dat$nbObs <- rep(1, ncol(dat$Y))
  dat$scale_factor <- scaleFactor
  dat$smoothGradient <- ifelse(smoothGradient, 1, 0)
  dat$weights <- c(curweight, rep((1 - curweight) / npoints, npoints))
  dat$zetaScale <- zetaScale

  dat$barrier_dist <- barrier_dist_mat
  dat$barrier_penalty <- barrier_pen
  dat <- c(dat, raster_data)

  return(dat)
}

# --- Helper: Extract and Format TMB Estimates ---
extract_tmb_estimates <- function(fit, obj, sdreport_out, re, map, data, scaleFactor, spatialCovs, getJointPrecision) {
  estimates <- list()
  covariance <- list()

  if (inherits(sdreport_out, "try-error")) {
    warning("TMB::sdreport failed to calculate standard errors. Trajectory and point estimates were recovered from the report, but SEs are unavailable.")

    # 1. working scale
    working_names <- names(fit$par)
    name_counts <- table(working_names)
    is_dup <- working_names %in% names(name_counts[name_counts > 1])
    if (any(is_dup)) {
      suffix <- ave(working_names, working_names, FUN = seq_along)
      working_names[is_dup] <- paste0(working_names[is_dup], "_", suffix[is_dup])
    }
    estimates$working <- data.frame("Estimate" = as.numeric(fit$par), "Std. Error" = NA_real_, check.names = FALSE)
    rownames(estimates$working) <- working_names

    # 2. natural scale
    rep_vals <- obj$report()
    nat_names <- c("beta", "sigma", "gamma", "rho_o", "tau", "psi")
    existing_names <- nat_names[nat_names %in% names(rep_vals)]
    nat_list <- lapply(existing_names, function(nm) rep_vals[[nm]])
    names(nat_list) <- existing_names
    nat_est <- unlist(nat_list, use.names=FALSE)

    nat_rn <- unlist(lapply(names(nat_list), function(x) {
      if (length((nat_list[[x]])) == 1) return(x)
      else return(paste0(x, "_", 1:length(nat_list[[x]])))
    }))

    estimates$natural <- data.frame("Estimate" = as.numeric(nat_est), "Std. Error" = NA_real_, check.names = FALSE)
    rownames(estimates$natural) <- nat_rn

    # 3. random effects (mu/vel)
    if(length(re)) {
      estimates$random <- list()
      for(i in seq_along(re)) {
        node <- re[i]
        estimates$random[[node]] <- list(
          est = data.frame(id = data$id, date = data$date, t(rep_vals[[node]]) * scaleFactor),
          se = NULL
        )
        colnames(estimates$random[[node]]$est) <- c("id", "date", paste0(node, ".", c("x", "y")))
      }
    }
  } else {

    estimates$natural <- as.data.frame(summary(sdreport_out, "report"))
    estimates$working <- as.data.frame(summary(sdreport_out, "fixed"))

    covariance$natural <- sdreport_out$cov
    covariance$working <- sdreport_out$cov.fixed

    for (est_type in c("natural", "working")) {
      rn <- rownames(estimates[[est_type]])
      name_counts <- table(rn)
      is_dup <- rn %in% names(name_counts[name_counts > 1])
      if (any(is_dup)) {
        suffix <- ave(rn, rn, FUN = seq_along)
        rn[is_dup] <- paste0(rn[is_dup], "_", suffix[is_dup])
      }
      rownames(estimates[[est_type]]) <- rn
    }

    # 3. random effects
    if(length(re)) {
      ran_est <- tryCatch(summary(sdreport_out, "random"), error = function(e) NULL)

      obj$fn(fit$par)
      rep_vals <- obj$report()
      estimates$random <- list()

      for(i in seq_along(re)) {
        node <- re[i]
        est_mat <- t(matrix(rep_vals[[node]], nrow = 2)) * scaleFactor
        se_full <- rep(0.0, length(rep_vals[[node]]))

        if (!is.null(ran_est) && nrow(ran_est) > 0 && node %in% rownames(ran_est)) {
          node_ran_est <- ran_est[rownames(ran_est) == node, , drop = FALSE]

          if (!is.null(map[[node]])) {
            map_int <- as.integer(map[[node]])
            valid_idx <- !is.na(map_int)
            if (any(valid_idx)) {
              se_full[valid_idx] <- node_ran_est[map_int[valid_idx], "Std. Error"]
            }
          } else {
            se_full <- node_ran_est[, "Std. Error"]
          }
        }

        se_mat <- t(matrix(se_full, nrow = 2)) * scaleFactor
        estimates$random[[node]] <- list()
        estimates$random[[node]]$est <- data.frame(id = data$id, date = data$date, est_mat)
        estimates$random[[node]]$se  <- data.frame(id = data$id, date = data$date, se_mat)
        colnames(estimates$random[[node]]$est) <- c("id", "date", paste0(node, ".", c("x", "y")))
        colnames(estimates$random[[node]]$se)  <- c("id", "date", paste0(node, ".", c("x", "y")))
      }
      if(getJointPrecision) covariance$random$jointPrecision <- sdreport_out$jointPrecision
    }
  }

  for (est_type in c("natural", "working")) {
    if (!is.null(estimates[[est_type]])) {
      rn <- rownames(estimates[[est_type]])

      beta_idx <- which(rn == "beta" | grepl("^beta(_|[.])", rn))
      if (length(beta_idx) > 0 && !is.null(names(spatialCovs))) {
        # Only map up to the number of covariates provided
        n_covs <- length(spatialCovs)
        n_to_map <- min(length(beta_idx), n_covs)
        for(i in 1:n_to_map) {
          rn[beta_idx[i]] <- paste0("beta_", names(spatialCovs)[i])
        }
      }

      # 2. Map Tau Names: tau_1 and tau_2
      tau_idx <- which(rn == "tau" | grepl("^tau(_|[.])", rn))
      if (length(tau_idx) == 2) {
        rn[tau_idx] <- c("tau_1", "tau_2")
      }

      rownames(estimates[[est_type]]) <- rn
    }
  }

  return(list(estimates = estimates, covariance = covariance))
}

#' Fit the habitat-driven Langevin diffusion
#'
#' Fit the underdamped or overdamped Langevin model to the location data provided, via numerical optimization of the log-likelihood
#' function using \code{TMB}. Location data can be temporally irregular and/or subject to measurement error.
#'
#' @param data \code{dataLangevin} object. See \code{\link{formatData}}.
#' @param model Character string indicating which Langevin diffusion model to fit (``underdamped'', or ``overdamped''). Default: ``underdamped''.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates.
#' @param par List containing the initial values for the parameters. These can include state process parameters for the habitat selection coefficients (``beta''), the ``speed'' parameter (``sigma''), and, for the underdamped model, the friction coefficient (``gamma'').
#' Observation process parameters include a scaling factor to account for uncertainty in the Argos error ellipse (``psi''), a scaling factor to account for uncertainty in the x- and y-axis errors for Argos least squares or GPS observations (``tau''), and a correlation term between the x- and y-axis errors for Argos least squares or GPS observations (``rho_o''). All parameter are specified on their natural scale and are converted to working scale internally.
#' Any missing state process parameters are generated using \code{\link{initialValues}}. Any missing observation process parameters are fixed to their default values (``psi'' = 1, ``tau'' = 1, and ``rho_o'' = 0) via \code{map}. See Details.
#' @param prior Optional 2-column data frame containing the mean (column 1) and standard deviation (column 2) for normally distributed priors on the working scale parameters. The row names must match the working scale parameter names (e.g., \code{"beta_cov1"}, \code{"log_sigma"}, \code{"log_gamma"}). Supplying the base name of a vector or matrix parameter (e.g., \code{"beta"}, \code{"mu"}, \code{"vel"}) will apply the prior to all of its elements. To target specific coordinates and time steps for the random effects, append the coordinate (\code{.x} or \code{.y}) and the row index of the observation to the base name (e.g., \code{"mu.x_1"} for the x-coordinate of the 1st observation in \code{data}, or \code{"vel.y_10"} for the y-velocity of the 10th observation in \code{data}). Parameters omitted from this data frame are assigned flat (improper) priors. Default: \code{NULL} (no priors).
#' @param map List defining how to optionally collect and fix parameters. See \code{\link[TMB]{MakeADFun}}.
#' @param coord Character vector identifying the coordinate names for the location data. Default: \code{c("x","y")}.
#' @param lambda Numeric. The penalty weight for the barrier constraint. Default: \code{NULL}. If \code{NULL}, the function will attempt to automatically extract a penalty value if \code{data} is a simulated object created by \code{\link{simLangevin}}. Otherwise, it must be provided. See Details and \code{\link{fitLangevin_barrier}} for a strategy to determine the optimal penalty.
#' @param scaleFactor Internal scaling factor for the coordinates and parameters. In some cases, setting \code{scaleFactor>1} can help with optimization.
#' @param smoothGradient Logical indicating whether or not to smooth the gradients. See Details. Default: \code{FALSE}.
#' @param npoints Number of smoothing points around current cell (4 = diagonal, 8 = queen neighborhood). Ignored unless \code{smoothGradient=TRUE}.
#' @param curweight Smoothing weight of current cell location (\code{0 <= curweight < 1}). Ignored unless \code{smoothGradient=TRUE}.
#' @param zetaScale Scale factor for smooth gradient neighborhood (\code{zetaScale>1} increases and \code{zetaScale<1} decreases the neighborhood). Ignored unless \code{smoothGradient=TRUE}.
#' @param hessian Logical indicating whether or not to calculate the Hessian at the optimum. See \code{\link[TMB]{MakeADFun}}. Default: \code{FALSE}.
#' @param silent Logical indicating whether or not to disable TMB tracing information. See \code{\link[TMB]{MakeADFun}}. Default: \code{FALSE}.
#' @param method Outer optimization method. Default: \code{"BFGS"}.
#' @param initialInner Logical indicating whether or not to first perform an inner optimization for the random effects (``mu'' and/or ``vel'') before optimizing over all parameters. Default: \code{TRUE}.
#' @param inner.control List controlling inner optimization. See \code{\link[TMB]{MakeADFun}}.
#' @param control A list of control parameters for outer optimization. See \code{\link[stats]{nlminb}}.
#' @param polishOptim Logical indicating whether or not to perform an additional ``polishing'' optimization after the initial optimization with \code{\link[stats]{nlminb}} has completed. Default: \code{FALSE}.
#' @param getJointPrecision Logical indicating whether or not to return the joint precision matrix for the random effects. Default: \code{FALSE}.
#'
#' @return \code{fitLangevin} object, i.e., a list of:
#' \item{par}{See \code{\link[stats]{nlminb}}}
#' \item{objective}{See \code{\link[stats]{nlminb}}}
#' \item{convergence}{See \code{\link[stats]{nlminb}}}
#' \item{message}{See \code{\link[stats]{nlminb}}}
#' \item{iterations}{See \code{\link[stats]{nlminb}}}
#' \item{evaluations}{See \code{\link[stats]{nlminb}}}
#' \item{elapsedTime}{Run time of the optimization}
#' \item{estimates}{List containing point estimates and standard errors for the natural scale parameters (``natural''), the working scale parameters (``working''), and the random effects (``random''), where ``random'' is itself a list containing point estimates (``est'') and standard errors (``se'') for the true locations (``mu'') and/or the true velocities (``vel'').}
#' \item{covariance}{List containing the covariance matrices for the natural scale parameters (``natural''), the working scale parameters (``working''), and, if \code{getJointPrecision=TRUE}, the joint covariance matrix for the random effects (``random'')}.
#' \item{conditions}{List containing the optimization settings}
#' \item{signatures}{List containing lightweight fingerprints for \code{data} and \code{spatialCovs} to protect downstream functions}
#' \item{tmb_setup}{Blueprint for reconstructing TMB objective function}
#'
#' @details
#' \strong{Measurement Error Models:}
#' The \code{fitLangevin} function accommodates various types of measurement error depending on the provided measurement error data. It handles measurement error through two distinct internal models with the same underlying structure, which differ in how the observation error covariance matrix is constructed based on the available error information:
#' \deqn{(x,y) \sim N((\mu_x,\mu_y), \Sigma),}
#' where \eqn{(x.y)} are the observed locations, \eqn{(\mu_x,\mu_y)} are the true (latent) locations, and \eqn{\Sigma} is the observation error covariance matrix.
#' \itemize{
#'   \item \strong{Error Ellipse Model (Argos Kalman Filter):} This model is used when error ellipse data (\code{smaj}, \code{smin}, and \code{eor}) are present for an observation. The (optional) observation process parameter \code{psi} acts as a scaling multiplier specifically on the semi-minor axis (\code{smin}) to account for uncertainty in the error ellipse.
#'   For the error ellipse model, \eqn{\Sigma} is derived from the semi-major axis, semi-minor axis, and error ellipse orientation (see McClintock et al. 2015), with \code{psi} scaling the semi-minor axis.
#'   \item \strong{Standard Deviation Model (Argos Least Squares, GPS, Generic Locations):} This model is used when x- and y-axis standard deviations (\code{x.err} and \code{y.err}) are present for an observation. The (optional) observation parameter \code{tau} (a 2-vector) scales the \code{x.err} and \code{y.err}, respectively, to account for uncertainty in the x- and y-axis standard deviations, while the (optional) parameter \code{rho_o} accounts for correlation between the x- and y-axis standard deviations.
#'   For the standard deviation model, \deqn{\Sigma = \begin{bmatrix} \tau_1^2 \sigma_x^2 & \rho_o \tau_1 \tau_2 \sigma_x \sigma_y \\ \rho_o \tau_1 \tau_2 \sigma_x \sigma_y & \tau_2^2 \sigma_y^2 \end{bmatrix},}
#'   with \eqn{\text{tau}=(\tau_1,\tau_2)} scaling the standard deviations \eqn{(\text{x.err}=\sigma_x,\text{y.err}=\sigma_y)} and rho_o \eqn{=\rho_o} accounting for the correlation between the x- and y-axis errors.
#' }
#'
#' \strong{Observation Process Parameters (\code{psi}, \code{tau}, \code{rho_o}):}
#' These parameters control the scaling and correlation of measurement errors and can be specified via \code{par} and \code{map}:
#' \itemize{
#'   \item \strong{1. Fixed at defaults:} If observation process parameters are omitted from \code{par}, the model automatically fixes them to neutral values: \code{psi = 1} (no scaling semi-minor axis), \code{tau = c(1, 1)} (no scaling of \code{x.err} or \code{y.err}), and \code{rho_o = 0} (no correlation between x- and y-axis errors).
#'   \item \strong{2. Estimated from data:} To estimate one or more of these parameters, include them in \code{par} with an initial starting value (e.g., \code{par = list(psi = 1.2)}).
#'   \item \strong{3. Fixed to custom values:} To fix a parameter to a specific value other than the default, include the custom value in \code{par} and explicitly map it to \code{NA} in \code{map} (e.g., \code{par = list(psi = 1.5), map = list(psi = factor(NA))}).
#' }
#'
#' \strong{Gradient smoothing:} Setting \code{smoothGradient=TRUE} applies the smoothing approach of Blackwell and Matthiopoulos (2024) to the gradients of the habitat selection covariates (\code{spatialCovs}), which can help reduce attenuation bias in habitat selection coefficients when observations are obtained at a coarser time resolution than the underlying continuous-time movement process.
#' This approach smooths the gradient at each location by taking a weighted average of the gradients at neighboring locations, where the weights are determined by a Gaussian kernel based on the distance between locations and the specified neighborhood parameters (\code{npoints}, \code{curweight}, and \code{zetaScale}).
#' The neighborhood for the underdamped model is \code{zetaScale * sqrt(2*pi) * sigma}, while the neightborhood for the overdamped model is \code{zetaScale * sigma / sqrt(2)}. See Blackwell and Matthiopoulos (2024) for more details.
#'
#' @template barrier_details
#'
#' @references
#' Dupont F, McClintock BT, Fischer J-O, Marcoux M, Hussey N, Auger-Methe M. 2025. Inferring resource selection and utilization distributions from irregular and error-prone animal tracking data using the habitat-driven Langevin diffusion.
#'
#' McClintock BT, London JM, Camneron MF, Boveng PL. 2015. Modelling animal movement using the Argos satellite telemetry location error ellipse. Methods Ecol Evol 6:266–277. doi: 10.1111/2041-210X.12286.
#'
#' Michelot T, Gloaguen P, Blackwell PG, Etienne M-P. 2019. The Langevin diffusion as a continuous-time model of animal movement and habitat selection. Methods Ecol Evol 10:1894–1907. doi: 10.1111/2041-210X.13275.
#'
#' Michelot T, Hanks E. 2025. Multiscale modelling of animal movement with persistent dynamics. arXiv doi: 10.48550/arXiv.2406.15195
#'
#' Blackwell PG, Matthiopoulos J. 2024. Joint inference for telemetry and spatial survey data. Ecology 105(12):e4457. doi: 10.1002/ecy.4457.
#'
#' @examples
#' # fit underdamped model with measurement error
#' # exampleDat included in package; see ?exampleDat for details
#' # exampleCovs included in package; see ?exampleCovs for details
#' fit <- fitLangevin(data = exampleDat,
#'                    spatialCovs = exampleCovs,
#'                    silent = TRUE)
#'
#' # fit overdamped model with measurement error
#' fit_od <- fitLangevin(data = exampleDat,
#'                       model = "overdamped",
#'                       spatialCovs = exampleCovs,
#'                       silent = TRUE)
#'
#' \dontrun{
#' # simulating with a barrier and passing the penalty to fitLangevin
#'
#' par <- list(beta = c(-4, 6, 5, -0.1),
#'             sigma = 5,
#'             gamma = 0.5)
#'
#' # create a dummy barrier mask (left half restricted = 0, right half allowed = 1)
#' coast_barrier <- exampleCovs[[1]]
#' terra::values(coast_barrier) <- ifelse(terra::crds(coast_barrier)[, "x"]
#'                                  >= mean(terra::crds(coast_barrier)[, "x"]), 1, 0)
#' names(coast_barrier) <- "coast_barrier"
#'
#' # convert mask to SDF and add to the spatial covariates list
#' exampleCovs_barrier <- exampleCovs
#' exampleCovs_barrier$coast_barrier <- prepBarrier(coast_barrier)
#'
#' # add a beta coefficient for the barrier to the parameter list
#' par_barrier <- par
#' par_barrier$beta <- c(par_barrier$beta, -0.2)
#'
#' # simulate the data
#' set.seed(123,kind="Mersenne-Twister",normal.kind="Inversion")
#' simDat_barrier <- simLangevin(par = par_barrier,
#'                               spatialCovs = exampleCovs_barrier,
#'                               measurementError = list(smaj.sd = 1.5,
#'                                                       smin.sd = 0.75,
#'                                                       eor.lim = c(0,180)))
#'
#' # Because simDat_barrier is a simLangevin object, fitLangevin will automatically
#' # detect and use the exact barrier penalty (lambda) that generated the data
#' fit_barrier <- fitLangevin(data = simDat_barrier,
#'                            spatialCovs = exampleCovs_barrier,
#'                            silent = TRUE)
#'
#' plot(fit_barrier,data=simDat_barrier,spatialCovs = exampleCovs_barrier)
#' }
#'
#' @rawNamespace useDynLib(langevinSSM, .registration=TRUE); useDynLib(langevinSSM_TMBExports)
#' @importFrom stats nlminb
#' @importFrom TMB MakeADFun sdreport oneStepPredict
#' @export
fitLangevin <- function(data, model = c("underdamped","overdamped"), spatialCovs, par, prior = NULL, map=NULL, coord = c("x", "y"), lambda = NULL, scaleFactor = 1, smoothGradient = FALSE, npoints = 4, curweight = 0.5, zetaScale = 1, hessian=FALSE, silent=FALSE, method="BFGS", initialInner = TRUE, inner.control=list(maxit=1000), control = list(trace=0,iter.max=1000,eval.max=1000), polishOptim = FALSE, getJointPrecision = FALSE){

  if(!inherits(data,"dataLangevin")) stop("'data' is not formatted as a 'dataLangevin' object. See ?formatData")
  model <- match.arg(model)

  if(!all(coord %in% colnames(data))) stop("coord not found in data.")
  if(smoothGradient & isFALSE(npoints %in% c(4,8))) stop("npoints must be 4 or 8")
  if(smoothGradient & isFALSE(curweight>=0 & curweight<1)) stop("curweight must be >=0 and <1")
  if(smoothGradient & isFALSE(zetaScale>0)) stop("zetaScale must be >0")

  checkErrorData(data, coord)

  barrier <- .find_barrier(spatialCovs)
  if (!is.null(barrier)) {
    if (is.null(lambda)) {
      sim_lambda <- attr(data, "lambda")
      if (!is.null(sim_lambda)) {
        lambda <- sim_lambda
        message("   Auto-detected barrier penalty (lambda) from simLangevin data: ", signif(lambda, 4))
      } else {
        stop("A numeric value for 'lambda' must be supplied when using a barrier constraint.\n",
             "  If you do not know what value to use, see ?fitLangevin_barrier.")
      }
    }
    .validate_lambda(lambda)
  }

  dat <- build_tmb_data(data, spatialCovs, model, coord, scaleFactor, smoothGradient, npoints, curweight, zetaScale, barrier_names = barrier, lambda = lambda)

  par <- initialValues(data, model, par, spatialCovs, coord)
  cp <- checkPar(par, model, map, dat=dat, spatialCovs=spatialCovs, prior=prior)
  par <- cp$par
  map <- cp$map
  re <- cp$re

  dat <- c(dat, cp$priors)
  map <- mapDuplicatedTimes(dat, map, par, re)

  # scale parameters internally for TMB
  par$log_sigma <- par$log_sigma - log(scaleFactor)
  par$mu <- par$mu / scaleFactor
  par$vel <- par$vel / scaleFactor

  message("   Fitting ",model," Langevin model...")

  if(initialInner){
    map_inner <- lapply(par[names(par)[!names(par) %in% c("mu","vel")]], function(x) factor(rep(NA,length(x))))
    if(!is.null(map$mu)) map_inner$mu <- map$mu
    if(!is.null(map$vel)) map_inner$vel <- map$vel

    obj1 <- try({
      TMB::MakeADFun(c(model="langevinSSM",dat), par, map = map_inner, random = re, DLL = "langevinSSM_TMBExports", hessian = hessian, method = method, silent = silent, inner.control = inner.control)
    }, silent = TRUE)

    if (inherits(obj1, "try-error")) stop("Initial inner optimization (obj1) failed during TMB::MakeADFun. Check parameter initial values or data.\nError details: ", attr(obj1, "condition")$message)

    obj1$fn()
    smoothed_pars <- obj1$env$parList()

    if("mu" %in% re) par$mu <- smoothed_pars$mu
    if("vel" %in% re) par$vel <- smoothed_pars$vel

    if(!is.null(map$mu)) par$mu[is.na(map$mu)] <- (cp$par$mu / scaleFactor)[is.na(map$mu)]
    if(!is.null(map$vel)) par$vel[is.na(map$vel)] <- (cp$par$vel / scaleFactor)[is.na(map$vel)]
  }

  obj2 <- try({
    TMB::MakeADFun(c(model="langevinSSM",dat), par, map = map, random = re, DLL = "langevinSSM_TMBExports", hessian = hessian, method = method, silent = silent, inner.control = inner.control)
  }, silent = TRUE)

  if (inherits(obj2, "try-error")) stop("TMB::MakeADFun failed to construct the objective function. Check parameter initial values or data constraints.\nError details: ", attr(obj2, "condition")$message)

  if (length(obj2$par) == 0) {
    start <- proc.time()
    fit <- list(par = obj2$par, objective = obj2$fn(obj2$par), convergence = 0, message = "All parameters fixed; no outer optimization required", iterations = 0, evaluations = c("function" = 1, "gradient" = 0))
    fit$elapsedTime <- proc.time() - start
  } else {
    start <- proc.time()
    fit <- try({
      do.call(stats::nlminb,args = list(start = obj2$par, objective = obj2$fn, gradient = obj2$gr, control=control))
    }, silent = TRUE)

    if (inherits(fit, "try-error") || !is.list(fit)) stop("Optimization via stats::nlminb failed. The model could not be fit.\nError details: ", attr(fit, "condition")$message)

    if(polishOptim==TRUE){
      obj2$fn(fit$par)
      message("   Polishing optimization with nlminb...")
      fit_polished <- try({
        do.call(stats::nlminb, args = list(start = fit$par, objective = obj2$fn, gradient = obj2$gr, control = control))
      }, silent = TRUE)
      if (!inherits(fit_polished, "try-error") && fit_polished$objective < fit$objective) fit <- fit_polished
    }

    fit$elapsedTime <- proc.time() - start
    if (!is.null(fit$convergence) && fit$convergence != 0) warning("nlminb optimization did not appear to converge. Code: ", fit$convergence, " - ", fit$message)
  }

  message("   Calculating SEs...")
  sdreport_out <- try({
    TMB::sdreport(obj2, getJointPrecision = getJointPrecision)
  }, silent = TRUE)

  extracted <- extract_tmb_estimates(fit, obj2, sdreport_out, re, map, data, scaleFactor, spatialCovs, getJointPrecision)
  fit$estimates <- extracted$estimates
  fit$covariance <- extracted$covariance

  out_of_bounds <- FALSE
  if (!is.null(fit$estimates$random$mu)) {
    mu_est <- fit$estimates$random$mu$est
    cov_ext <- as.vector(terra::ext(spatialCovs[[1]]))
    cov_res <- terra::res(spatialCovs[[1]])

    safe_xmin <- cov_ext["xmin"] + cov_res[1]; safe_xmax <- cov_ext["xmax"] - cov_res[1]
    safe_ymin <- cov_ext["ymin"] + cov_res[2]; safe_ymax <- cov_ext["ymax"] - cov_res[2]

    out_of_bounds <- any(mu_est$mu.x < safe_xmin | mu_est$mu.x > safe_xmax |
                           mu_est$mu.y < safe_ymin | mu_est$mu.y > safe_ymax, na.rm = TRUE)
  }

  fit$conditions <- list(hessian = hessian, method = method, silent = silent, initialInner = initialInner, inner.control = inner.control, control = control, scaleFactor = scaleFactor, model = model, smoothGradient = smoothGradient, npoints = npoints, curweight = curweight, zetaScale = zetaScale, coord = coord, out_of_bounds = out_of_bounds, barrier = barrier, lambda = lambda)

  boundsWarning(fit)

  fit$signatures <- list(data = get_data_signature(data, coord), covs = get_covs_signature(spatialCovs))
  fit$tmb_setup <- list(parList = obj2$env$parList(fit$par), map = map, random = re, priors = cp$priors)

  fit <- class_fitLangevin(fit)
  return(fit)
}
