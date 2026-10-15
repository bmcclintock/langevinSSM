build_tmb_data <- function(data, spatialCovs, model, coord, scaleFactor, smoothGradient, npoints, curweight, zetaScale, obs_dist, barrier_sdf = NULL, lambda = NULL, lambda_mode = 0L, grad = "quintic", smoothPar = NULL, barrierOnset = 0) {
  if (!is.null(barrier_sdf)) {
    barrier_dist_mat <- terra::as.matrix(barrier_sdf, wide = TRUE) / scaleFactor
    # lambda is on the natural scale; the template works in coordinates / scaleFactor (lambda * scaleFactor^2)
    barrier_pen <- if (lambda_mode == 0L) lambda * scaleFactor^2 else 0
  } else {
    lambda_mode <- 0L
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
    obs_dist = ifelse(obs_dist == "t", 1, 0),
    Y = t(data[, coord]) / scaleFactor,
    times = track_times,
    delta_t = data$dt
  )

  dat$skip_step <- as.integer(dat$delta_t < 1.e-6)
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
  # fixed (sigma, gamma) for the smoothing radius, natural scale -> working scale
  if (is.null(smoothPar)) smoothPar <- list(sigma = 1, gamma = 1)
  dat$smooth_sigma <- smoothPar$sigma / scaleFactor
  dat$smooth_gamma <- smoothPar$gamma

  dat$barrier_dist <- barrier_dist_mat
  dat$barrier_penalty <- barrier_pen
  dat$lambda_mode <- as.integer(lambda_mode)
  # onset length of the barrier potential (natural units -> working units)
  dat$barrier_onset <- if (is.null(barrier_sdf)) 0 else .resolve_onset(barrierOnset, barrier_sdf) / scaleFactor
  pos_dt <- dat$delta_t[dat$delta_t > 1.e-6]
  dat$median_dt <- if (length(pos_dt)) stats::median(pos_dt) else 1
  dat$spline_degree <- .grad_degree(grad)
  dat <- c(dat, raster_data)

  return(dat)
}

# --- Helper: Extract and Format TMB Estimates ---
extract_tmb_estimates <- function(fit, obj, sdreport_out, re, map, data, scaleFactor, spatialCovs, getJointPrecision) {
  estimates <- list()
  covariance <- list()

  if (inherits(sdreport_out, "try-error") || inherits(sdreport_out, "skip_se")) {
    if (inherits(sdreport_out, "try-error")) {
      warning("TMB::sdreport failed to calculate standard errors. Trajectory and point estimates were recovered from the report, but SEs are unavailable.")
    }

    # 1. working scale
    working_names <- names(fit$par)
    name_counts <- table(working_names)
    is_dup <- working_names %in% names(name_counts[name_counts > 1])
    if (any(is_dup)) {
      suffix <- ave(working_names, working_names, FUN = seq_along)
      working_names[is_dup] <- paste0(working_names[is_dup], "_", suffix[is_dup])
    }
    estimates$working <- data.frame("Estimate" = as.numeric(fit$par), "Std. Error" = rep(NA_real_, length(fit$par)), check.names = FALSE)
    rownames(estimates$working) <- working_names

    # 2. natural scale
    rep_vals <- obj$report()
    nat_names <- c("beta", "sigma", "gamma", "rho_o", "tau", "psi","df")
    if (isTRUE(obj$env$data$lambda_mode > 0)) nat_names <- c(nat_names, "lambda")
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

    if("All parameters fixed; no outer optimization required" %in% fit$message){
      estimates$working <- as.data.frame(suppressWarnings(summary(sdreport_out, "fixed")))
    } else {
      estimates$working <- as.data.frame(summary(sdreport_out, "fixed"))
    }

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
#' @param obs_dist Character string specifying the measurement error observation distribution. Can be ``Gaussian'' (default) or ``t''.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates.
#' @param barrier Character string. The name of the barrier in \code{spatialCovs} that is represented as a signed distance field (see \code{\link{prepBarrier}}). If provided, this raster is exclusively used for the barrier penalty and is not included in the habitat selection covariates. Default: \code{NULL} (no barrier).
#' @param par List containing the initial values for the parameters. These can include state process parameters for the habitat selection coefficients (``beta''), the ``speed'' parameter (``sigma''), and, for the underdamped model, the friction coefficient (``gamma''). If \code{lambda = "estimate"}, \code{par} can also include the initial barrier strength (``lambda'', natural scale; default: the stability limit used by \code{lambda = "tied"} at the median time step, see \code{\link{initialValues}}), which can be fixed using \code{map = list(lambda = factor(NA))}. ``lambda'' cannot be included in \code{par} otherwise.
#' Observation process parameters include a scaling factor to account for uncertainty in the Argos error ellipse (``psi''), a scaling factor to account for uncertainty in the x- and y-axis errors for Argos least squares or GPS observations (``tau''), and a correlation term between the x- and y-axis errors for Argos least squares or GPS observations (``rho_o''). If \code{obs_dist = "t"}, observation process parameters also include the degrees of freedom (``df''). All parameters are specified on their natural scale and are converted to working scale internally.
#' Any missing state process parameters are generated using \code{\link{initialValues}}. Any missing observation process parameters are fixed to their default values (``psi'' = 1, ``tau'' = 1, ``rho_o'' = 0, and ``df'' = 5) via \code{map}. See Details.
#' @param lambda Barrier strength. One of: \code{"tied"} (default), a non-negative number (fixed), or \code{"estimate"}. With \code{"tied"} the barrier is treated as a hard wall: at every time step \eqn{t} the spring constant is set to the largest value for which the one-step transition is numerically stable, \eqn{\lambda_t = \lambda_{max}(\sigma, \gamma, \Delta t_t)}, computed inside the likelihood from the current \code{sigma} and \code{gamma} (see Details), so \eqn{\lambda} is not a free parameter and the wall becomes harder automatically as the time step shrinks. A number fixes \eqn{\lambda} (a soft wall of known stiffness, e.g. data simulated by \code{\link{simLangevin}}); \code{"estimate"} treats \eqn{\log\lambda} as a free parameter (only identifiable when animals actually enter the restricted zone), with its starting value given by \code{par$lambda}. Numeric values, \code{par$lambda}, and the reported estimates are all on the natural scale (units of 1/distance^2 in the units of the coordinates) and do not depend on \code{scaleFactor}. If \code{NULL}, the value stored by \code{\link{simLangevin}} is used when available, otherwise \code{"tied"}.
#' @param barrierOnset Onset length \eqn{\ell} of the barrier potential, in the units of the coordinates (or \code{"res"} for one cell of the barrier raster). With \code{0} the potential is \eqn{B(d) = \min(d,0)^2/2} (the penalty of Delporte et al. 2025), whose force is continuous at the boundary but whose curvature jumps from 0 to \eqn{\lambda}. With \eqn{\ell > 0} the force rises smoothly (\eqn{C^4}) from zero at the boundary over the first \eqn{\ell} of the restricted side and is thereafter the quadratic wall shifted by \eqn{\ell/2}; the allowed side is unaffected. A small positive onset (e.g., one cell) makes the Laplace approximation better behaved when many (smoothed) locations lie close to the boundary, and is a natural way to encode positional uncertainty of the boundary itself. It is not a fixed quantity and should be chosen for the species and data: for animals that should never be on the restricted side (e.g., whales at a coastline) it is a numerical device of the order of the raster resolution or coastline accuracy; for animals that legitimately use the restricted side (e.g., pinnipeds hauling out on land) a larger value (or a soft, estimated \code{lambda}) is more appropriate. Ignored without a barrier. If \code{NULL} (default), the value stored by \code{\link{simLangevin}} is used when available, otherwise \code{0}. See Details.
#' @param prior Optional 2-column data frame containing the mean (column 1) and standard deviation (column 2) for normally distributed priors on the working scale parameters. The row names must match the working scale parameter names (e.g., \code{"beta_cov1"}, \code{"log_sigma"}, \code{"log_gamma"}, or \code{"log_lambda"} when \code{lambda = "estimate"}). Supplying the base name of a vector or matrix parameter (e.g., \code{"beta"}, \code{"mu"}, \code{"vel"}) will apply the prior to all of its elements. To target specific coordinates and time steps for the random effects, append the coordinate (\code{.x} or \code{.y}) and the row index of the observation to the base name (e.g., \code{"mu.x_1"} for the x-coordinate of the 1st observation in \code{data}, or \code{"vel.y_10"} for the y-velocity of the 10th observation in \code{data}). Parameters omitted from this data frame are assigned flat (improper) priors. Default: \code{NULL} (no priors).
#' @param map List defining how to optionally collect and fix parameters. See \code{\link[TMB]{MakeADFun}}.
#' @param coord Character vector identifying the coordinate names for the location data. Default: \code{c("x","y")}.
#' @param scaleFactor Internal scaling factor for the coordinates and parameters. In some cases, setting \code{scaleFactor>1} can help with optimization.
#' @param grad Interpolation used for the covariate surfaces (and the barrier distance) whose gradients give the habitat force. One of \code{"quintic"} (default), \code{"cubic"} or \code{"bilinear"}. See Details.
#' @param smoothGradient Logical indicating whether or not to smooth the gradients. See Details. Default: \code{FALSE}.
#' @param npoints Number of smoothing points around current cell (4 = diagonal, 8 = queen neighborhood). Ignored unless \code{smoothGradient=TRUE}.
#' @param curweight Smoothing weight of current cell location (\code{0 <= curweight < 1}). Ignored unless \code{smoothGradient=TRUE}.
#' @param zetaScale Multiplier of the smoothing radius (\code{zetaScale>1} increases and \code{zetaScale<1} decreases the neighborhood). \code{zetaScale = 1} (default) is the moment-matched radius for both models, evaluated at values of \code{sigma} and \code{gamma} derived from the data before fitting (see Details). Ignored unless \code{smoothGradient=TRUE}.
#' @param hessian Logical indicating whether or not to calculate the Hessian at the optimum. See \code{\link[TMB]{MakeADFun}}. Default: \code{FALSE}.
#' @param silent Logical indicating whether or not to disable TMB tracing information. See \code{\link[TMB]{MakeADFun}}. Default: \code{TRUE}.
#' @param method Character string indicating the type of algorithm used for the outer optimization (e.g., \code{"BFGS"}). When \code{optMethod = "nloptr"}, \code{"BFGS"} is internally mapped to \code{"NLOPT_LD_LBFGS"} for compatibility. Default: \code{"BFGS"}.
#' @param optMethod Character string indicating the optimization function (engine) to use for the outer optimization. Options are \code{"nlminb"} or \code{"nloptr"}. If \code{"nloptr"} is selected, the \code{nloptr} package must be installed. Default: \code{"nlminb"}.
#' @param inner.control List controlling inner optimization. See \code{\link[TMB]{MakeADFun}}.
#' @param control A list of control parameters for the outer optimization. For \code{optMethod="nlminb"}, see \code{\link[stats]{nlminb}}. For \code{optMethod="nloptr"}, these arguments are directly mapped to the \code{opts} list in \code{\link[nloptr]{nloptr}}.
#' @param getJointPrecision Logical indicating whether or not to return the joint precision matrix for the random effects. Default: \code{FALSE}.
#' @param calcSE Logical indicating whether or not to calculate standard errors for the parameter estimates. Default: \code{TRUE}.
#'
#' @return \code{fitLangevin} object, i.e., a list of:
#' \item{par}{See \code{\link[stats]{nlminb}} or \code{\link[nloptr]{nloptr}}}
#' \item{objective}{See \code{\link[stats]{nlminb}} or \code{\link[nloptr]{nloptr}}}
#' \item{convergence}{See \code{\link[stats]{nlminb}} or \code{\link[nloptr]{nloptr}}}
#' \item{message}{See \code{\link[stats]{nlminb}} or \code{\link[nloptr]{nloptr}}}
#' \item{iterations}{See \code{\link[stats]{nlminb}} or \code{\link[nloptr]{nloptr}}}
#' \item{evaluations}{See \code{\link[stats]{nlminb}} or \code{\link[nloptr]{nloptr}}}
#' \item{elapsedTime}{Run time of the optimization}
#' \item{estimates}{List containing point estimates and standard errors for the natural scale parameters (``natural''), the working scale parameters (``working''), and the random effects (``random''), where ``random'' is itself a list containing point estimates (``est'') and standard errors (``se'') for the true locations (``mu'') and/or the true velocities (``vel'').}
#' \item{covariance}{List containing the covariance matrices for the natural scale parameters (``natural''), the working scale parameters (``working''), and, if \code{getJointPrecision=TRUE}, the joint covariance matrix for the random effects (``random''). If \code{hessian=TRUE}, it also returns the exact \code{hessian} matrix computed via automatic differentiation.}
#' \item{conditions}{List containing the optimization settings}
#' \item{signatures}{List containing lightweight fingerprints for \code{data} and \code{spatialCovs} to protect downstream functions}
#' \item{tmb_setup}{Blueprint for reconstructing TMB objective function}
#'
#' @details
#' \strong{Measurement Error Models:}
#' The \code{fitLangevin} function accommodates various types of measurement error depending on the provided measurement error data. It handles measurement error through two distinct internal models with the same underlying structure, which differ in how the observation error covariance matrix is constructed based on the available error information. The observations can be modeled using either a bivariate Gaussian distribution (\code{obs_dist = "Gaussian"}):
#' \deqn{(x,y) \sim N((\mu_x,\mu_y), \Sigma),}
#' or a bivariate Student's t-distribution (\code{obs_dist = "t"}) to better handle heavy-tailed observation outliers:
#' \deqn{(x,y) \sim t_{\nu}((\mu_x,\mu_y), \Sigma),}
#' where \eqn{(x,y)} are the observed locations, \eqn{(\mu_x,\mu_y)} are the true (latent) locations, \eqn{\Sigma} is the observation error covariance matrix, and \eqn{\nu} represents the degrees of freedom (\code{df} parameter).
#' \itemize{
#'   \item \strong{Error Ellipse Model (Argos Kalman Filter):} This model is used when error ellipse data (\code{smaj}, \code{smin}, and \code{eor}) are present for an observation. The (optional) observation process parameter \code{psi} acts as a scaling multiplier specifically on the semi-minor axis (\code{smin}) to account for uncertainty in the error ellipse.
#'   For the error ellipse model, \eqn{\Sigma} is derived from the semi-major axis, semi-minor axis, and error ellipse orientation (see McClintock et al. 2015), with \eqn{\psi} (\code{psi}) scaling the semi-minor axis:
#'   \deqn{\Sigma = \begin{bmatrix} (M/\sqrt{2})^2 \sin^2(c) + (\psi m/\sqrt{2})^2 \cos^2(c) & \frac{1}{2}\left(M^2 - (\psi m)^2\right) \cos(c) \sin(c) \\ \frac{1}{2}\left(M^2 - (\psi m)^2\right) \cos(c) \sin(c) & (M/\sqrt{2})^2 \cos^2(c) + (\psi m/\sqrt{2})^2 \sin^2(c) \end{bmatrix}}
#'   where \eqn{M} is the semi-major axis (\code{smaj}), \eqn{m} is the semi-minor axis (\code{smin}), and \eqn{c} is the ellipse orientation in radians from north (\code{eor}).
#'   \item \strong{Standard Deviation Model (Argos Least Squares, GPS, Generic Locations):} This model is used when x- and y-axis standard deviations (\code{x.err} and \code{y.err}) are present for an observation. The (optional) observation parameter \code{tau} (a 2-vector) scales the \code{x.err} and \code{y.err}, respectively, to account for uncertainty in the x- and y-axis standard deviations, while the (optional) parameter \code{rho_o} accounts for correlation between the x- and y-axis standard deviations.
#'   For the standard deviation model:
#'   \deqn{\Sigma = \begin{bmatrix} \tau_1^2 \sigma_x^2 & \rho_o \tau_1 \tau_2 \sigma_x \sigma_y \\ \rho_o \tau_1 \tau_2 \sigma_x \sigma_y & \tau_2^2 \sigma_y^2 \end{bmatrix}}
#'   with \eqn{\text{tau}=(\tau_1,\tau_2)} scaling the standard deviations \eqn{(\text{x.err}=\sigma_x,\text{y.err}=\sigma_y)} and rho_o \eqn{=\rho_o} accounting for the correlation between the x- and y-axis errors.
#' }
#' Whenever a coordinate is \code{NA} (either present in the original dataset or introduced via the \code{predTimes} argument in \code{\link{formatData}}), \code{fitLangevin} skips the evaluation of the observation error model for that time step. Instead, it relies solely on the Langevin diffusion state process model to predict the true (latent) location.
#'
#' \strong{Observation Process Parameters (\code{psi}, \code{tau}, \code{rho_o}, \code{df}):}
#' These parameters control the scaling and correlation of measurement errors and can be specified via \code{par} and \code{map}:
#' \itemize{
#'   \item \strong{1. Fixed at defaults:} If observation process parameters are omitted from \code{par}, the model automatically fixes them to neutral values: \code{psi = 1} (no scaling semi-minor axis), \code{tau = c(1, 1)} (no scaling of \code{x.err} or \code{y.err}), and \code{rho_o = 0} (no correlation between x- and y-axis errors). When \code{obs_dist="t"}, the default is a robust, heavy-tailed \code{df = 5}.
#'   \item \strong{2. Estimated from data:} To estimate one or more of these parameters, include them in \code{par} with an initial starting value (e.g., \code{par = list(psi = 1.2)}).
#'   \item \strong{3. Fixed to custom values:} To fix a parameter to a specific value other than the default, include the custom value in \code{par} and explicitly map it to \code{NA} in \code{map} (e.g., \code{par = list(psi = 1.5), map = list(psi = factor(NA))}).
#' }
#'
#' \strong{Gradient interpolation (\code{grad}):} the force is the gradient of the covariate surfaces, which are interpolated from the raster cell-centre values by an interpolating tensor-product B-spline of degree 5 (\code{"quintic"}), 3 (\code{"cubic"}) or 1 (\code{"bilinear"}). All three pass exactly through the cell values, and in all three the raster cell is looked up afresh from the current location at every evaluation of the objective. They differ in smoothness, which matters for the Laplace approximation: the inner (latent-path) problem uses the 2nd and 3rd derivatives of the surface, and TMB's gradient of the marginal likelihood the 4th. With \code{"quintic"} all of these are continuous, so the marginal likelihood and its gradient are smooth. With \code{"cubic"} the force and its first derivative are continuous but the 3rd derivative of the surface jumps at cell edges, so the marginal likelihood is slightly non-smooth; the inner Hessian is cheaper. With \code{"bilinear"} the force itself jumps at every cell edge and its derivatives within a cell are zero, which makes each inner step cheap but the inner problem non-smooth; the outer optimizer can then report false convergence and stop at different points from different starting values. \code{"bilinear"} is mainly useful for quick exploratory fits or comparison with earlier versions of the package; check its results against \code{"cubic"} or \code{"quintic"}. The estimates (in particular \code{beta}) depend somewhat on the interpolation, because the gradient of a rough raster is not uniquely defined; \code{simLangevin}, \code{predLangevin}, \code{cAIC} and \code{residuals} use the setting stored in the fit.
#'
#' \strong{Gradient smoothing:} Setting \code{smoothGradient=TRUE} applies the smoothing approach of Blackwell and Matthiopoulos (2024) to the gradients of the habitat selection covariates (\code{spatialCovs}), which can help reduce attenuation bias in habitat selection coefficients when observations are obtained at a coarser time resolution than the underlying continuous-time movement process.
#' This approach smooths the gradient at each location by taking a weighted average of the gradients at neighboring locations, with weights and neighborhood determined by the parameters (\code{npoints}, \code{curweight}, and \code{zetaScale}).
#' The \code{npoints} neighbours lie on a ring of radius \eqn{r} around the current location (equally spaced; weight \code{(1-curweight)/npoints} each, \code{curweight} at the centre), so the smoothed gradient is a discrete average of the gradient over the positions the animal may occupy during the interval.
#' The transition density holds the force fixed at its value at the start of each interval \eqn{\Delta}, whereas in the continuous-time process it acts along the path; to second order the force averaged over the interval (weighted as in the transition mean) is \eqn{h(x) + \frac{1}{2} m \nabla^2 h(x)}, where \eqn{m} is the correspondingly weighted mean squared displacement per coordinate.
#' The ring has second moment \eqn{(1-c) r^2/2} per coordinate (\eqn{c} = \code{curweight}) for both \code{npoints = 4} and \code{npoints = 8}, so it reproduces this correction when \eqn{r^* = \sqrt{2m/(1-c)}}.
#' For the underdamped model, with \eqn{a = \gamma\Delta}, \deqn{m = \frac{2\sigma^2}{\gamma^2}\,\frac{a^2/2 - 2a + 3 - (3+a)e^{-a}}{a - 1 + e^{-a}},}
#' so that \eqn{r^* \approx \sigma\Delta/\sqrt{3(1-c)}} when \eqn{\gamma\Delta \ll 1} (ballistic motion within the interval) and \eqn{r^* \approx \sigma\sqrt{2\Delta/(\gamma(1-c))}} when \eqn{\gamma\Delta \gg 1}. For the overdamped model \eqn{m = \sigma^2\Delta/2} and \eqn{r^* = \sigma\sqrt{\Delta/(1-c)}}; for the default \code{curweight = 0.5} this is \eqn{\sigma\sqrt{2\Delta}}, the root-mean-square two-dimensional step length, i.e. exactly the radius \eqn{\zeta\sqrt{2\delta t}} of Blackwell and Matthiopoulos (2024, eqs. 29--30).
#' The radius used is \eqn{r = } \code{zetaScale} \eqn{\times\, r^*} (\code{zetaScale = 1}, the default, is the moment-matched radius for both models), computed for each interval \eqn{\Delta}. As in Blackwell and Matthiopoulos (2024), where \eqn{\zeta} is computed from the observed step lengths before fitting, \eqn{\sigma} and \eqn{\gamma} in \eqn{r^*} are fixed values derived from the data alone: the measurement-error-corrected variogram estimates of \code{\link{initialValues}} (\code{type = "variogram"}), or, if the variogram fit fails, the moment estimators \eqn{\sigma = \sqrt{\overline{|\Delta x|^2/(2\Delta)}}} (B&M's \eqn{\zeta}) and \eqn{\gamma = 1/\mathrm{median}(\Delta)}. They do not depend on \code{par} or on the parameters being estimated, so the smoothing is a fixed property of the likelihood. The values used are stored in \code{fit$conditions$smoothPar}. For the underdamped model the radius is insensitive to \eqn{\gamma} when \eqn{\gamma\Delta \ll 1}; to adjust the radius, use \code{zetaScale}.
#' The ring averages over the stationary velocity distribution rather than conditioning on the latent velocity, and it corrects only the mean of the transition (not its variance), so it mainly addresses attenuation from abrupt changes in the gradient over the interval; it cannot remove the bias that arises when the observation interval is long relative to the spatial scale of the covariates. See Blackwell and Matthiopoulos (2024) for the smoothing approach.
#'
#' @template barrier_details
#'
#' @references
#' Dupont F, McClintock BT, Fischer J-O, Marcoux M, Hussey N, Auger-Methe M. 2025. Inferring resource selection and utilization distributions from irregular and error-prone animal tracking data using the habitat-driven Langevin diffusion.
#'
#' Delporte, A., Ditlevsen, S. and Samson, A., 2025. Spatial constraints improve filtering of measurement noise from animal tracks. arXiv preprint arXiv:2511.22430.
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
#'                    spatialCovs = exampleCovs)
#'
#' # fit overdamped model with measurement error
#' fit_od <- fitLangevin(data = exampleDat,
#'                       model = "overdamped",
#'                       spatialCovs = exampleCovs)
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
#' # Because simDat_barrier is a simLangevin object, fitLangevin will automatically
#' # detect and use the exact barrier penalty (lambda) that generated the data
#' fit_barrier <- fitLangevin(data = simDat_barrier,
#'                            spatialCovs = exampleCovs_barrier)
#'
#' plot(fit_barrier, data = simDat_barrier,
#'                   spatialCovs = exampleCovs_barrier,
#'                   maskRast = coast_barrier)
#' }
#'
#' @rawNamespace useDynLib(langevinSSM, .registration=TRUE); useDynLib(langevinSSM_TMBExports)
#' @importFrom stats nlminb
#' @importFrom TMB MakeADFun sdreport oneStepPredict
#' @export
fitLangevin <- function(data, model = c("underdamped","overdamped"), obs_dist = c("Gaussian", "t"), spatialCovs, barrier = NULL, par, lambda = "tied", barrierOnset = NULL, prior = NULL, map=NULL, coord = c("x", "y"), scaleFactor = 1, grad = c("quintic", "cubic", "bilinear"), smoothGradient = FALSE, npoints = 4, curweight = 0.5, zetaScale = 1, hessian=FALSE, silent=TRUE, method="BFGS", optMethod=c("nlminb", "nloptr"), inner.control=list(maxit=1000), control = list(trace=0,iter.max=1000,eval.max=1000), getJointPrecision = FALSE, calcSE = TRUE){

  if(!inherits(data,"dataLangevin")) stop("'data' is not formatted as a 'dataLangevin' object. See ?formatData")

  data <- class_dataLangevin(data)

  model <- match.arg(model)
  obs_dist <- match.arg(obs_dist)
  optMethod <- match.arg(optMethod)

  if (optMethod == "nloptr" && !requireNamespace("nloptr", quietly = TRUE)) {
    stop("Package 'nloptr' is required when optMethod = 'nloptr'. Please install it using install.packages('nloptr').")
  }

  if(!all(coord %in% colnames(data))) stop("coord not found in data.")
  grad <- match.arg(grad)
  if(smoothGradient & isFALSE(npoints %in% c(4,8))) stop("npoints must be 4 or 8")
  if(smoothGradient & isFALSE(curweight>=0 & curweight<1)) stop("curweight must be >=0 and <1")
  if(smoothGradient & isFALSE(zetaScale>0)) stop("zetaScale must be >0")

  checkErrorData(data, coord)

  if(is.null(barrier)){
    sim_barrier <- attr(data,"barrier")
    if(!is.null(sim_barrier)){
      barrier <- sim_barrier
      message("   Auto-detected barrier from simLangevin data: ", barrier)
    }
  }

  orig_spatialCovs <- spatialCovs

  if (!is.null(barrier)) {
    if (!(barrier %in% names(spatialCovs))) stop(sprintf("Barrier raster '%s' not found in spatialCovs.", barrier))

    barrier_sdf <- spatialCovs[[barrier]]
    if(!isTRUE(attr(barrier_sdf,"barLangevin"))) stop("barrier is not a 'barLangevin' object created by prepBarrier.")
    spatialCovs[[barrier]] <- NULL

    lambda_mode <- 0L
    if (is.null(lambda)) {
      sim_lambda <- attr(data, "lambda")
      if (!is.null(sim_lambda)) {
        lambda <- if (is.character(sim_lambda) || is.infinite(sim_lambda)) "tied" else sim_lambda
        message("   Auto-detected barrier penalty (lambda) from simLangevin data: ", if (is.character(lambda)) lambda else signif(lambda, 4))
      } else lambda <- "tied"
    }
    if (is.character(lambda)) {
      lambda <- match.arg(lambda, c("tied", "estimate"))
      lambda_mode <- if (lambda == "tied") 1L else 2L
    }
    if (lambda_mode == 0L) .validate_lambda(lambda)
    if (is.null(barrierOnset)) {
      sim_onset <- attr(data, "barrierOnset")
      barrierOnset <- if (!is.null(sim_onset)) sim_onset else 0
      if (!is.null(sim_onset) && sim_onset > 0) message("   Auto-detected barrier onset length from simLangevin data: ", signif(sim_onset, 4))
    }
    barrierOnset <- .resolve_onset(barrierOnset, barrier_sdf)
  } else {
    barrier_sdf <- NULL
    lambda_mode <- 0L
    barrierOnset <- 0
  }

  if (length(spatialCovs) == 0) {
    stop("At least one habitat covariate must remain in spatialCovs after isolating the barrier constraint.")
  }

  # Smoothing radius (Blackwell & Matthiopoulos 2024): fixed before optimisation, derived from the data only
  smoothPar <- if (!smoothGradient) NULL else if (!is.null(attr(data, "smoothPar"))) attr(data, "smoothPar") else
    .smooth_par_data(data, model, obs_dist, spatialCovs, coord)   # attribute: set internally by predLangevin()
  par <- initialValues(data, model, obs_dist, par, spatialCovs, coord, lambda = if (lambda_mode == 2L) "estimate" else NULL)

  dat <- build_tmb_data(data, spatialCovs, model, coord, scaleFactor, smoothGradient, npoints, curweight, zetaScale, obs_dist, barrier_sdf = barrier_sdf, lambda = lambda, lambda_mode = lambda_mode, grad = grad, smoothPar = smoothPar, barrierOnset = barrierOnset)

  if (lambda_mode != 2L && !is.null(prior) && "log_lambda" %in% rownames(prior)) {
    stop("A prior on 'log_lambda' can only be specified when a barrier is specified and lambda = \"estimate\".")
  }
  origprior <- prior
  cp <- checkPar(par, model, map, dat=dat, spatialCovs=spatialCovs, prior=prior, lambda_mode=lambda_mode)
  tmbpar <- cp$par
  tmbmap <- cp$map
  re <- cp$re

  dat <- c(dat, cp$priors)
  tmbmap <- mapDuplicatedTimes(dat, tmbmap, tmbpar, re)

  # scale parameters internally for TMB; par$lambda (natural scale) -> log_lambda (working scale)
  tmbpar$log_sigma <- tmbpar$log_sigma - log(scaleFactor)
  if (lambda_mode == 2L) tmbpar$log_lambda <- tmbpar$log_lambda + 2 * log(scaleFactor)
  tmbpar$mu <- tmbpar$mu / scaleFactor
  tmbpar$vel <- tmbpar$vel / scaleFactor

  message("   Constructing objective function...")

  obj2 <- try({
    TMB::MakeADFun(c(model="langevinSSM",dat), tmbpar, map = tmbmap, random = re, DLL = "langevinSSM_TMBExports", hessian = hessian, method = method, silent = silent, inner.control = inner.control)
  }, silent = TRUE)

  if (inherits(obj2, "try-error")) stop("TMB::MakeADFun failed to construct the objective function. Check parameter initial values or data constraints.\nError details: ", attr(obj2, "condition")$message)

  message("   Fitting ",model," Langevin model...")

  if (length(obj2$par) == 0) {
    start <- proc.time()
    fit <- list(par = obj2$par, objective = obj2$fn(obj2$par), convergence = 0, message = "All parameters fixed; no outer optimization required", iterations = 0, evaluations = c("function" = 1, "gradient" = 0))
    fit$elapsedTime <- proc.time() - start
  } else {
    start <- proc.time()

    if (optMethod == "nlminb") {
      fit <- try({
        do.call(stats::nlminb,args = list(start = obj2$par, objective = obj2$fn, gradient = obj2$gr, control=control))
      }, silent = TRUE)
    } else if (optMethod == "nloptr") {

      # map standard BFGS to the nloptr equivalent
      algo <- method
      if (algo == "BFGS") algo <- "NLOPT_LD_LBFGS"

      # control list for nloptr mapping from the default nlminb list when possible
      opts <- list(algorithm = algo,
                   maxeval = if (!is.null(control$eval.max)) control$eval.max else 1000)

      # Overwrite with any user-provided controls mapping directly to nloptr (e.g., control$algorithm = "NLOPT_GD_MLSL")
      for (n in names(control)) {
        if (n %in% c("trace", "iter.max", "eval.max", "algorithm")) next # skip standard nlminb arguments and already set algorithms
        opts[[n]] <- control[[n]]
      }

      fit_try <- try({
        nloptr::nloptr(x0 = obj2$par, eval_f = function(x) obj2$fn(x), eval_grad_f = function(x) obj2$gr(x), opts = opts)
      }, silent = TRUE)

      if (!inherits(fit_try, "try-error")) {
        # Format nloptr output to match expected nlminb structure
        fit <- list(
          par = fit_try$solution,
          objective = fit_try$objective,
          convergence = if (fit_try$status < 0) fit_try$status else 0, # Map negative nloptr statuses to an error code
          message = fit_try$message,
          iterations = fit_try$iterations,
          evaluations = c("function" = fit_try$iterations, "gradient" = fit_try$iterations)
        )
        names(fit$par) <- names(obj2$par)
      } else {
        fit <- fit_try
      }
    }

    if (inherits(fit, "try-error") || !is.list(fit)) stop("Outer optimization failed. The model could not be fit.\nError details: ", attr(fit, "condition")$message)

    fit$elapsedTime <- proc.time() - start

    if (!is.null(fit$convergence) && fit$convergence != 0) warning("Optimization did not appear to converge. Code: ", fit$convergence, " - ", fit$message)

  }

  if (calcSE) {
    message("   Calculating SEs...")
    sdreport_out <- try({
      TMB::sdreport(obj2, getJointPrecision = getJointPrecision)
    }, silent = TRUE)
  } else {
    sdreport_out <- structure(list(), class = "skip_se")
  }

  extracted <- extract_tmb_estimates(fit, obj2, sdreport_out, re, tmbmap, data, scaleFactor, spatialCovs, getJointPrecision)
  fit$estimates <- extracted$estimates
  fit$covariance <- extracted$covariance

  if (hessian) {
    fit$covariance$hessian <- tryCatch(obj2$he(fit$par), error = function(e) NULL)
    if (!is.null(fit$covariance$hessian)) {
      rownames(fit$covariance$hessian) <- names(fit$par)
      colnames(fit$covariance$hessian) <- names(fit$par)
    }
  }

  out_of_bounds <- FALSE
  if (!is.null(fit$estimates$random$mu)) {
    mu_est <- fit$estimates$random$mu$est
    cov_ext <- as.vector(terra::ext(spatialCovs[[1]]))
    cov_res <- terra::res(spatialCovs[[1]])

    safe_xmin <- cov_ext[1] + cov_res[1]; safe_xmax <- cov_ext[2] - cov_res[1]
    safe_ymin <- cov_ext[3] + cov_res[2]; safe_ymax <- cov_ext[4] - cov_res[2]

    out_of_bounds <- any(mu_est$mu.x < safe_xmin | mu_est$mu.x > safe_xmax |
                           mu_est$mu.y < safe_ymin | mu_est$mu.y > safe_ymax, na.rm = TRUE)
  }

  lambda_report <- if (lambda_mode == 0L) lambda else as.numeric(obj2$report(obj2$env$last.par.best)$lambda)
  # "tied": lambda is a numerical device that grows without bound as the time step shrinks, so the
  # continuous-time model is a hard wall; downstream UD calculations use lambda = Inf (zero use of land).
  lambda_est <- if (lambda_mode == 1L) Inf else lambda_report
  if (is.null(barrier)) { lambda_est <- NULL; lambda_report <- NULL }
  if (lambda_mode == 2L) {
    # lambda is only informed by (smoothed) locations on the restricted side of the barrier
    n_restr <- tryCatch({
      n_obs <- sum(terra::extract(barrier_sdf, cbind(data$x, data$y))[, 1] < 0, na.rm = TRUE)
      n_mu <- if (!is.null(fit$estimates$random$mu)) sum(terra::extract(barrier_sdf, cbind(fit$estimates$random$mu$est$mu.x, fit$estimates$random$mu$est$mu.y))[, 1] < 0, na.rm = TRUE) else 0
      max(n_obs, n_mu)
    }, error = function(e) NA)
    lse <- tryCatch(fit$estimates$natural["lambda", "Std. Error"] / lambda_report, error = function(e) NA)
    if ((!is.na(n_restr) && n_restr < 5) || (calcSE && (is.na(lse) || !is.finite(lse) || lse > 2)))
      warning("lambda = \"estimate\": the barrier strength is poorly identified (", ifelse(is.na(n_restr), "?", n_restr),
              " observed/estimated locations on the restricted side of the barrier",
              if (calcSE) paste0("; relative SE ", ifelse(is.finite(lse), signif(lse, 2), "NaN")) else "",
              "). Consider fixing lambda or using lambda = \"tied\".", call. = FALSE)
  }
  fit$conditions <- list(hessian = hessian, method = method, optMethod = optMethod, silent = silent, inner.control = inner.control, control = control, scaleFactor = scaleFactor, model = model, obs_dist = obs_dist, par = par, map = map, grad = grad, smoothGradient = smoothGradient, npoints = npoints, curweight = curweight, zetaScale = zetaScale, smoothPar = smoothPar, prior = origprior, coord = coord, out_of_bounds = out_of_bounds, barrier = barrier, lambda = lambda_est, lambda_mode = if (is.null(barrier)) NULL else c("fixed","tied","estimate")[lambda_mode + 1L], lambda_report = lambda_report, barrierOnset = if (is.null(barrier)) NULL else barrierOnset, calcSE = calcSE)

  boundsWarning(fit)

  fit$signatures <- list(data = get_data_signature(data, coord), covs = get_covs_signature(orig_spatialCovs))
  fit$tmb_setup <- list(parList = obj2$env$parList(fit$par), map = tmbmap, random = re, priors = cp$priors)
  # Downstream functions rebuild the objective with lambda fixed at its estimate (lambda = "estimate"),
  # or with the same tied rule as the fit (lambda = "tied"; see .lambda_args()).
  fit$tmb_setup$map$log_lambda <- factor(NA)


  fit <- class_fitLangevin(fit)
  return(fit)
}

# Natural-scale (sigma, gamma) at which the smoothing radius is evaluated, derived from the data
# alone (never from 'par'), as zeta is in Blackwell & Matthiopoulos (2024): the measurement-error-
# corrected variogram estimates of initialValues(type = "variogram"), falling back to the moment
# estimators sigma = sqrt(mean(|dx|^2 / (2 dt))) (= B&M's zeta) and gamma = 1 / median(dt).
.smooth_par_data <- function(data, model, obs_dist, spatialCovs, coord) {
  iv <- suppressMessages(initialValues(data, model, obs_dist, spatialCovs = spatialCovs, coord = coord, type = "variogram"))
  sp <- list(sigma = as.numeric(iv$sigma), gamma = if (model == "underdamped") as.numeric(iv$gamma) else 1)
  if (!is.finite(sp$sigma) || sp$sigma <= 0 || !is.finite(sp$gamma) || sp$gamma <= 0)
    stop("Could not derive the gradient smoothing radius from the data.")
  sp
}

# (sigma, gamma) used for the smoothing radius of a fit, for downstream functions. Fits made
# before this was stored (whose radius tracked the current sigma and gamma) use the estimates.
.cond_smoothPar <- function(cond, fit = NULL) {
  if (!is.null(cond$smoothPar)) return(cond$smoothPar)
  pl <- fit$tmb_setup$parList
  if (!is.null(pl$log_sigma))
    return(list(sigma = exp(pl$log_sigma) * cond$scaleFactor,
                gamma = if (cond$model == "underdamped") exp(pl$log_gamma) else 1))
  list(sigma = 1, gamma = 1)
}
