#' Suggest the Optimal Boundary Penalty (lambda)
#'
#' Calculates the mathematical floor and ceiling for the barrier penalty parameter (\code{lambda})
#' based on observed measurement errors, theoretical SDE stability limits, and the kinetic energy of the track.
#'
#' @param data A \code{dataLangevin} object containing the tracking data.
#' @param fit A \code{fitLangevin} object fit WITHOUT a barrier constraint to estimate baseline parameters.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates.
#' @param barrier Character string specifying the name of the barrier mask within \code{spatialCovs}.
#' @param tolerance Numeric. The maximum acceptable distance (in \code{coord} units) that the estimated track is allowed to "leak" into the restricted zone. Default: \code{1}.
#' @param coord Character vector identifying the coordinate names. Default: \code{c("x", "y")}.
#'
#' @details
#' The optimal boundary penalty parameter (\code{lambda}) is determined by analyzing the signal-to-noise ratio between the observation error and the true movement process, operating within two strict mathematical boundaries:
#'
#' \strong{1. The Stability Ceiling (\code{lambda_max}):} If the penalty spring is too stiff relative to the animal's diffusion rate (\code{sigma}) and the maximum sampling interval (\code{dt}), the discrete-time SDE solver will experience severe numerical instability. The theoretical maximum stability limits are:
#' \itemize{
#'   \item Overdamped: \eqn{\lambda_{max} = \frac{2}{\sigma^2 \Delta t_{max}}}
#'   \item Underdamped: \eqn{\lambda_{max} = \frac{1}{\sigma^2 \Delta t_{max}^2}}
#' }
#'
#' \strong{2. The Physical Floor (\code{lambda_min}):} If the tracking data contain measurement error, observed locations may fall inside restricted zones. To guarantee the estimated true track (\eqn{\mu}) never leaks farther into restricted zones than the user-defined \code{tolerance} (\eqn{\epsilon}) when accounting for the worst-case observed distance within restricted zones (\eqn{Y_{max}}), the required minimum stiffness is:
#' \deqn{\lambda_{min} = \frac{\frac{Y_{max}}{\epsilon} - 1}{\sigma_{err}^2}}
#' where \eqn{\sigma_{err}^2} is the maximum eigenvalue of the point-specific measurement error covariance matrix for that deepest observation within the restricted zone (accounting for \code{psi}, \code{tau}, and \code{rho_o}).
#'
#' \strong{Recommendation Strategy (The Signal-to-Noise Ratio):}
#' The function calculates the ratio \eqn{R = \frac{\sigma_{err}}{\sigma \sqrt{\Delta t_{max}}}} to evaluate whether the data is dominated by the kinetic energy of the animal or the noise of the GPS.
#' \itemize{
#'   \item \strong{\eqn{R < 0.5} (Movement Dominates):} Recommends \code{lambda_max} to forcefully deflect the high-speed track.
#'   \item \strong{\eqn{R > 2.0} (Noise Dominates):} Recommends \code{lambda_min} to act as a statistical ``shock absorber''.
#'   \item \strong{\eqn{0.5 \le R \le 2.0} (Balanced):} Recommends the geometric mean of the floor and ceiling.
#' }
#' Note: If \code{lambda_min} exceeds \code{lambda_max} (often due to large time steps), the function defaults to \code{lambda_max} to help prevent numerical issues during model fitting.
#'
#' @return A list containing \code{lambda_min}, \code{lambda_max}, \code{y_max}, and the recommended \code{lambda}.
#' @export
suggestLambda <- function(data, fit, spatialCovs, barrier, tolerance = 1, coord = c("x", "y")) {

  if (!inherits(fit, "fitLangevin")) stop("'fit' must be a fitted langevinSSM object.")
  if (!is.null(fit$conditions$barrier)) stop("The provided 'fit' object must be run WITHOUT a barrier constraint to estimate baseline parameters.")

  .validate_barrier(barrier, spatialCovs)

  get_est <- function(name, default) {
    if (name %in% rownames(fit$estimates$natural)) {
      val <- fit$estimates$natural[name, "Estimate"]
      if (!is.na(val)) return(val)
    }
    return(default)
  }

  sigma_est <- get_est("sigma", NA)
  if (is.na(sigma_est)) stop("Model fit does not contain a valid 'sigma' estimate.")

  max_dt <- max(data$dt, na.rm = TRUE)
  if (max_dt <= 0) max_dt <- 1

  model_type <- fit$conditions$model
  if (model_type == "overdamped") {
    lambda_max <- 2 / (sigma_est^2 * max_dt)
  } else {
    lambda_max <- 1 / (sigma_est^2 * max_dt^2)
  }

  sdf_rast <- .get_barrier_sdf(barrier, spatialCovs)
  pts <- cbind(data[[coord[1]]], data[[coord[2]]])

  dist_vals <- suppressWarnings(terra::extract(sdf_rast, pts))[, 1]

  nogo_idx <- which(dist_vals <= 0)

  if (length(nogo_idx) == 0) {
    message("No points observed inside the restricted barrier zone.")
    lambda_min <- 0
    y_max <- 0
    sigma_err <- NA
    snr <- NA
    sigma_err_sq <- NA
  } else {
    # Find the deepest nogo point
    nogo_dists <- abs(dist_vals[nogo_idx])
    deepest_local_idx <- which.max(nogo_dists)
    global_idx <- nogo_idx[deepest_local_idx]
    y_max <- nogo_dists[deepest_local_idx]

    tau1_est <- get_est("tau_1", 1)
    tau2_est <- get_est("tau_2", 1)
    rho_o_est <- get_est("rho_o", 0)
    psi_est <- get_est("psi", 1)

    # calculate the maximum variance (largest eigenvalue) pulling the point into restricted zone
    sigma_err_sq <- NA

    if ("smaj" %in% names(data) && !is.na(data$smaj[global_idx])) {
      M2 <- (data$smaj[global_idx]^2) / 2
      m2 <- ((data$smin[global_idx] * psi_est)^2) / 2
      sigma_err_sq <- max(M2, m2)

    } else if ("x.err" %in% names(data) && !is.na(data$x.err[global_idx])) {
      C11 <- (data$x.err[global_idx] * tau1_est)^2
      C22 <- (data$y.err[global_idx] * tau2_est)^2
      C12 <- data$x.err[global_idx] * data$y.err[global_idx] * tau1_est * tau2_est * rho_o_est

      trace_val <- C11 + C22
      det_val <- (C11 * C22) - (C12^2)
      sigma_err_sq <- (trace_val + sqrt(max(0, trace_val^2 - 4 * det_val))) / 2

    } else if (tau1_est != 1 || tau2_est != 1) {
      sigma_err_sq <- max(tau1_est^2, tau2_est^2)
    }

    if (is.na(sigma_err_sq) || sigma_err_sq == 0) {
      stop("The farthest observation into the restricted zone has no point-specific measurement error.\nIf true locations are inside restricted zones, the barrier is invalid.")
    }

    # signal-to-noise ratio
    sigma_err <- sqrt(sigma_err_sq)
    snr <- sigma_err / (sigma_est * sqrt(max_dt))

    lambda_min <- ((y_max / tolerance) - 1) / sigma_err_sq
    if (lambda_min < 0) lambda_min <- 0
  }

  cat("Boundary Penalty Diagnostics:\n")
  cat("-----------------------------\n")
  cat(sprintf("Baseline speed (sigma): %.3f\n", sigma_est))
  cat(sprintf("Maximum time step (dt):     %.3f\n", max_dt))
  cat(sprintf("Max observed leak into %s (Y_max):  %.3f units\n", barrier, y_max))

  if (lambda_min > 0 || y_max > 0) {
    cat(sprintf("Chosen leak tolerance:      %.3f units\n", tolerance))
    cat(sprintf("Point-specific error (VAR): %.3f\n", sigma_err_sq))
    cat(sprintf("Signal-to-noise ratio (R):  %.3f\n", snr))
  }

  cat("\nConstraints:\n")
  cat("-----------------------------\n")
  cat(sprintf("Minimum stiffness required: %.4f\n", lambda_min))
  cat(sprintf("Maximum stability limit:    %.4f\n", lambda_max))

  if (lambda_min > lambda_max) {
    recommended <- lambda_max
    cat("WARNING: Minimum stiffness exceeds theoretical stability limit.\n")
    cat("Defaulting to ceiling to prevent SDE solver failure.\n")
  } else if (lambda_min == 0) {
    recommended <- lambda_max / 10
    if (y_max == 0) {
      cat("Observed data respect the barrier. Recommending safely rigid value (lambda_max / 10).\n")
    } else {
      cat("Maximum leak does not exceed requested tolerance. Recommending safely rigid value (lambda_max / 10).\n")
    }
  } else {
    if (snr < 0.5) {
      recommended <- lambda_max
      cat("Movement dominates noise (R < 0.5). Recommending stability ceiling (lambda_max).\n")
    } else if (snr > 2.0) {
      recommended <- lambda_min
      cat("Noise dominates movement (R > 2.0). Recommending physical floor (lambda_min).\n")
    } else {
      # Use 1e-8 as a safe log floor just in case lambda_min evaluates exactly to 0
      safe_min <- max(lambda_min, 1e-8)
      recommended <- exp(mean(log(c(safe_min, lambda_max))))
      cat("Balanced signal-to-noise (0.5 <= R <= 2.0). Recommending geometric mean.\n")
    }
  }

  cat(sprintf("\n--> Recommended Lambda: %.4f\n", recommended))

  return(invisible(list(
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    recommended = recommended,
    Y_max = y_max
  )))
}
