#' Evaluate data and suggest barrier penalty (lambda)
#'
#' Evaluates raw tracking data against a barrier mask, calculates the mathematical floor and ceiling for the barrier penalty parameter (\code{lambda}), and prepares a cleaned dataset by converting statistically implausible locations into NA observations based on the reported measurement error.
#'
#' @param data A \code{dataLangevin} object containing the raw tracking data.
#' @param model Character string indicating which Langevin diffusion model to fit internally (``underdamped'', or ``overdamped'') to estimate baseline parameters. Default: ``underdamped''.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates.
#' @param barrier Character string specifying the name of the barrier mask within \code{spatialCovs}.
#' @param tolerance Numeric. The maximum acceptable distance (in \code{coord} units) that the estimated track is allowed to "leak" into the restricted zone. Default: \code{1}.
#' @param z_threshold Numeric. The maximum allowable statistical Z-score (leak / observation standard error) for a point inside the restricted zone. If a point exceeds this, it is considered overly confident about being in a restricted zone and is flagged as implausible. Default: \code{3}.
#' @param coord Character vector identifying the coordinate names. Default: \code{c("x", "y")}.
#' @param ... Additional arguments passed to \code{\link{fitLangevin}} during the internal baseline estimation (e.g., \code{map}, \code{control}, \code{scaleFactor}).
#'
#' @details
#' The suggested boundary penalty hyperparameter (\code{lambda}) is determined by balancing the physical requirement of the barrier against the numerical stability limits of the continuous-time movement model. The function executes a two-pass algorithm:
#'
#' \strong{1. Baseline Estimation & Pre-Filtering:}
#' To evaluate the model's stability, the animal's baseline diffusion rate (speed, \eqn{\sigma}) must be estimated. Before fitting this unconstrained baseline model, the function pre-filters any "statistically implausible" locations so they do not artificially inflate \eqn{\sigma}. A point is deemed implausible if its raw Z-score exceeds the \code{z_threshold}:
#' \deqn{Z = \frac{y}{\sigma_{obs}} > z_{threshold}}
#' where \eqn{y} is the distance leaked into the restricted zone and \eqn{\sigma_{obs}} is the marginal standard error of the observation pointing toward the barrier.
#'
#' \strong{2. The Stability Ceiling (\eqn{\lambda_{max}}):}
#' In continuous-time models, boundaries are enforced via a continuous penalty potential. When discretized for numerical optimization, if the penalty spring is too stiff relative to the animal's baseline diffusion rate (\eqn{\sigma}) and the maximum sampling interval (\eqn{\Delta t_{max}}), the model fitting algorithm can experience severe numerical instability. The mathematical ceilings for stability are:
#' \deqn{\lambda_{max} = \frac{2}{\sigma^2 \Delta t_{max}} \quad \text{(Overdamped)}}
#' \deqn{\lambda_{max} = \frac{1}{\sigma^2 \Delta t_{max}^2} \quad \text{(Underdamped)}}
#'
#' \strong{3. The Physical Floor (\eqn{\lambda_{min}}):}
#' The model estimates true latent locations (\eqn{\mu}) by balancing the data likelihood (pulling toward the observation) against the boundary penalty (pushing out of the barrier). To guarantee the estimated true track never leaks farther into the restricted zone than the user-defined \code{tolerance} (\eqn{\epsilon}), a minimum spring stiffness is required. For each observation, this is calculated as:
#' \deqn{\lambda_{req} = \frac{(y / \epsilon) - 1}{\sigma_{obs}^2}}
#' The global physical floor (\eqn{\lambda_{min}}) is defined as the maximum \eqn{\lambda_{req}} across all mathematically plausible points.
#'
#' \strong{4. Recommendation Logic:}
#' The recommended boundary penalty balances physical forces and numerical stability based on the signal-to-noise ratio (\eqn{R}) between the observation error and the movement process.
#' \itemize{
#'   \item \strong{Clean Data (\eqn{R < 0.5}):} The movement signal dominates the noise. The function recommends the stability ceiling (\eqn{\lambda_{max}}).
#'   \item \strong{Noisy Data (\eqn{R > 2.0}):} The noise dominates the movement. To prevent the model from artificially fighting GPS scatter, the function recommends the physical floor (\eqn{\lambda_{min}}).
#'   \item \strong{Balanced Data (\eqn{0.5 \le R \le 2.0}):} The function recommends the geometric mean of the floor and ceiling, placing the penalty exactly halfway between the two boundaries in log-space.
#'   \item \strong{No Restriction Conflict:} If all observed locations respect the barrier (or leaks are within the acceptable tolerance without requiring extra stiffness), the function defaults to a standard safe rigid wall defined as \eqn{\lambda_{max} / 10}.
#' }
#' If the required physical floor (\eqn{\lambda_{min}}) exceeds the stability ceiling (\eqn{\lambda_{max}}), the function caps the recommendation at the ceiling to prevent model fitting failure.
#'
#' \strong{Implausible Locations:}
#' Any location that fails the initial Z-score check, OR requires a point-specific stiffness (\eqn{\lambda_{req}}) greater than the stability ceiling (\eqn{\lambda_{max}}) is removed by setting its spatial coordinates and error parameters to \code{NA} in the returned \code{$filtered_data}. When fitting these filtered data, \code{\link{fitLangevin}} will treat these as missing observations and predict the latent locations based on the surrounding track, preserving the original observation times.
#'
#' @return A list containing:
#' \item{lambda_min}{The minimum stiffness required to satisfy the leak tolerance.}
#' \item{lambda_max}{The theoretical stability ceiling of the continuous-time movement model.}
#' \item{recommended}{The recommended \code{lambda} penalty.}
#' \item{Y_max}{The deepest observed distance into the restricted zone.}
#' \item{implausible_locations}{A data frame containing the row index, ID, date, leak distance into the restricted area, the specific error type ("KF" or "LS"), and the minimum required native error parameters (\code{req_val1} and \code{req_val2}, corresponding to \code{smaj}/\code{smin} or \code{x.err}/\code{y.err}) for the observation to be mathematically plausible. \code{NULL} if none exist.}
#' \item{filtered_data}{A \code{dataLangevin} object ready for final model fitting. If implausible locations were found, their spatial coordinates and error parameters are set to \code{NA}.}
#' @export
checkBarrier <- function(data, model = c("underdamped", "overdamped"), spatialCovs, barrier, tolerance = 1, z_threshold = 3, coord = c("x", "y"), ...) {

  if (!inherits(data, "dataLangevin")) stop("'data' must be a dataLangevin object.")
  model <- match.arg(model)
  .validate_barrier(barrier, spatialCovs)

  message("Evaluating tracking data against barrier geometry...")

  sdf_rast <- .get_barrier_sdf(barrier, spatialCovs)
  pts <- cbind(data[[coord[1]]], data[[coord[2]]])
  dist_vals <- suppressWarnings(terra::extract(sdf_rast, pts))[, 1]

  nogo_idx <- which(dist_vals <= 0)
  z_fail_idx <- integer(0)

  # --- PRE-PASS: Filter ONLY statistically implausible points for baseline fit ---
  if (length(nogo_idx) > 0) {
    for (idx in nogo_idx) {
      y_i <- abs(dist_vals[idx])
      err_sq_i <- NA

      if ("smaj" %in% names(data) && !is.na(data$smaj[idx])) {
        M2 <- (data$smaj[idx]^2) / 2
        m2 <- (data$smin[idx]^2) / 2
        err_sq_i <- max(M2, m2)
      } else if ("x.err" %in% names(data) && !is.na(data$x.err[idx])) {
        err_sq_i <- max(data$x.err[idx]^2, data$y.err[idx]^2)
      }

      if (is.na(err_sq_i) || err_sq_i == 0) {
        stop(sprintf("Observation at row %d is in restricted zone but lacks measurement error.", idx))
      }

      # Calculate raw Z-score to catch egregiously confident points
      if ((y_i / sqrt(err_sq_i)) > z_threshold) {
        z_fail_idx <- c(z_fail_idx, idx)
      }
    }
  }

  if (length(nogo_idx) > 0) {
    if (length(z_fail_idx) > 0) {
      message(sprintf("Found %d restricted locations. Pre-filtering %d statistically implausible points to estimate baseline speed (sigma)...", length(nogo_idx), length(z_fail_idx)))
      keep_idx <- setdiff(1:nrow(data), z_fail_idx)
      data_filter <- data[keep_idx, ]

      class(data_filter) <- class(data)
      attr(data_filter, "time.unit") <- attr(data, "time.unit")
      attr(data_filter, "lambda") <- attr(data, "lambda")
      attr(data_filter, "proj") <- attr(data, "proj")

      if (inherits(data_filter$date, "POSIXt") || inherits(data_filter$date, "Date")) {
        track_times <- as.numeric(difftime(data_filter$date, as.POSIXct("1970-01-01 00:00:00", tz = "UTC"), units = attr(data, "time.unit")))
      } else {
        track_times <- as.numeric(data_filter$date)
      }

      dt_calc <- c(0, diff(track_times))
      id_vec <- as.character(data_filter$id)
      dt_calc[id_vec != c(NA, id_vec[-nrow(data_filter)])] <- 0
      data_filter$dt <- dt_calc

    } else {
      message(sprintf("Found %d restricted locations (all statistically plausible). Estimating baseline speed (sigma)...", length(nogo_idx)))
      data_filter <- data
    }
  } else {
    message("No restricted locations found. Estimating baseline speed (sigma)...")
    data_filter <- data
  }

  # baseline fit
  fit_uncon <- tryCatch({
    suppressMessages(fitLangevin(data = data_filter, model = model, spatialCovs = spatialCovs, ...))
  }, error = function(e) e)

  if (inherits(fit_uncon, "error")) {
    stop("Internal baseline model fit failed. Ensure initial parameters and data are valid: ", fit_uncon$message)
  }

  get_est <- function(name, default) {
    if (name %in% rownames(fit_uncon$estimates$natural)) {
      val <- fit_uncon$estimates$natural[name, "Estimate"]
      if (!is.na(val)) return(val)
    }
    return(default)
  }

  sigma_est <- get_est("sigma", NA)
  if (is.na(sigma_est)) stop("Baseline model fit did not produce a valid 'sigma' estimate.")

  max_dt <- max(data$dt, na.rm = TRUE)
  if (max_dt <= 0) max_dt <- 1

  if (model == "overdamped") {
    lambda_max <- 2 / (sigma_est^2 * max_dt)
  } else {
    lambda_max <- 1 / (sigma_est^2 * max_dt^2)
  }

  lambda_min <- 0
  y_max <- 0
  sigma_err <- NA
  snr <- NA
  sigma_err_sq <- NA

  implausible_pts <- data.frame(
    row_index = integer(), id = character(), date = data$date[0],
    leak = numeric(), error_type = character(),
    req_val1 = numeric(), req_val2 = numeric(),
    stringsAsFactors = FALSE
  )

  message("Calculating theoretical model stability limits...")

  # --- RIGOROUS PASS: Final Stability Limit Check ---
  if (length(nogo_idx) > 0) {
    tau1_est <- get_est("tau_1", 1)
    tau2_est <- get_est("tau_2", 1)
    rho_o_est <- get_est("rho_o", 0)
    psi_est <- get_est("psi", 1)

    for (idx in nogo_idx) {
      y_i <- abs(dist_vals[idx])
      if (y_i > y_max) y_max <- y_i

      err_sq_i <- NA
      if ("smaj" %in% names(data) && !is.na(data$smaj[idx])) {
        M2 <- (data$smaj[idx]^2) / 2
        m2 <- ((data$smin[idx] * psi_est)^2) / 2
        err_sq_i <- max(M2, m2)
      } else if ("x.err" %in% names(data) && !is.na(data$x.err[idx])) {
        C11 <- (data$x.err[idx] * tau1_est)^2
        C22 <- (data$y.err[idx] * tau2_est)^2
        C12 <- data$x.err[idx] * data$y.err[idx] * tau1_est * tau2_est * rho_o_est
        trace_val <- C11 + C22
        det_val <- (C11 * C22) - (C12^2)
        err_sq_i <- (trace_val + sqrt(max(0, trace_val^2 - 4 * det_val))) / 2
      } else if (tau1_est != 1 || tau2_est != 1) {
        err_sq_i <- max(tau1_est^2, tau2_est^2)
      }

      sigma_obs_i <- sqrt(err_sq_i)
      z_obs_i <- y_i / sigma_obs_i

      # required stiffness for this specific point
      l_req_i <- ((y_i / tolerance) - 1) / err_sq_i
      if (l_req_i < 0) l_req_i <- 0

      is_implausible <- FALSE
      req_err_sd <- NA

      # Dual Check: Is it mathematically implausible?
      if (z_obs_i > z_threshold) {
        is_implausible <- TRUE
        req_err_sd <- y_i / z_threshold
      } else if (l_req_i > lambda_max) {
        is_implausible <- TRUE
        v_req <- max(0, ((y_i / tolerance) - 1) / lambda_max)
        req_err_sd <- sqrt(v_req)
      }

      if (is_implausible) {
        # Proportional scaling multiplier
        k_scale <- req_err_sd / sigma_obs_i

        if ("smaj" %in% names(data) && !is.na(data$smaj[idx])) {
          e_type <- "KF"
          r1 <- data$smaj[idx] * k_scale
          r2 <- data$smin[idx] * k_scale
        } else {
          e_type <- "LS"
          r1 <- data$x.err[idx] * k_scale
          r2 <- data$y.err[idx] * k_scale
        }

        implausible_pts <- rbind(implausible_pts, data.frame(
          row_index = idx,
          id = as.character(data$id[idx]),
          date = data$date[idx],
          leak = round(y_i, 3),
          error_type = e_type,
          req_val1 = round(r1, 3),
          req_val2 = round(r2, 3),
          stringsAsFactors = FALSE
        ))
      } else {
        if (l_req_i > lambda_min) {
          lambda_min <- l_req_i
          sigma_err_sq <- err_sq_i
        }
      }
    }

    if (!is.na(sigma_err_sq)) {
      sigma_err <- sqrt(sigma_err_sq)
      snr <- sigma_err / (sigma_est * sqrt(max_dt))
    }
  }

  filtered_data <- data

  cat("\n==================================================\n")
  cat("           BARRIER PENALTY DIAGNOSTICS\n")
  cat("==================================================\n")
  cat(sprintf("%-27s : %8.4f\n", "Baseline Speed (sigma)", sigma_est))
  cat(sprintf("%-27s : %8.4f\n", "Maximum Time Step (dt)", max_dt))
  cat(sprintf("%-27s : %8.4f units\n", "Max Leakage Observed", y_max))

  cat(sprintf("%-27s : %8.4f units\n", "Chosen Leak Tolerance", tolerance))
  if (lambda_min > 0) {
    cat(sprintf("%-27s : %8.4f\n", "Max Point-Specific Var", sigma_err_sq))
    cat(sprintf("%-27s : %8.4f\n", "Signal-to-Noise Ratio (R)", snr))
  }

  cat("\n--------------------------------------------------\n")
  cat("                 CONSTRAINTS\n")
  cat("--------------------------------------------------\n")
  cat(sprintf("%-27s : %8.4f\n", "Min Stiffness Required", lambda_min))
  cat(sprintf("%-27s : %8.4f\n", "Max Stability Limit", lambda_max))

  cat("\n")

  if (lambda_min > lambda_max) {
    recommended <- lambda_max
    cat("! WARNING: Required stiffness exceeds stability limit.\n")
    cat("! Defaulting to ceiling (lambda_max) to prevent model fitting failure.\n")
  } else if (lambda_min == 0) {
    recommended <- lambda_max / 10
    if (y_max == 0) {
      cat("* Observed data fully respect the restricted zones.\n")
      cat("* Recommending standard safe rigid wall (lambda_max / 10).\n")
    } else {
      cat("* All leaking locations were either within tolerance or flagged as implausible.\n")
      cat("* Recommending standard safe rigid wall (lambda_max / 10).\n")
    }
  } else {
    if (snr < 0.5) {
      recommended <- lambda_max
      cat("* Movement dominates noise (R < 0.5).\n")
      cat("* Recommending stability ceiling (lambda_max).\n")
    } else if (snr > 2.0) {
      recommended <- lambda_min
      cat("* Noise dominates movement (R > 2.0).\n")
      cat("* Recommending physical floor (lambda_min).\n")
    } else {
      safe_min <- max(lambda_min, 1e-8)
      recommended <- exp(mean(log(c(safe_min, lambda_max))))
      cat("* Balanced signal-to-noise (0.5 <= R <= 2.0).\n")
      cat("* Recommending geometric mean.\n")
    }
  }

  cat(sprintf("\n=> RECOMMENDED LAMBDA       : %8.4f\n", recommended))
  cat("==================================================\n")

  if (nrow(implausible_pts) > 0) {
    cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("! IMPLAUSIBLE LOCATIONS DETECTED\n")
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("The following observations reside so deeply in the restricted\n")
    cat("zone relative to their measurement error that they appear statistically\n")
    cat("implausible. To prevent numerical issues, their\n")
    cat("coordinates and errors have been set to NA in the returned\n")
    cat("$filtered_data data frame\n\n")

    print_limit <- min(5, nrow(implausible_pts))
    for (i in 1:print_limit) {
      if (implausible_pts$error_type[i] == "KF") {
        req_str <- sprintf("Needs smaj >= %.2f, smin >= %.2f",
                           implausible_pts$req_val1[i], implausible_pts$req_val2[i])
      } else {
        req_str <- sprintf("Needs x.err >= %.2f, y.err >= %.2f",
                           implausible_pts$req_val1[i], implausible_pts$req_val2[i])
      }
      cat(sprintf("  * Row %d (ID %s): %.2f units leak (%s)\n",
                  implausible_pts$row_index[i], implausible_pts$id[i],
                  implausible_pts$leak[i], req_str))
    }
    if (nrow(implausible_pts) > 5) {
      cat(sprintf("  ... and %d more. (See $implausible_locations)\n", nrow(implausible_pts) - 5))
    }
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

    bad_idx <- implausible_pts$row_index

    # Safely convert implausible observations to NA
    filtered_data[[coord[1]]][bad_idx] <- NA
    filtered_data[[coord[2]]][bad_idx] <- NA

    if ("smaj" %in% names(filtered_data)) {
      filtered_data$smaj[bad_idx] <- NA
      filtered_data$smin[bad_idx] <- NA
      filtered_data$eor[bad_idx] <- NA
    }

    if ("x.err" %in% names(filtered_data)) {
      filtered_data$x.err[bad_idx] <- NA
      filtered_data$y.err[bad_idx] <- NA
    }
  }

  return(invisible(list(
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    recommended = recommended,
    Y_max = y_max,
    implausible_locations = if (nrow(implausible_pts) > 0) implausible_pts else NULL,
    filtered_data = filtered_data
  )))
}
