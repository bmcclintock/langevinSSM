#' Iteratively fit the Langevin model with a barrier penalty
#'
#' Fits a Langevin diffusion model constrained by a spatial barrier using a two-stage grid search.
#' It evaluates penalties by using posterior predictive checks (PPCs) to calculate the Kolmogorov-Smirnov (KS)
#' distance between the observed and simulated distances from the barrier. It uses a linearly spaced top-down coarse
#' search grid to identify the most optimal penalty (\code{lambda}), followed by a finer linear search around the winner from the coarse grid search.
#'
#' @param data A \code{dataLangevin} object containing the tracking data.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates.
#' @param lambda_max Numeric. The maximum barrier penalty to start the top-down search. If \code{NULL}, it is estimated from the empirical step variance.
#' @param n_coarse Integer. The number of linearly spaced penalties to evaluate in the coarse grid. Default: \code{3}.
#' @param n_fine Integer. The number of linearly spaced penalties to evaluate in the fine grid around the coarse winner. Default: \code{4}.
#' @param n_sims Integer. The number of repeated simulations to run per penalty to calculate the overall (median) KS score. Default: \code{5}.
#' @param timeStep Time step to use for the simulation(s). Determines the resolution of the discrete-time approximation of the continuous-time process. The smaller the \code{timeStep}, the more accurate the approximation. Default: 0.01.
#' @param ... Additional arguments passed to \code{\link{fitLangevin}}.
#'
#' @template barrier_details
#' @details
#' \strong{Why a Grid Search?}
#' The barrier penalty (\code{lambda}) cannot be estimated at the same time as the other movement parameters. To find the most appropriate penalty, this function uses a two-phase search to try and balance boundary enforcement with numerical stability.
#'
#' \strong{Two-phase search strategy:}
#' \itemize{
#'   \item \strong{Phase I (Coarse Search):} If \code{lambda_max} is not provided, the algorithm calculates the maximum penalty the numerical solver can process. It then evaluates a set of \code{n_coarse} penalties, spaced evenly from \code{lambda_max} down to a lower boundary.
#'   \item \strong{Phase II (Fine Search):} After finding the most effective penalty from Phase I, the algorithm creates a narrower search window around that value. It then evaluates \code{n_fine} evenly spaced penalties to refine the estimate.
#' }
#'
#' \strong{Calculating the Search Boundaries}
#' If \code{lambda_max} is \code{NULL}, the algorithm calculates the upper and lower boundaries for the grid search using empirical properties of the tracking data:
#' \itemize{
#'   \item \strong{Upper Boundary (\code{lambda_max}):} This is estimated using the maximum observed time step between locations, combined with the animal's empirical speed (and friction, if using the underdamped model). Higher values beyond this limit will likely lead to numerical issues during model fitting.
#'   \item \strong{Lower Boundary (\code{lambda_min}):} This represents the weakest plausible penalty. It is calculated by finding the maximum distance the animal was observed inside the restricted area. The algorithm sets \code{lambda_min} just high enough to ensure the barrier's restoring force could counteract the animal's empirical kinetic energy at that specific depth.
#' }
#'
#' \strong{The Scoring Metric: Comparing Distance Distributions}
#' To determine which \code{lambda} is best, the function compares simulated tracks to the real data. An ideal penalty will generate simulated tracks that interact with the boundary in the same way the real observed tracks do.
#'
#' For each candidate penalty in the search grid, the function:
#' \enumerate{
#'   \item Fits the Langevin model.
#'   \item Simulates \code{n_sims} new tracks from the fitted model.
#'   \item Calculates the distance from the boundary for every simulated location.
#'   \item Compares the distribution of these simulated distances to the actual distances observed in the tracking data.
#' }
#'
#' The similarity between the simulated and real distances is measured using the two-sample Kolmogorov-Smirnov (KS) test. The KS test provides a score between 0 and 1, where lower scores indicate a closer match to the real data. Because simulations involve randomness, the algorithm runs \code{n_sims} simulations per penalty and uses the \strong{median KS score}.
#'
#' @return A \code{fitLangevin} object that includes the most optimal penalty from the grid search.
#' @export
fitLangevin_barrier <- function(data, spatialCovs, lambda_max = NULL, n_coarse = 3, n_fine = 4, n_sims = 5, timeStep = 0.01, ...) {

  if (!inherits(data, "dataLangevin")) stop("'data' must be a dataLangevin object.")
  barrier <- .find_barrier(spatialCovs)
  if (is.null(barrier)) stop("No barrier found in 'spatialCovs'. Did you run prepBarrier()?")
  sdf_rast <- spatialCovs[[barrier]]

  args <- list(...)
  model_type <- if (!is.null(args$model)) match.arg(args$model, c("underdamped", "overdamped")) else "underdamped"
  coord <- if (!is.null(args$coord)) args$coord else c("x", "y")

  checkErrorData(data, coord)
  max_dt <- max(data$dt, na.rm = TRUE)
  if (max_dt <= 0) max_dt <- 1

  # --- Pre-calculate Observed Distances for PPC Scoring ---
  obs_pts <- cbind(data[[coord[1]]], data[[coord[2]]])

  # Extract raster metadata once for the fast Rcpp extractor
  r_ext <- as.numeric(as.vector(terra::ext(sdf_rast)))
  r_res <- as.numeric(as.vector(terra::res(sdf_rast)))
  b_mat <- terra::as.matrix(sdf_rast, wide = TRUE)

  obs_dist <- extract_sdf_rcpp(obs_pts, b_mat, r_ext, r_res)
  obs_dist <- obs_dist[!is.na(obs_dist)]

  if (length(obs_dist) == 0) stop("No valid spatial observations to calculate distance.")

  # --- 1. Top-Down Grid Initialization ---
  message("   Initializing predictive top-down search grids...")

  dx <- diff(data[[coord[1]]])
  dy <- diff(data[[coord[2]]])
  dt_steps <- data$dt[-1]
  valid <- which(dt_steps > 0 & !is.na(dx) & !is.na(dy))

  if(length(valid) < 3) {
    stop("Insufficient valid steps to auto-calculate empirical parameters. The dataset must contain at least 4 valid, consecutive locations.")
  }

  mean_dt <- mean(dt_steps[valid])

  if (model_type == "underdamped") {
    # 1. Empirical Gamma via velocity autocorrelation (OU mean-reversion)
    vx <- dx[valid] / dt_steps[valid]
    vy <- dy[valid] / dt_steps[valid]

    rho_x <- suppressWarnings(stats::cor(vx[-length(vx)], vx[-1], use = "complete.obs"))
    rho_y <- suppressWarnings(stats::cor(vy[-length(vy)], vy[-1], use = "complete.obs"))
    rho_mean <- mean(c(rho_x, rho_y), na.rm = TRUE)

    if (is.na(rho_mean)) rho_mean <- 0.01

    # Bound rho between high friction (0.001) and highly ballistic (0.999)
    rho_bounded <- max(0.001, min(0.999, rho_mean))
    empirical_gamma <- -log(rho_bounded) / mean_dt

    # 2. Empirical Sigma (Underdamped formulation)
    msd <- mean(dx[valid]^2 + dy[valid]^2)
    num_sig <- msd * (empirical_gamma^2)
    den_sig <- 2 * (empirical_gamma * mean_dt - 1 + exp(-empirical_gamma * mean_dt))
    empirical_sigma <- sqrt(num_sig / den_sig)

  } else {
    # Empirical Sigma (Overdamped formulation)
    empirical_sigma <- sqrt(mean((dx[valid]^2 + dy[valid]^2) / (2 * mean_dt)))
  }

  # Set lambda_max
  if (is.null(lambda_max)) {
    if (model_type == "overdamped") {
      lambda_max <- 2 / (empirical_sigma^2 * max_dt)
    } else {
      # Use exact analytical stability limit for the underdamped SDE solver
      num <- (empirical_gamma^2) * (1 - exp(-empirical_gamma * max_dt))
      den <- (empirical_sigma^2) * (1 - exp(-empirical_gamma * max_dt) - (empirical_gamma * max_dt * exp(-empirical_gamma * max_dt)))
      lambda_max <- num / den
    }

    #lambda_max <- lambda_max * 0.99
    message("     'lambda_max' was not provided. Using data-driven empirical stability limit: ",round(lambda_max,4),".\n     Note that measurement error can break autocorrelation and cause the empirical limit to be inflated.\n")
  }

  # Set lambda_min
  valid_dists <- obs_dist[obs_dist > 0]
  D_max <- if (length(valid_dists) > 0) max(valid_dists, na.rm = TRUE) else max(diff(range(data[[coord[1]]], na.rm=TRUE)), diff(range(data[[coord[2]]], na.rm=TRUE)))
  lambda_min <- max(1e-6, (empirical_sigma^2) / (2 * D_max^2))

  if (lambda_min >= lambda_max) lambda_min <- lambda_max * 0.01

  # --- Construct Linear Coarse Grid ---
  coarse_grid <- seq(lambda_max, lambda_min, length.out = n_coarse)

  message(sprintf("    Upper boundary: %.4f", lambda_max))
  message(sprintf("    Lower boundary: %.6f", lambda_min))
  message(sprintf("    Linear coarse grid: [%.4f down to %.6f] with %d steps\n", lambda_max, lambda_min, n_coarse))

  # --- Helper: Model Fitting & Scoring Function ---
  score_lambda <- function(lam) {
    message(sprintf("  -> Testing lambda = %.4f...", lam))

    fit_args <- args
    fit_args$data <- data
    fit_args$spatialCovs <- spatialCovs
    fit_args$lambda <- lam
    fit_args$silent <- TRUE

    fit <- tryCatch(do.call(fitLangevin, fit_args), error = function(e) e)

    if (inherits(fit, "error") || (!is.null(fit$convergence) && fit$convergence != 0)) {
      err_msg <- if(inherits(fit, "error")) fit$message else paste("Convergence code", fit$convergence)
      message(sprintf("     [Fit failed: %s]", err_msg))
      return(list(score = Inf, fit = NULL))
    }

    sim_args <- list(model = fit, data = data, spatialCovs = spatialCovs, timeStep = timeStep, conditional = FALSE)

    ks_scores <- c()
    err_msgs <- c()

    for (s in 1:n_sims) {
      if(s==1) message("   Simulating tracks forward using movement parameters fixed at point estimates...")
      sim_data <- tryCatch(suppressMessages(do.call(simLangevin, sim_args)), error = function(e) e)

      if (inherits(sim_data, "error")) {
        err_msgs <- c(err_msgs, sim_data$message)
        next
      }

      sim_pts <- cbind(sim_data[[coord[1]]], sim_data[[coord[2]]])
      sim_dist <- suppressWarnings(terra::extract(sdf_rast, sim_pts, method = "bilinear"))[, 1]
      sim_dist <- sim_dist[!is.na(sim_dist)]

      if(length(sim_dist) == 0) {
        err_msgs <- c(err_msgs, "No valid spatial observations extracted")
        next
      }

      ks_scores <- c(ks_scores, suppressWarnings(stats::ks.test(obs_dist, sim_dist)$statistic))
    }

    if (length(ks_scores) == 0) {
      message(sprintf("     [All simulations failed. First error: %s]", err_msgs[1]))
      return(list(score = Inf, fit = fit))
    }

    overall_ks <- stats::median(ks_scores)
    ks_details <- paste(sprintf("%.4f", ks_scores), collapse = ", ")
    message(sprintf("     [Median KS Score: %.4f (over %d valid sims) | Individual KS: %s]", overall_ks, length(ks_scores), ks_details))
    return(list(score = unname(overall_ks), fit = fit))
  }

  # --- 2. Phase I: Top-Down Coarse Linear Search ---
  message("==================================================")
  message(" PHASE I: Top-Down Coarse Linear Search")
  message("==================================================")

  coarse_results <- list()
  for (i in seq_along(coarse_grid)) {
    coarse_results[[i]] <- score_lambda(coarse_grid[i])
  }

  coarse_scores <- sapply(coarse_results, function(x) x$score)
  valid_coarse_idx <- which(!is.infinite(coarse_scores))

  if (length(valid_coarse_idx) == 0) {
    stop("All models in the coarse grid failed. Check your data, initialization parameters, or provide a lower lambda_max.")
  }

  best_coarse_idx <- valid_coarse_idx[which.min(coarse_scores[valid_coarse_idx])]
  best_coarse_lam <- coarse_grid[best_coarse_idx]
  best_coarse_score <- coarse_scores[best_coarse_idx]

  message(sprintf("\n=> Coarse Winner: lambda = %.4f (Median KS Score: %.4f)\n", best_coarse_lam, best_coarse_score))

  # --- 3. Phase II: Fine Linear Search ---
  message("==================================================")
  message(" PHASE II: Fine Linear Search")
  message("==================================================")

  # Determine bounds for the fine grid centrally around the winner
  step_size <- (lambda_max - lambda_min) / (n_coarse - 1)
  upper_bound <- best_coarse_lam + step_size
  lower_bound <- max(lambda_min, best_coarse_lam - step_size)

  # Create linear fine grid (descending)
  fine_grid <- seq(upper_bound, lower_bound, length.out = n_fine + 2)[-c(1, n_fine + 2)]

  # Exclude the coarse winner to prevent redundant fitting and simulation noise
  fine_grid <- fine_grid[abs(fine_grid - best_coarse_lam) > 1e-8]

  fine_results <- list()
  for (i in seq_along(fine_grid)) {
    fine_results[[i]] <- score_lambda(fine_grid[i])
  }

  fine_scores <- sapply(fine_results, function(x) x$score)
  valid_fine_idx <- which(!is.infinite(fine_scores))

  # Combine coarse winner and valid fine results
  all_lams <- best_coarse_lam
  all_scores <- best_coarse_score
  all_fits <- list(coarse_results[[best_coarse_idx]]$fit)

  if (length(valid_fine_idx) > 0) {
    all_lams <- c(all_lams, fine_grid[valid_fine_idx])
    all_scores <- c(all_scores, fine_scores[valid_fine_idx])
    all_fits <- c(all_fits, lapply(fine_results[valid_fine_idx], function(x) x$fit))
  }

  best_final_idx <- which.min(all_scores)
  best_final_lam <- all_lams[best_final_idx]
  best_final_score <- all_scores[best_final_idx]
  best_final_fit <- all_fits[[best_final_idx]]

  message(sprintf("\n=================================================="))
  message(sprintf(" OPTIMAL LAMBDA FOUND: %.4f (Median KS Score: %.4f)", best_final_lam, best_final_score))
  message(sprintf("=================================================="))

  return(best_final_fit)

}
