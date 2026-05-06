#' Iteratively fit the Langevin model with a barrier penalty
#'
#' Fits a Langevin diffusion model constrained by a spatial barrier using a two-stage grid search.
#' It evaluates penalties by using posterior predictive checks (PPCs) to calculate the Kolmogorov-Smirnov (KS)
#' distance between the observed and simulated distances from the barrier. It uses a linearly spaced top-down coarse
#' search grid to identify the most optimal penalty (\code{lambda}), followed by a finer linear search around the winner from the coarse grid search.
#'
#' @param data A \code{dataLangevin} object containing the tracking data.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates.
#' @param barrier Character string. The name of the barrier in \code{spatialCovs} that is represented as a signed distance field (see \code{\link{prepBarrier}}). If provided, this raster is exclusively used for the barrier penalty and is not included in the habitat selection covariates. Default: \code{NULL} (no barrier).
#' @param lambda_max Numeric. The maximum barrier penalty to start the top-down search. If \code{NULL}, it is estimated automatically using a "baseline 0" model fit that assumes \code{lambda}=0 (i.e., no barrier).
#' @param n_coarse Integer. The number of linearly spaced penalties to evaluate in the coarse grid. Default: \code{3}.
#' @param n_fine Integer. The number of linearly spaced penalties to evaluate in the fine grid around the coarse winner. Default: \code{4}.
#' @param n_sims Integer. The number of repeated simulations to run per penalty to calculate the overall (median) KS score. Default: \code{5}.
#' @param timeStep Time step to use for the simulation(s). Determines the resolution of the discrete-time approximation of the continuous-time process. The smaller the \code{timeStep}, the more accurate the approximation. Default: 0.01.
#' @param ncores Integer. Number of cores to use for parallel processing of the search grids. Default: \code{1} (sequential).
#' @param ... Additional arguments passed to \code{\link{fitLangevin}}.
#'
#' @template barrier_details
#' @details
#' \strong{Why a Grid Search?}
#' The barrier penalty (\code{lambda}) cannot be estimated at the same time as the other movement parameters. To find the most appropriate penalty, this function uses a two-phase search to try and balance boundary enforcement with numerical stability.
#'
#' \strong{Two-phase search strategy:}
#' \itemize{
#'   \item \strong{Phase I (Coarse Search):} It first evaluates a set of \code{n_coarse} penalties, spaced evenly from \code{lambda_max} down to a lower boundary.
#'   \item \strong{Phase II (Fine Search):} After finding the ``best'' penalty from Phase I, the algorithm creates a narrower search window around that value. It then evaluates \code{n_fine} evenly spaced penalties to refine the penalty.
#' }
#'
#' \strong{Calculating the Search Boundaries}
#' If \code{lambda_max} is \code{NULL}, the algorithm calculates the upper and lower boundaries for the grid search using a "baseline 0" approach:
#' \itemize{
#'   \item \strong{Upper Boundary (\code{lambda_max}):} The algorithm first fits a baseline model with the barrier penalty temporarily turned off (i.e., \code{lambda = 0}). This allows the model to approximate the movement speed and friction. The algorithm then uses these rough estimates to calculate the maximum barrier penalty the model appears to be able to handle.
#'   Note that this is just a ballpark estimate based on a model that ignores the barrier and can be corrupted by measurement error. It may therefore not be optimal, and users are encouraged to explore other (larger) values for \code{lambda_max}.
#'   \item \strong{Lower Boundary (\code{lambda_min}):} This represents the weakest plausible penalty. It is calculated by finding the maximum distance an animal was observed inside the restricted area. The algorithm sets \code{lambda_min} just high enough to ensure the barrier's restoring force could push the animal back from that specific depth.
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
#' @return A \code{fitLangevin} object that includes the ``best'' penalty from the grid search.
#' @importFrom terra wrap unwrap
#' @export
tuneBarrier <- function(data, spatialCovs, barrier, lambda_max = NULL, n_coarse = 3, n_fine = 4, n_sims = 5, timeStep = 0.01, ncores = 1, ...) {

  if (!inherits(data, "dataLangevin")) stop("'data' must be a dataLangevin object.")
  if (missing(spatialCovs) || !is.list(spatialCovs)) stop("'spatialCovs' list must be provided.")

  if (missing(barrier) || is.null(barrier)) {
    sim_barrier <- attr(data,"barrier")
    if(!is.null(sim_barrier)){
      barrier <- sim_barrier
    } else stop("'barrier' argument (the name of the barrier raster in 'spatialCovs') must be provided.")
  }
  if (!(barrier %in% names(spatialCovs))) stop(sprintf("Barrier raster '%s' not found in spatialCovs.", barrier))

  args <- list(...)
  model_type <- if (!is.null(args$model)) match.arg(args$model, c("underdamped", "overdamped")) else "underdamped"
  coord <- if (!is.null(args$coord)) args$coord else c("x", "y")
  scaleFactor <- if (!is.null(args$scaleFactor)) args$scaleFactor else 1

  sdf_rast <- spatialCovs[[barrier]]
  if(!isTRUE(attr(sdf_rast,"barLangevin"))) stop("barrier is not a 'barLangevin' object created by prepBarrier.")
  checkErrorData(data, coord)

  global_warnings <- character()

  # --- POSIXt Type Enforcement & Parsing ---
  time.unit <- attr(data, "time.unit")
  if (inherits(data$date, c("POSIXt", "Date"))) {
    if (!is.character(timeStep)) {
      stop("When the 'date' column is POSIXt or Date, 'timeStep' must be a character string specifying the time interval (e.g., '1 sec', '30 mins', '1 hour', '6 hours').")
    }
    if (is.null(time.unit)) time.unit <- "secs"

    t0 <- as.POSIXct("1970-01-01", tz = "UTC")
    t1 <- tryCatch(seq(t0, by = timeStep, length.out = 2)[2], error = function(e) NA)

    if (is.na(t1)) stop("Invalid 'timeStep' string provided. Use valid base R formats like '1 sec', '30 mins', or '1 hour'.")

    timeStep_num <- as.numeric(difftime(t1, t0, units = time.unit))
    if (timeStep_num <= 0) stop("Invalid 'timeStep' string provided. Must result in a positive duration.")
  } else {
    if (!is.numeric(timeStep)) stop("When the 'date' column is numeric, 'timeStep' must also be numeric.")
  }

  max_dt <- max(data$dt, na.rm = TRUE)
  if (max_dt <= 0) max_dt <- 1

  # --- Warning Deduplication Architecture ---
  seen_warnings <- new.env(parent = emptyenv())
  dedup_warning_handler <- function(w) {
    msg <- conditionMessage(w)
    if (grepl("cannot compute exact p-value with ties", msg)) {
      invokeRestart("muffleWarning")
    }
    if (is.null(seen_warnings[[msg]])) {
      seen_warnings[[msg]] <- TRUE
      warning(msg, call. = FALSE, immediate. = TRUE)
    }
    invokeRestart("muffleWarning")
  }

  # --- Setup Parallel Processing ---
  if (ncores > 1) {
    if (!requireNamespace("foreach", quietly = TRUE) || !requireNamespace("doFuture", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)  || !requireNamespace("doRNG", quietly = TRUE)) {
      stop("Packages 'foreach', 'future', 'doFuture', and 'doRNG' are required for multicore processing. Please install them.")
    } else {
      oldDoPar <- doFuture::registerDoFuture()
      on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)
      future::plan(future::multisession, workers = ncores)
      `%loop%` <- doRNG::`%dorng%`
    }
  } else {
    if (!requireNamespace("foreach", quietly = TRUE)) stop("Package 'foreach' is required. Please install it.")
    `%loop%` <- foreach::`%do%`
  }

  obs_pts <- cbind(data[[coord[1]]], data[[coord[2]]])
  obs_dist <- terra::extract(sdf_rast, obs_pts)[, 1]
  obs_dist <- obs_dist[!is.na(obs_dist)]

  if (length(obs_dist) == 0) stop("No valid spatial observations to calculate distance.")

  # --- 1. Top-Down Grid Initialization (Phase 0) ---
  message("   Initializing Top-Down Search Grids...")

  if (is.null(lambda_max)) {
    message("     Method: Baseline 0 fit (lambda = 0)")

    fit0_args <- args
    fit0_args$data <- data
    fit0_args$spatialCovs <- spatialCovs
    fit0_args$barrier <- barrier
    fit0_args$lambda <- 0
    fit0_args$silent <- TRUE

    withCallingHandlers({
      fit0 <- suppressMessages(tryCatch(do.call(fitLangevin, fit0_args), error = function(e) e))
    }, warning = dedup_warning_handler)

    if (inherits(fit0, "error")) {
      stop("Baseline 0 fit failed. Please provide 'lambda_max' manually.")
    }

    est_sigma <- fit0$estimates$natural["sigma", "Estimate"]
    est_sigma_work <- est_sigma / scaleFactor

    valid_dists <- obs_dist[obs_dist > 0]
    D_max <- if (length(valid_dists) > 0) max(valid_dists, na.rm = TRUE) else max(diff(range(data[[coord[1]]], na.rm=TRUE)), diff(range(data[[coord[2]]], na.rm=TRUE))) / scaleFactor

    # heuristic formula for lambda_min is invariant and does not scale with 1/L^2 like lambda_max does
    # normalize the spatial parameters, calculate the boundaries, and then apply the 1/L^2 transformation
    RS <- 10^(floor(log10(est_sigma_work)))

    sigma_rs <- est_sigma_work / RS
    D_max_rs <- D_max / RS

    if (model_type == "underdamped") {
      est_gamma <- fit0$estimates$natural["gamma", "Estimate"]

      num <- (est_gamma^2) * (1 - exp(-est_gamma * max_dt))
      den <- (sigma_rs^2) * (1 - exp(-est_gamma * max_dt) - (est_gamma * max_dt * exp(-est_gamma * max_dt)))
      lambda_max_rs <- num / den
    } else {
      lambda_max_rs <- 2 / (sigma_rs^2 * max_dt)
    }

    lambda_min_rs <- max(1e-12, (sigma_rs^2) / (2 * D_max_rs^2))

    if (lambda_min_rs >= lambda_max_rs) lambda_min_rs <- lambda_max_rs * 0.01

    lambda_max <- (lambda_max_rs / (RS^2)) * 0.99
    lambda_min <- lambda_min_rs / (RS^2)

    message(sprintf("     Calculated Maximum Penalty: %g\n", lambda_max))
  } else {
    message("     Method: User-provided limit\n")
    var_x <- stats::var(data[[coord[1]]], na.rm = TRUE)
    var_y <- stats::var(data[[coord[2]]], na.rm = TRUE)
    time_num <- as.numeric(data$date)
    T_total <- max(time_num, na.rm = TRUE) - min(time_num, na.rm = TRUE)
    D_macro <- (var_x + var_y) / max(1, T_total)

    est_sigma <- if (model_type == "underdamped") sqrt(D_macro) else sqrt(2 * D_macro)
    est_sigma_work <- est_sigma / scaleFactor

    valid_dists <- obs_dist[obs_dist > 0]
    D_max <- if (length(valid_dists) > 0) max(valid_dists, na.rm = TRUE) else max(diff(range(data[[coord[1]]], na.rm=TRUE)), diff(range(data[[coord[2]]], na.rm=TRUE))) / scaleFactor

    RS <- 10^(floor(log10(est_sigma_work)))

    sigma_rs <- est_sigma_work / RS
    D_max_rs <- D_max / RS

    lambda_max_rs <- lambda_max * (RS^2)
    lambda_min_rs <- max(1e-12, (sigma_rs^2) / (2 * D_max_rs^2))

    if (lambda_min_rs >= lambda_max_rs) lambda_min_rs <- lambda_max_rs * 0.01

    lambda_min <- lambda_min_rs / (RS^2)
  }

  coarse_grid <- seq(lambda_max, lambda_min, length.out = n_coarse)

  message("   Grid Parameters:")
  message(sprintf("     Upper boundary : %g", lambda_max))
  message(sprintf("     Lower boundary : %g", lambda_min))
  message(sprintf("     Coarse grid    : %d steps [%g -> %g]\n", n_coarse, lambda_max, lambda_min))

  spatialCovs_wrapped <- lapply(spatialCovs, terra::wrap)

  # --- Helper: Model Fitting & Scoring Function ---
  score_lambda <- function(lam, step_idx, total_steps) {
    prefix_msg <- sprintf("     [%d/%d] Testing lambda = %-8.4g ... ", step_idx, total_steps, lam)
    local_warnings <- character()

    res <- withCallingHandlers({
      spatialCovs_worker <- lapply(spatialCovs_wrapped, terra::unwrap)
      attr(spatialCovs_worker[[barrier]],"barLangevin") <- TRUE

      fit_args <- args
      fit_args$data <- data
      fit_args$spatialCovs <- spatialCovs_worker
      fit_args$barrier <- barrier
      fit_args$lambda <- lam
      fit_args$silent <- TRUE

      fit <- suppressMessages(tryCatch(do.call(fitLangevin, fit_args), error = function(e) e))

      if (inherits(fit, "error")) {
        err_msg <- if(inherits(fit, "error")) fit$message else paste("Convergence code", fit$convergence)
        message(paste0(prefix_msg, "Fit Failed (", err_msg, ")"))
        list(score = Inf, sd = Inf, fit = NULL)
      } else {
        sim_args <- list(model = fit, data = data, spatialCovs = spatialCovs_worker, timeStep = timeStep, conditional = FALSE)

        ks_scores <- c()
        err_msgs <- c()

        for (s in 1:n_sims) {
          sim_data <- suppressMessages(tryCatch(do.call(simLangevin, sim_args), error = function(e) e))

          if (inherits(sim_data, "error")) {
            err_msgs <- c(err_msgs, sim_data$message)
            next
          }

          sim_pts <- cbind(sim_data[[coord[1]]], sim_data[[coord[2]]])
          sim_dist <- terra::extract(spatialCovs_worker[[barrier]], sim_pts)[, 1]
          sim_dist <- sim_dist[!is.na(sim_dist)]

          if(length(sim_dist) == 0) {
            err_msgs <- c(err_msgs, "No valid spatial observations extracted")
            next
          }

          ks_scores <- c(ks_scores, suppressWarnings(stats::ks.test(obs_dist, sim_dist)$statistic))
        }

        if (length(ks_scores) == 0) {
          message(paste0(prefix_msg, "Simulations Failed (", err_msgs[1], ")"))
          list(score = Inf, sd = Inf, fit = fit)
        } else {
          overall_ks <- stats::median(ks_scores)
          sd_ks <- if(length(ks_scores) > 1) stats::sd(ks_scores) else 0.0

          message(paste0(prefix_msg, sprintf("Median KS = %.4f, SD = %.4f (%d/%d valid sims)", overall_ks, sd_ks, length(ks_scores), n_sims)))
          list(score = unname(overall_ks), sd = unname(sd_ks), fit = fit)
        }
      }
    }, warning = function(w) {
      local_warnings <<- c(local_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

    res$warnings <- local_warnings
    return(res)
  }

  # --- 2. Phase I: Top-Down Coarse Linear Search ---
  message("   --------------------------------------------------")
  message("    PHASE I: Coarse Linear Search")
  message("   --------------------------------------------------")

  coarse_results <- foreach::foreach(i = seq_along(coarse_grid), .packages = c("terra", "stats", "langevinSSM")) %loop% {
    score_lambda(coarse_grid[i], i, n_coarse)
  }

  coarse_scores <- sapply(coarse_results, function(x) x$score)
  valid_coarse_idx <- which(!is.infinite(coarse_scores))

  if (length(valid_coarse_idx) == 0) {
    stop("All models in the coarse grid failed. Check your data, initialization parameters, or provide a lower lambda_max.")
  }

  best_coarse_idx <- valid_coarse_idx[which.min(coarse_scores[valid_coarse_idx])]
  best_coarse_lam <- coarse_grid[best_coarse_idx]
  best_coarse_score <- coarse_scores[best_coarse_idx]
  best_coarse_sd <- coarse_results[[best_coarse_idx]]$sd

  global_warnings <- c(global_warnings, unlist(lapply(coarse_results, function(x) x$warnings)))

  message(sprintf("\n   => Phase I Winner: lambda = %g (Median KS = %.4f, SD = %.4f)\n", best_coarse_lam, best_coarse_score, best_coarse_sd))

  # --- 3. Phase II: Fine Linear Search ---
  message("   --------------------------------------------------")
  message("    PHASE II: Fine Linear Search")
  message("   --------------------------------------------------")

  step_size <- (lambda_max - lambda_min) / (n_coarse - 1)
  upper_bound <- best_coarse_lam + step_size
  lower_bound <- max(lambda_min, best_coarse_lam - step_size)

  fine_grid <- seq(upper_bound, lower_bound, length.out = n_fine + 2)[-c(1, n_fine + 2)]
  fine_grid <- fine_grid[abs(fine_grid - best_coarse_lam) > 1e-8]

  n_fine_actual <- length(fine_grid)
  if (n_fine_actual > 0) {
    fine_results <- foreach::foreach(i = seq_along(fine_grid), .packages = c("terra", "stats", "langevinSSM")) %loop% {
      score_lambda(fine_grid[i], i, n_fine_actual)
    }

    global_warnings <- c(global_warnings, unlist(lapply(fine_results, function(x) x$warnings)))
  } else {
    fine_results <- list()
  }

  fine_scores <- sapply(fine_results, function(x) x$score)
  valid_fine_idx <- which(!is.infinite(fine_scores))

  all_lams <- best_coarse_lam
  all_scores <- best_coarse_score
  all_sds <- best_coarse_sd
  all_fits <- list(coarse_results[[best_coarse_idx]]$fit)

  if (length(valid_fine_idx) > 0) {
    all_lams <- c(all_lams, fine_grid[valid_fine_idx])
    all_scores <- c(all_scores, fine_scores[valid_fine_idx])
    all_sds <- c(all_sds, sapply(fine_results[valid_fine_idx], function(x) x$sd))
    all_fits <- c(all_fits, lapply(fine_results[valid_fine_idx], function(x) x$fit))
  }

  best_final_idx <- which.min(all_scores)
  best_final_lam <- all_lams[best_final_idx]
  best_final_score <- all_scores[best_final_idx]
  best_final_sd <- all_sds[best_final_idx]
  best_final_fit <- all_fits[[best_final_idx]]

  message(sprintf("\n   =================================================="))
  message(sprintf("    OPTIMAL LAMBDA FOUND: %g (Median KS = %.4f, SD = %.4f)", best_final_lam, best_final_score, best_final_sd))
  message(sprintf("   =================================================="))

  global_warnings <- unique(global_warnings)
  if (length(global_warnings) > 0) {
    for (w in global_warnings) {
      warning(w, call. = FALSE, immediate. = TRUE)
    }
  }

  return(best_final_fit)
}
