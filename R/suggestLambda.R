#' Suggest a barrier penalty
#'
#' Fits a baseline unconstrained model to estimate the animal's speed and uses the maximum time step
#' to suggest a stable upper limit for the barrier penalty (\code{lambda}). This is just a suggestion and is not necessarily the optimal value, but it can help guide users in selecting a reasonable starting point for \code{lambda} when fitting models with barriers.
#'
#' @param data A \code{dataLangevin} object containing the tracking data.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates.
#' @param barrier Character string. The name of the barrier in \code{spatialCovs} that is represented as a signed distance field (see \code{\link{prepBarrier}}). If missing, it will attempt to extract the barrier name from the data attributes.
#' @param ... Additional arguments passed to \code{\link{fitLangevin}}.
#'
#' @return A numeric value representing the suggested maximum barrier penalty (\code{lambda}).
#' @export
suggestLambda <- function(data, spatialCovs, barrier, ...) {

  if (!inherits(data, "dataLangevin")) stop("'data' must be a dataLangevin object.")
  if (missing(spatialCovs) || !is.list(spatialCovs)) stop("'spatialCovs' list must be provided.")

  if (missing(barrier) || is.null(barrier)) {
    sim_barrier <- attr(data,"barrier")
    if(!is.null(sim_barrier)){
      barrier <- sim_barrier
    } else stop("'barrier' argument (the name of the barrier raster in 'spatialCovs') must be provided.")
  }

  max_dt <- max(data$dt, na.rm = TRUE)

  args <- list(...)
  fit0_args <- args
  fit0_args$data <- data
  fit0_args$spatialCovs <- spatialCovs
  fit0_args$barrier <- barrier
  fit0_args$lambda <- 0

  if(is.null(fit0_args$model)) model <- "underdamped" else model <- fit0_args$model
  fit0_args$model <- model

  message("   Fitting baseline ",model," Langevin model (lambda = 0) to estimate movement parameters...")

  fit0 <- tryCatch(suppressMessages(do.call(fitLangevin, fit0_args)), error = function(e) e)

  if (inherits(fit0, "error")) {
    stop("Baseline fit failed.\n    Error details: ", fit0$message)
  }

  print(fit0)

  scaleFactor <- fit0$conditions$scaleFactor

  est_sigma <- fit0$estimates$natural["sigma", "Estimate"]
  est_sigma_work <- est_sigma / scaleFactor

  # Normalize the spatial parameters, calculate the boundaries, and then apply the 1/L^2 transformation
  RS <- 10^(floor(log10(est_sigma_work)))

  sigma_rs <- est_sigma_work / RS

  if (model == "underdamped") {
    est_gamma <- fit0$estimates$natural["gamma", "Estimate"]

    num <- (est_gamma^2) * (1 - exp(-est_gamma * max_dt))
    den <- (sigma_rs^2) * (1 - exp(-est_gamma * max_dt) - (est_gamma * max_dt * exp(-est_gamma * max_dt)))
    lambda_max_rs <- num / den
  } else {
    lambda_max_rs <- 2 / (sigma_rs^2 * max_dt)
  }

  lambda_max <- (lambda_max_rs / (RS^2))

  message(sprintf("   => Suggested maximum barrier penalty (lambda): %g\n", lambda_max))

  return(unname(lambda_max))
}
