#' Suggest a barrier penalty
#'
#' Uses the estimated speed/variance from a baseline unconstrained model and a user-specified maximum time step
#' to suggest a stable upper limit for the barrier penalty (\code{lambda}). This is just a suggestion and is not necessarily the optimal value, but it can help guide users in selecting a reasonable starting point for \code{lambda} when fitting models with barriers.
#'
#' @param fit A \code{fitLangevin} object representing the baseline unconstrained model (fitted without a barrier constraint, or with \code{lambda = 0}).
#' @param max_dt Numeric. A representative maximum time step to use for the calculation. This must be provided by the user. Supplying a value representing the typical sampling rate, a high quantile of the observed time steps, or the intended prediction grid resolution is recommended over the absolute maximum observed time step, which can be heavily skewed by extreme outlier gaps.
#'
#' @details
#' The maximum stable penalty is derived from the animal's movement characteristics (estimated by the baseline \code{fitLangevin} model) and the largest time gap between observations (\code{max_dt}). It attempts to balance the movement variance (\code{sigma}), the friction coefficient (\code{gamma}, for underdamped models), and the maximum time step. Fundamentally, if the time gaps between observations are large, or if the animal's movement variance is high, the barrier penalty must be kept relatively small to prevent the trajectory from overshooting when it encounters the barrier.
#'
#' The suggested maximum penalty \eqn{\lambda_{max}} is calculated based on the chosen movement model. For the underdamped model:
#'
#' \deqn{\lambda_{max}=\frac{\gamma^2(1-e^{-\gamma\Delta t_{max}})}{\sigma_{work}^2(1-e^{-\gamma\Delta t_{max}}-\gamma\Delta t_{max}e^{-\gamma\Delta t_{max}})}}
#'
#' For the overdamped model:
#'
#' \deqn{\lambda_{max}=\frac{2}{\sigma_{work}^2\Delta t_{max}}}
#'
#' where \eqn{\gamma} is the estimated friction coefficient, \eqn{\Delta t_{max}} is the user-specified maximum time step (\code{max_dt}), and \eqn{\sigma_{work}} is the estimated movement variance adjusted to the internal working scale (\eqn{\sigma/\text{scaleFactor}}).
#'
#' @return A numeric value representing the suggested maximum barrier penalty (\code{lambda}).
#' @export
suggestLambda <- function(fit, max_dt) {

  if (missing(fit) || !inherits(fit, "fitLangevin")) stop("'fit' must be provided and must be a fitLangevin object.")

  # Ensure the provided model was actually unconstrained
  if (!is.null(fit$conditions$lambda) && fit$conditions$lambda != 0) {
    stop("The provided 'fit' object must be an unconstrained model fitted with lambda = 0 (or no barrier).")
  }

  if (missing(max_dt) || !is.numeric(max_dt) || length(max_dt) != 1 || max_dt <= 0) {
    stop("'max_dt' must be provided as a single positive numeric value.")
  }

  model <- fit$conditions$model
  scaleFactor <- fit$conditions$scaleFactor

  est_sigma <- fit$estimates$natural["sigma", "Estimate"]
  est_sigma_work <- est_sigma / scaleFactor

  # Normalize the spatial parameters, calculate the boundaries, and then apply the 1/L^2 transformation
  RS <- 10^(floor(log10(est_sigma_work)))

  sigma_rs <- est_sigma_work / RS

  if (model == "underdamped") {
    est_gamma <- fit$estimates$natural["gamma", "Estimate"]

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
