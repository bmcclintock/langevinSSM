#' Extract natural scale parameters from a fitted Langevin model
#'
#' @description
#' Extracts the fitted fixed and random effect parameter estimates from a \code{fitLangevin} object. The parameters are returned on their natural, unscaled scale in a list format perfectly structured to be supplied to the \code{par} argument of \code{\link{fitLangevin}} for warm-starting or refitting models.
#'
#' @param fit A \code{fitLangevin} object.
#'
#' @return A list of natural scale parameter estimates (e.g., \code{beta}, \code{sigma}, \code{gamma}, \code{mu}, \code{vel}, and any observation error parameters like \code{psi}, \code{tau}, or \code{rho_o} if they were included in the model).
#' @export
getPar <- function(fit) {

  if (!inherits(fit, "fitLangevin")) {
    stop("'fit' must be a 'fitLangevin' object.")
  }

  scaleFactor <- fit$conditions$scaleFactor
  if (is.null(scaleFactor)) scaleFactor <- 1

  parList <- fit$tmb_setup$parList
  model <- fit$conditions$model

  par_out <- list()

  par_out$beta <- parList$beta
  par_out$sigma <- exp(parList$log_sigma) * scaleFactor

  if (model == "underdamped") {
    par_out$gamma <- exp(parList$log_gamma)
  }

  parNames <- names(fit$conditions$par)

  if ("l_psi" %in% parNames) {
    par_out$psi <- exp(parList$l_psi)
  }
  if ("l_tau" %in% parNames) {
    par_out$tau <- exp(parList$l_tau)
  }
  if ("l_rho_o" %in% parNames) {
    par_out$rho_o <- 2 / (1 + exp(-parList$l_rho_o)) - 1
  }

  if ("mu" %in% names(parList)) {
    par_out$mu <- unname(as.matrix(t(parList$mu) * scaleFactor))
  }

  if (model == "underdamped" && "vel" %in% names(parList)) {
    par_out$vel <- unname(as.matrix(t(parList$vel) * scaleFactor))
  }

  return(par_out)
}
