#' Extract natural scale parameters from a fitted Langevin model
#'
#' @description
#' Extracts the fitted fixed and random effect parameter estimates from a \code{fitLangevin} object. The parameters are returned on their natural scale in the list format required for the \code{par} argument of \code{\link{fitLangevin}}.
#'
#' @param fit A \code{fitLangevin} object.
#'
#' @return A list of natural scale parameter estimates (e.g., \code{beta}, \code{sigma}, \code{gamma}, \code{mu}, \code{vel} and any observation error parameters (e.g., \code{psi}, \code{tau}, \code{rho_o}) that were estimated).
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

  # Extract names of actively estimated parameters from the outer optimizer
  active_pars <- names(fit$par)

  if ("beta" %in% active_pars) {
    par_out$beta <- parList$beta
  }

  if ("log_sigma" %in% active_pars) {
    par_out$sigma <- exp(parList$log_sigma) * scaleFactor
  }

  if (model == "underdamped" && "log_gamma" %in% active_pars) {
    par_out$gamma <- exp(parList$log_gamma)
  }

  if ("l_psi" %in% active_pars) {
    par_out$psi <- exp(parList$l_psi)
  }
  if ("l_tau" %in% active_pars) {
    par_out$tau <- exp(parList$l_tau)
  }
  if ("l_rho_o" %in% active_pars) {
    par_out$rho_o <- 2 / (1 + exp(-parList$l_rho_o)) - 1
  }

  par_out$mu <- unname(as.matrix(t(parList$mu) * scaleFactor))

  if (model == "underdamped") {
    par_out$vel <- unname(as.matrix(t(parList$vel) * scaleFactor))
  }

  return(par_out)
}
