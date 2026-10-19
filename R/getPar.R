#' Extract natural scale parameters from a fitted Langevin model
#'
#' @description
#' Extracts the fitted fixed and random effect parameter estimates from a \code{fitLangevin} object. The parameters are returned on their natural scale in the list format required for the \code{par} argument of \code{\link{fitLangevin}}.
#'
#' @param fit A \code{fitLangevin} object.
#'
#' @return A list of natural scale parameter estimates (e.g., \code{beta}, \code{sigma}, \code{gamma}, \code{lambda} (if the barrier strength was estimated, i.e. \code{lambda = "estimate"}), \code{mu}, \code{vel} and any observation error parameters (e.g., \code{psi}, \code{psi_maj}, \code{b_maj}, \code{b_min}, \code{tau}, \code{rho_o}, \code{df}, \code{k}) that were estimated). Parameters specified by location class (or observation type) in \code{\link{fitLangevin}} are returned in the same form: \code{tau} as a 2-row matrix with one column per location class, \code{df} or \code{k} as a named vector.
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

  # barrier strength (natural scale), only a parameter when fitted with lambda = "estimate"
  if ("log_lambda" %in% active_pars || identical(fit$conditions$lambda_mode, "estimate")) {
    par_out$lambda <- exp(parList$log_lambda) / scaleFactor^2
  }

  spec <- .cond_obs_spec(fit$conditions)
  if ("l_psi" %in% active_pars) {
    par_out$psi <- exp(parList$l_psi)
  }
  for (nm in .ell_par_names) {
    if (paste0("l_", nm) %in% active_pars) par_out[[nm]] <- exp(parList[[paste0("l_", nm)]])
  }
  if ("l_tau" %in% active_pars) {
    par_out$tau <- .label_obs_par("tau", exp(parList$l_tau), spec)
  }
  if ("l_rho_o" %in% active_pars) {
    par_out$rho_o <- 2 / (1 + exp(-parList$l_rho_o)) - 1
  }
  if ("l_df" %in% active_pars) {
    par_out$df <- .label_obs_par("df", exp(parList$l_df), spec)
  }
  if ("l_k" %in% active_pars) {
    par_out$k <- .label_obs_par("k", exp(parList$l_k), spec)
  }

  par_out$mu <- unname(as.matrix(t(parList$mu) * scaleFactor))

  if (model == "underdamped") {
    par_out$vel <- unname(as.matrix(t(parList$vel) * scaleFactor))
  }

  return(par_out)
}
