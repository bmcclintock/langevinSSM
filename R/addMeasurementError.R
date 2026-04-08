#' @importFrom dplyr %>% mutate
addMeasurementError <- function(model, out, par, measurementError = NULL, exact = FALSE) {

  # Initialize x and y to true locations if they don't exist yet
  if (!"x" %in% names(out)) out$x <- out$mu.x
  if (!"y" %in% names(out)) out$y <- out$mu.y

  if (is.null(measurementError) && !exact) {
    # No error: observations equal true locations
    out <- out %>% dplyr::mutate(smaj = NA, smin = NA, eor = NA, x.sd = NA, y.sd = NA)
    return(out)
  }

  out$lc <- "G"

  # Robustly extract parameters (handles both natural and working scale inputs)
  psi <- if (!is.null(par$psi)) par$psi else if (!is.null(par$l_psi)) exp(par$l_psi) else 1
  tau <- if (!is.null(par$tau)) par$tau else if (!is.null(par$l_tau)) exp(par$l_tau) else c(1, 1)
  rho_o <- if (!is.null(par$rho_o)) par$rho_o else if (!is.null(par$l_rho_o)) (2 / (1 + exp(-par$l_rho_o)) - 1) else 0

  if (exact) {
    # EXACT MODE: Data already contains smaj/x.sd columns. We apply the math row-by-row.
    has_KF <- "smaj" %in% names(out) && any(!is.na(out$smaj))
    has_LS <- "x.sd" %in% names(out) && any(!is.na(out$x.sd))

    if (has_KF) {
      res <- measurementError_rcpp(out, 0, 0, c(0,0), psi, ifelse(model == "underdamped", 1, 0), TRUE)
      out$x <- res$x
      out$y <- res$y
    }

    if (has_LS) {
      res <- measurementError_LS_rcpp(out, 0, 0, tau[1], tau[2], rho_o, ifelse(model == "underdamped", 1, 0), TRUE)
      out$x <- res$x
      out$y <- res$y
    }

  } else {
    # DE NOVO MODE: Generate completely new errors
    if (!is.list(measurementError)) stop("'measurementError' must be a list.")

    has_KF <- all(c("smaj.sd", "smin.sd") %in% names(measurementError))
    has_LS <- all(c("x.sd", "y.sd") %in% names(measurementError))

    if (has_KF && has_LS) {
      stop("Cannot provide both Argos KF (smaj.sd, smin.sd, eor) and LS/GPS (x.sd, y.sd) error parameters.")
    }

    if (has_KF) {
      smaj.sd <- measurementError$smaj.sd
      smin.sd <- measurementError$smin.sd
      eor_range <- if (!is.null(measurementError$eor)) measurementError$eor else c(0, 180)

      res <- measurementError_rcpp(out, smaj.sd, smin.sd, eor_range, psi, ifelse(model == "underdamped", 1, 0), FALSE)

      out$x <- res$x
      out$y <- res$y
      out$smaj <- res$smaj
      out$smin <- res$smin
      out$eor <- res$eor
      out$x.sd <- NA
      out$y.sd <- NA

    } else if (has_LS) {
      x_sd_val <- measurementError$x.sd
      y_sd_val <- measurementError$y.sd

      res <- measurementError_LS_rcpp(out, x_sd_val, y_sd_val, tau[1], tau[2], rho_o, ifelse(model == "underdamped", 1, 0), FALSE)

      out$x <- res$x
      out$y <- res$y
      out$x.sd <- res$x.sd
      out$y.sd <- res$y.sd
      out$smaj <- NA
      out$smin <- NA
      out$eor <- NA
    }
  }
  return(out)
}
