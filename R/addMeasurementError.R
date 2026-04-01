#' @importFrom dplyr %>% mutate
addMeasurementError <- function(model, out, par, measurementError) {
  if (is.null(measurementError)) {
    # No error: observations equal true locations
    out <- out %>% dplyr::mutate(x = mu.x, y = mu.y, smaj = NA, smin = NA,
                                 eor = NA, x.sd = NA, y.sd = NA)
    return(out)
  }

  if (!is.list(measurementError)) stop("'measurementError' must be a list.")

  # Check for mutual exclusivity
  has_KF <- all(c("smaj.sd", "smin.sd") %in% names(measurementError))
  has_LS <- all(c("x.sd", "y.sd") %in% names(measurementError))

  if (has_KF && has_LS) {
    stop("Cannot provide both Argos KF (smaj.sd, smin.sd, eor) and LS/GPS (x.sd, y.sd) error parameters.")
  }

  n <- nrow(out)
  out$lc <- "G"

  if (has_KF) {
    # --- Argos KF ---
    smaj.sd <- measurementError$smaj.sd
    smin.sd <- measurementError$smin.sd
    eor_range <- if (!is.null(measurementError$eor)) measurementError$eor else c(0, 180)
    psi <- if (!is.null(par$l_psi)) exp(par$l_psi) else 1

    # Use your existing Rcpp function for KF sampling
    out <- measurementError_rcpp(out, smaj_sd = smaj.sd, smin_sd = smin.sd, eor = eor_range, psi = psi,
                                 model = ifelse(model == "underdamped", 1, 0))
    out$x.sd <- NA
    out$y.sd <- NA

  } else if (has_LS) {
    # --- LS/GPS ---
    x_sd_val <- measurementError$x.sd
    y_sd_val <- measurementError$y.sd
    tau <- if (!is.null(par$l_tau)) exp(par$l_tau) else c(1, 1)
    rho_o <- if (!is.null(par$l_rho_o)) (2 / (1 + exp(-par$l_rho_o)) - 1) else 0

    out <- measurementError_LS_rcpp(data = out,
                                    x_sd = x_sd_val,
                                    y_sd = y_sd_val,
                                    tau_x = tau[1],
                                    tau_y = tau[2],
                                    rho_o = rho_o,
                                    model = ifelse(model == "underdamped", 1, 0))

    out$smaj <- NA
    out$smin <- NA
    out$eor <- NA
  }
  return(out)
}
