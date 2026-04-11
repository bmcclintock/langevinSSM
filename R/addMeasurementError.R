#' Add measurement error to true locations
#'
#' This function adds measurement error to the true locations in the data frame, either by using provided measurement error parameters or by using existing error columns in the data. It supports both Argos Kalman Filter (KF) and Least Squares (LS)/GPS error models
#' @param data A data frame containing the true locations (columns specified by `coord`) and optionally measurement error data (e.g., `smaj`, `smin`, `eor` for KF or `x.sd`, `y.sd` for LS/GPS).
#' @param par A list of parameters for the error model. For KF, this can include `psi`. For LS/GPS, this can include `tau` (as a vector of length 2) and `rho_o`. See \code{\link{fitLangevin}}.
#' @param measurementError A list of measurement error parameters. For KF, this should include `smaj.sd`, `smin.sd`, and optionally `eor`. For LS/GPS, this should include `x.sd` and `y.sd`. If not provided, the function will look for appropriate columns in `data`.
#' @param coord A character vector of length 2 specifying the column names in `data` that contain the true x and y locations (default is `c("mu.x", "mu.y")`).
#' @return A data frame with observed locations (`x`, `y`). and measurement error terms added (if applicable).
#' @importFrom dplyr %>% mutate
#' @export
addMeasurementError <- function(data, par = NULL, measurementError = NULL, coord = c("mu.x", "mu.y")) {

  # Initialize x and y to true locations if they don't exist yet
  if (!"x" %in% names(data)) data$x <- data[,coord[1]]
  if (!"y" %in% names(data)) data$y <- data[,coord[2]]

  knownError <- is.null(measurementError)

  if(knownError & !("smaj" %in% names(data) && "smin" %in% names(data) && "eor" %in% names(data)) &
     !("x.sd" %in% names(data) && "y.sd" %in% names(data))) {
    stop("No measurement error parameters provided.\nPlease provide either 'measurementError' or appropriate measurement error columns in 'data' (i.e., 'smaj', 'smin', 'eor', 'x.sd' and 'y.sd').")
  }

  if(!is.null(par$psi) & !is.null(par$l_psi)) stop("Cannot provide both 'psi' and 'l_psi' in 'par'.")
  if(!is.null(par$tau) & !is.null(par$l_tau)) stop("Cannot provide both 'tau' and 'l_tau' in 'par'.")
  if(!is.null(par$rho_o) & !is.null(par$l_rho_o)) stop("Cannot provide both 'rho_o' and 'l_rho_o' in 'par'.")

  # Robustly extract parameters (handles both natural and working scale inputs)
  psi <- if (!is.null(par$psi)) par$psi else if (!is.null(par$l_psi)) exp(par$l_psi) else 1
  tau <- if (!is.null(par$tau)) par$tau else if (!is.null(par$l_tau)) exp(par$l_tau) else c(1, 1)
  rho_o <- if (!is.null(par$rho_o)) par$rho_o else if (!is.null(par$l_rho_o)) (2 / (1 + exp(-par$l_rho_o)) - 1) else 0

  if(knownError & ("smaj" %in% names(data) && "smin" %in% names(data) && "eor" %in% names(data))) {
    if(!all(is.na(data$eor)) && max(data$eor, na.rm = TRUE) < pi) {
     warning("eor values were converted to radians, but they appear to have been provided in radians rather than degrees from north. Please ensure that the eor column was provided in degrees from north.")
    }
    data$eor <- data$eor * pi / 180 # convert from degrees to radians
  }

  checkErrorData(data, coord, measurementError, knownError)

  if (knownError) {
    # Data already contains smaj/x.sd columns. We apply the math row-by-row.
    has_KF <- "smaj" %in% names(data) && any(!is.na(data$smaj))
    has_LS <- "x.sd" %in% names(data) && any(!is.na(data$x.sd))

    if (has_KF) {
      res <- measurementError_rcpp(data, 0, 0, c(0,0), psi, TRUE)
      data$x <- res$x
      data$y <- res$y
    } else {
      data$smaj <- NA
      data$smin <- NA
      data$eor <- NA
    }

    if (has_LS) {
      res <- measurementError_LS_rcpp(data, 0, 0, tau[1], tau[2], rho_o, TRUE)
      data$x <- res$x
      data$y <- res$y
    } else {
      data$x.sd <- NA
      data$y.sd <- NA
    }

  } else {
    # Generate completely new errors
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

      res <- measurementError_rcpp(data, smaj.sd, smin.sd, eor_range, psi, FALSE)

      data$x <- res$x
      data$y <- res$y
      data$smaj <- res$smaj
      data$smin <- res$smin
      data$eor <- res$eor
      data$x.sd <- NA
      data$y.sd <- NA

    } else if (has_LS) {
      x_sd_val <- measurementError$x.sd
      y_sd_val <- measurementError$y.sd

      res <- measurementError_LS_rcpp(data, x_sd_val, y_sd_val, tau[1], tau[2], rho_o, FALSE)

      data$x <- res$x
      data$y <- res$y
      data$x.sd <- res$x.sd
      data$y.sd <- res$y.sd
      data$smaj <- NA
      data$smin <- NA
      data$eor <- NA
    }
  }
  return(data)
}
