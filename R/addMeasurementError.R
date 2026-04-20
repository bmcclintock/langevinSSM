#' Add measurement error to true locations
#'
#' This function adds measurement error to the true locations in the data frame, either by simulating new errors using provided measurement error parameters or by applying existing error columns natively found in the data. It supports both Argos Kalman Filter (KF) and Least Squares (LS)/GPS error models.
#' @param data A data frame containing the true locations (columns specified by `coord`) and optionally measurement error data (e.g., `smaj`, `smin`, `eor` for KF or `x.err`, `y.err` for LS/GPS).
#' @param par A list of parameters for the error model. For KF, this can include `psi`. For LS/GPS, this can include `tau` (as a vector of length 2) and `rho_o`. See \code{\link{fitLangevin}}.
#' @param measurementError A list or data frame used to simulate observation error. The structure determines the error model used:
#' \itemize{
#'   \item \strong{Argos Kalman Filter (KF):} A list containing \code{smaj.sd} (numeric; the SD to generate the semi-major axis), \code{smin.sd} (numeric; the SD to generate the semi-minor axis), and optionally \code{eor.lim} (numeric vector of length 2; boundaries in degrees for uniform orientation, defaulting to \code{c(0, 180)}).
#'   \item \strong{Least Squares (LS) or GPS:} A list containing \code{x.sd} and \code{y.sd} (numeric; the SDs of the half-normal distributions used to randomly generate \code{x.err} and \code{y.err}).
#'   \item \strong{Location Quality Class (EMF):} A data frame (e.g., as returned by \code{\link{getEMF}}) with an additional column \code{prob} that sums to 1. The function will randomly assign a location class (\code{lc}) to each observation based on these probabilities and simulate \code{x.err} and \code{y.err} using the corresponding \code{emf.x} and \code{emf.y} values.
#' }
#' You cannot provide parameters for multiple error models simultaneously. If \code{NULL}, the function assumes known error magnitudes already exist within \code{data} and applies the error covariance matrix row-by-row.
#' @param coord A character vector of length 2 specifying the column names in `data` that contain the true x and y locations (default is `c("mu.x", "mu.y")`).
#' @return A data frame with observed locations (`x`, `y`). and measurement error terms added.
#' @importFrom dplyr %>% mutate
#' @importFrom stats rnorm
#' @export
addMeasurementError <- function(data, par = NULL, measurementError = NULL, coord = c("mu.x", "mu.y")) {

  # Initialize x and y to true locations if they don't exist yet
  if (!"x" %in% names(data)) data$x <- data[,coord[1]]
  if (!"y" %in% names(data)) data$y <- data[,coord[2]]

  if (is.data.frame(measurementError)) {
    if (!all(c("lc", "emf.x", "emf.y", "prob") %in% names(measurementError))) {
      stop("If 'measurementError' is a data frame, it must contain columns 'lc', 'emf.x', 'emf.y', and 'prob'.")
    }
    if (abs(sum(measurementError$prob) - 1) > 1e-6) {
      stop("The 'prob' column in the 'measurementError' data frame must sum to 1.")
    }

    n_obs <- nrow(data)
    drawn_lc <- sample(measurementError$lc, size = n_obs, replace = TRUE, prob = measurementError$prob)

    match_idx <- match(drawn_lc, measurementError$lc)
    drawn_emf_x <- measurementError$emf.x[match_idx]
    drawn_emf_y <- measurementError$emf.y[match_idx]

    data$lc <- drawn_lc
    data$x.err <- abs(stats::rnorm(n_obs, 0, drawn_emf_x))
    data$y.err <- abs(stats::rnorm(n_obs, 0, drawn_emf_y))

    data$smaj <- NA
    data$smin <- NA
    data$eor <- NA

    # Nullify measurementError so it behaves as knownError from here on out
    measurementError <- NULL
  }
  # -----------------------------------

  knownError <- is.null(measurementError)

  if(knownError & !("smaj" %in% names(data) && "smin" %in% names(data) && "eor" %in% names(data)) &
     !("x.err" %in% names(data) && "y.err" %in% names(data))) {
    stop("No measurement error parameters provided.\nPlease provide either 'measurementError' or appropriate measurement error columns in 'data' (i.e., 'smaj', 'smin', 'eor', 'x.err' and 'y.err').")
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
    # Data already contains smaj/x.err columns. We apply the math row-by-row.
    has_KF <- "smaj" %in% names(data) && any(!is.na(data$smaj))
    has_LS <- "x.err" %in% names(data) && any(!is.na(data$x.err))

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
      data$x.err <- NA
      data$y.err <- NA
    }

  } else {
    # Generate completely new errors
    if (!is.list(measurementError)) stop("'measurementError' must be a list.")

    has_KF <- all(c("smaj.sd", "smin.sd") %in% names(measurementError))
    has_LS <- all(c("x.sd", "y.sd") %in% names(measurementError))

    if (has_KF && has_LS) {
      stop("Cannot provide both Argos KF (smaj.sd, smin.sd, eor.lim) and LS/GPS (x.sd, y.sd) error parameters.")
    }

    if (has_KF) {
      smaj.sd <- measurementError$smaj.sd
      smin.sd <- measurementError$smin.sd
      eor_range <- if (!is.null(measurementError$eor.lim)) measurementError$eor.lim else c(0, 180)

      res <- measurementError_rcpp(data, smaj.sd, smin.sd, eor_range, psi, FALSE)

      data$x <- res$x
      data$y <- res$y
      data$smaj <- res$smaj
      data$smin <- res$smin
      data$eor <- res$eor
      data$x.err <- NA
      data$y.err <- NA

    } else if (has_LS) {
      x_sd_val <- measurementError$x.sd
      y_sd_val <- measurementError$y.sd

      res <- measurementError_LS_rcpp(data, x_sd_val, y_sd_val, tau[1], tau[2], rho_o, FALSE)

      data$x <- res$x
      data$y <- res$y

      data$x.err <- res$x.err
      data$y.err <- res$y.err

      data$smaj <- NA
      data$smin <- NA
      data$eor <- NA
    }
  }
  return(data)
}
