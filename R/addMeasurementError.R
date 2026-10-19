#' Add measurement error to true locations
#'
#' This function adds measurement error to the true locations in the data frame, either by simulating new errors using provided measurement error parameters or by applying existing error columns natively found in the data. It supports both Argos Kalman Filter (KF) and Least Squares (LS)/GPS error models.
#' @param data A \code{dataLangevin} object containing the true locations (columns specified by `coord`) and optionally measurement error data (e.g., `smaj`, `smin`, `eor` for KF or `x.err`, `y.err` for LS/GPS).
#' @param par A list of parameters for the error model. For KF, this can include `psi` and the error ellipse power law parameters `psi_maj`, `b_maj` and `b_min`. For LS/GPS or EMF, this can include `tau` (as a vector of length 2, or by location class (and/or "sderr") as a 2-row matrix or named list) and `rho_o`. If (some) observations have \code{obs_dist="t"}, then \code{par} can also include the degrees of freedom (`df`) parameter. If (some) observations have \code{obs_dist="Huber"}, then \code{par} can also include the threshold (`k`) parameter. `df` and `k` can be a single value or a named vector by error type ("epar", "sderr") and/or location class. See \code{\link{fitLangevin}}.
#' @param measurementError A list or data frame used to simulate observation error. The structure determines the error model used:
#' \itemize{
#'   \item \strong{Argos Kalman Filter (KF):} A list containing \code{smaj.sd} (numeric; the SD to generate the semi-major axis), \code{smin.sd} (numeric; the SD to generate the semi-minor axis), and optionally \code{eor.lim} (numeric vector of length 2; boundaries in degrees for uniform orientation, defaulting to \code{c(0, 180)}).
#'   \item \strong{Least Squares (LS) or GPS:} A list containing \code{x.sd} and \code{y.sd} (numeric; the SDs of the half-normal distributions used to randomly generate \code{x.err} and \code{y.err}).
#'   \item \strong{Location Quality Class (EMF):} A data frame (e.g., as returned by \code{\link{getEMF}}) with an additional column \code{prob} that sums to 1. The function will randomly assign a location class (\code{lc}) to each observation based on these probabilities and simulate \code{x.err} and \code{y.err} using the corresponding \code{emf.x} and \code{emf.y} values.
#' }
#' You cannot provide parameters for multiple error models simultaneously. If \code{NULL}, the function assumes known error magnitudes already exist within \code{data} and applies the error covariance matrix row-by-row.
#' @param coord A character vector of length 2 specifying the column names in `data` that contain the true x and y locations (default is `c("mu.x", "mu.y")`).
#' @param obs_dist The measurement error distribution: ``Gaussian'' (default), ``t'', or ``Huber'' (pseudo-Huber) for all observations, or a named character vector by error type (\code{"epar"}, \code{"sderr"}) and/or location class, e.g. \code{c(epar = "Huber", sderr = "Gaussian")}. See \code{\link{fitLangevin}}.
#' @return A data frame with observed locations (`x`, `y`). and measurement error terms added.
#' @importFrom dplyr %>% mutate
#' @importFrom stats rnorm
#' @export
addMeasurementError <- function(data, par = NULL, measurementError = NULL, coord = c("mu.x", "mu.y"), obs_dist = c("Gaussian", "t", "Huber")) {

  obs_dist <- .norm_obs_dist(obs_dist)

  if (!inherits(data, "dataLangevin")) {
    stop("'data' must be a 'dataLangevin' object (as returned by formatData or simLangevin.")
  }

  # Initialize x and y to true locations if they don't exist yet
  if (!"x" %in% names(data)) data$x <- data[,coord[1]]
  if (!"y" %in% names(data)) data$y <- data[,coord[2]]

  if(!is.null(par$psi) & !is.null(par$l_psi)) stop("Cannot provide both 'psi' and 'l_psi' in 'par'.")
  if(!is.null(par$tau) & !is.null(par$l_tau)) stop("Cannot provide both 'tau' and 'l_tau' in 'par'.")
  if(!is.null(par$rho_o) & !is.null(par$l_rho_o)) stop("Cannot provide both 'rho_o' and 'l_rho_o' in 'par'.")
  if(!is.null(par$df) & !is.null(par$l_df)) stop("Cannot provide both 'df' and 'l_df' in 'par'.")
  if(!is.null(par$k) & !is.null(par$l_k)) stop("Cannot provide both 'k' and 'l_k' in 'par'.")
  for (nm in .ell_par_names) if (!is.null(par[[nm]]) & !is.null(par[[paste0("l_", nm)]])) stop("Cannot provide both '", nm, "' and 'l_", nm, "' in 'par'.")

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

  # Extract parameters (handles both natural and working scale inputs)
  nat <- function(nm, default, inv = exp) if (!is.null(par[[nm]])) par[[nm]] else if (!is.null(par[[paste0("l_", nm)]])) inv(par[[paste0("l_", nm)]]) else default
  psi <- nat("psi", 1)
  ell <- vapply(.ell_par_names, function(nm) as.numeric(nat(nm, 1)), numeric(1))
  ell_power <- any(c(.ell_par_names, paste0("l_", .ell_par_names)) %in% names(par))
  if (any(!is.finite(ell)) || any(ell <= 0) || length(psi) != 1) stop("'psi', 'psi_maj', 'b_maj' and 'b_min' must be single positive values.")
  ell_ref <- attr(par, "ell_ref")   # reference lengths of the power law (internal; default: geometric means of the axes)
  if (is.null(ell_ref)) ell_ref <- c(maj = 0, min = 0)
  tau <- .norm_tau(nat("tau", c(1, 1)))
  rho_o <- nat("rho_o", 0, function(x) 2 / (1 + exp(-x)) - 1)
  df <- nat("df", 5); if (is.list(df)) df <- unlist(df)
  k <- nat("k", 2); if (is.list(k)) k <- unlist(k)
  if (any(!is.finite(tau)) || any(tau <= 0)) stop("'tau' must be greater than zero.")
  if (any(!is.finite(df)) || any(df <= 0)) stop("'df' must be greater than zero.")
  if (any(!is.finite(k)) || any(k <= 0)) stop("'k' must be greater than zero.")
  grouped <- function(x) !is.null(names(x)) || length(x) > 1
  if (grouped(df)) .check_labels(names(df), "par$df")
  if (grouped(k)) .check_labels(names(k), "par$k")
  if (is.matrix(tau)) .check_labels(colnames(tau), "par$tau", types = "sderr")

  # per-observation values of grouped entries (by error type and/or location class)
  by_group <- function(x, labels, rows, what) {
    n <- nrow(data)
    if (is.null(labels)) return(rep(as.numeric(x[1]), n))
    lc <- .data_lc(data)
    type <- .obs_type(data, coord = c("x", "y"))
    g <- .group_index(labels, lc, type)
    if (any(rows & is.na(g))) stop("'", what, "' is specified by error type or location class, but some observations (",
                                   paste(unique(paste0(type[rows & is.na(g)], ifelse(is.na(lc[rows & is.na(g)]), "", paste0(", lc ", lc[rows & is.na(g)])))), collapse = "; "),
                                   ") have no matching entry.")
    v <- x[g]
    v[is.na(v)] <- x[1]
    v
  }
  me_rows <- function() !is.na(.obs_type(data, coord = c("x", "y")))
  tau_vec <- function() {
    rows <- .obs_type(data, coord = c("x", "y")) %in% "sderr"
    if (!is.matrix(tau)) return(list(x = tau[1], y = tau[2]))
    list(x = as.numeric(by_group(tau[1, ], colnames(tau), rows, "tau")), y = as.numeric(by_group(tau[2, ], colnames(tau), rows, "tau")))
  }
  # distribution code (0 Gaussian, 1 t, 2 Huber) and shape parameter (df or k) of each observation
  dist_vec <- function() {
    if (is.null(names(obs_dist))) return(rep(match(obs_dist, .dist_choices) - 1L, nrow(data)))
    d <- by_group(match(unname(obs_dist), .dist_choices) - 1L, names(obs_dist), me_rows(), "obs_dist")
    d[is.na(d)] <- 0L
    as.integer(d)
  }
  shape_vec <- function(code) {
    rows <- me_rows()
    out <- rep(1, nrow(data))
    if (any(code == 1L)) out[code == 1L] <- as.numeric(by_group(df, if (grouped(df)) names(df) else NULL, rows & code == 1L, "df"))[code == 1L]
    if (any(code == 2L)) out[code == 2L] <- as.numeric(by_group(k, if (grouped(k)) names(k) else NULL, rows & code == 2L, "k"))[code == 2L]
    out
  }
  dist_args <- function() { code <- dist_vec(); list(code = code, shape = shape_vec(code)) }

  checkErrorData(data, coord, measurementError, knownError)

  if (knownError) {
    # Data already contains smaj/x.err columns. We apply the math row-by-row.
    has_KF <- "smaj" %in% names(data) && any(!is.na(data$smaj))
    has_LS <- "x.err" %in% names(data) && any(!is.na(data$x.err))

    if (has_KF) {
      if (!"x.err" %in% names(data)) { data$x.err <- NA; data$y.err <- NA }
      da <- dist_args()
      res <- measurementError_rcpp(data, 0, 0, c(0,0), psi, da$code, da$shape, TRUE,
                                   ell_power, ell[["psi_maj"]], ell[["b_maj"]], ell[["b_min"]], ell_ref[["maj"]], ell_ref[["min"]])
      data$x <- res$x
      data$y <- res$y
    } else {
      data$smaj <- NA
      data$smin <- NA
      data$eor <- NA
    }

    if (has_LS) {
      if (!"smaj" %in% names(data)) { data$smaj <- NA; data$smin <- NA; data$eor <- NA }
      tv <- tau_vec()
      da <- dist_args()
      res <- measurementError_LS_rcpp(data, 0, 0, tv$x, tv$y, rho_o, da$code, da$shape, TRUE)
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

      if (is.matrix(tau)) stop("'tau' by location class cannot be used with error ellipses ('smaj.sd', 'smin.sd').")
      data$smaj <- 1; data$smin <- 1; data$eor <- 0; data$x.err <- NA; data$y.err <- NA   # error type 'epar' for dist_args()
      da <- dist_args()
      res <- measurementError_rcpp(data, smaj.sd, smin.sd, eor_range, psi, da$code, da$shape, FALSE,
                                   ell_power, ell[["psi_maj"]], ell[["b_maj"]], ell[["b_min"]], ell_ref[["maj"]], ell_ref[["min"]])

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

      if (ell_power) stop("'psi_maj', 'b_maj' and 'b_min' can only be used with error ellipses.")
      data$x.err <- 1; data$y.err <- 1; data$smaj <- NA; data$smin <- NA; data$eor <- NA   # error type 'sderr'
      tv <- tau_vec()
      da <- dist_args()
      res <- measurementError_LS_rcpp(data, x_sd_val, y_sd_val, tv$x, tv$y, rho_o, da$code, da$shape, FALSE)

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
