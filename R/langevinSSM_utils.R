#' @importFrom utils globalVariables
utils::globalVariables(c("mu.x", "mu.y", "vel.x", "vel.y", "id", "psi", "tau", "lon", "lat", "dt", "x", "y", "smaj", "smin", "eor", "x.sd", "y.sd", "val"))

#' Example Spatial Covariates
#'
#' A list of SpatRaster objects used as example covariates.
#' @name exampleCovs
#' @docType data
#' @export
NULL

.onLoad <- function(libname, pkgname) {
  env <- asNamespace(pkgname)

  makeActiveBinding("exampleCovs", function() {

    cov_names <- c("exampleCov1", "exampleCov2", "exampleCov3", "exampleCov4")

    cov_list <- lapply(cov_names, function(cov) {

      file_path <- system.file("extdata", paste0(cov, ".tif"), package = "langevinSSM")

      # Safety check
      if (file_path == "") {
        stop(paste("Could not find", paste0(cov, ".tif"), "in the extdata/ folder."))
      }

      terra::rast(file_path)
    })

    names(cov_list) <- c("cov1","cov2","cov3","d2c")

    return(cov_list)

  }, env)
}

#' Example Tracking Data
#'
#' A \code{dataLangevin} object containing example simulated movement tracks.
#'
#' @name exampleDat
#' @docType data
NULL

rasterList <- function (rast)
{
  lim <- as.vector(terra::ext(rast))
  res <- terra::res(rast)
  xgrid <- seq(lim[1] + res[1]/2, lim[2] - res[1]/2, by = res[1])
  ygrid <- seq(lim[3] + res[2]/2, lim[4] - res[2]/2, by = res[2])

  # Add wide = TRUE so terra returns [nrow, ncol] instead of [ncells, nlayers]
  z_mat <- terra::as.matrix(rast, wide = TRUE)
  z <- t(apply(z_mat, 2, rev))

  return(list(x = xgrid, y = ygrid, z = z))
}

## modified from aniMotum version 1.2-15
#' @importFrom sf st_as_sf st_crs
#' @importFrom dplyr tibble
format_data <- function(x, id = "id", date = "date", lc = "lc", coord = c("x", "y"), epar = c("smaj", "smin", "eor"), sderr = c("x.sd", "y.sd"), tz = "UTC") {

  if (id %in% names(x))
    stopifnot(`id must be a character string` = is.character(id))
  else stop("An 'id' variable must be included in the input data\n")

  stopifnot(`date must be a character string` = is.character(date))
  stopifnot(`lc must be a character string` = is.character(lc))
  stopifnot(`coord must be a character vector with 2 elements` = all(is.character(coord)) && length(coord) == 2)
  stopifnot(`epar must be a character vector with 3 elements` = all(is.character(epar)) && length(epar) == 3)
  stopifnot(`sderr must be a character vector with 2 elements` = all(is.character(sderr)) && length(sderr) == 2)

  stopifnot(`An id variable must be included in the input data` = id %in% names(x))
  stopifnot(`A date/time variable must be included in the input data` = date %in% names(x))
  stopifnot(`Coordinate variables must be included in the input data` = all(coord %in% names(x)))

  # --- DYNAMIC COLUMN INJECTION ---
  if (all(!epar %in% names(x))) {
    x[[epar[1]]] <- as.double(NA)
    x[[epar[2]]] <- as.double(NA)
    x[[epar[3]]] <- as.double(NA)
  }
  if (all(!sderr %in% names(x))) {
    x[[sderr[1]]] <- as.double(NA)
    x[[sderr[2]]] <- as.double(NA)
  }

  xt.vars <- names(x)[!names(x) %in% c(id, date, lc, coord, epar, sderr)]

  # Subset cleanly
  xx <- x[, c(id, date, lc, coord, epar, sderr, xt.vars)]

  # Force standard names
  new_names <- c("id", "date", "lc", coord, "smaj", "smin", "eor", "x.sd", "y.sd", xt.vars)
  names(xx) <- new_names

  if (is.factor(xx$id)) xx$id <- droplevels(xx$id)
  xx$id <- as.character(xx$id)

  if (!inherits(xx$date, "POSIXt")) {
    xx$date <- try(as.POSIXct(xx$date, tz = tz), silent = TRUE)
    if (inherits(xx$date, "try-error"))
      stop("dates must be in a standard format: YYYY-MM-DD HH:MM:SS")
  }

  xx <- xx[order(xx$date), ]
  return(xx)
}

# modfied from aniMotum version 1.2-15
#' Error multiplication factors
#'
#' A function to generate a data frame of error multiplication factors (EMF) for different location classes, which can be used to account for measurement error for observations that lack error information but have a known location quality class. The default values are based on the EMF values for Argos satellite telemetry data, but users can specify their own EMF values for different location classes as needed. It is a modified version of the \code{emf} function from the \href{https://ianjonsen.github.io/aniMotum/}{aniMotum} package.
#'
#' @param gps A numeric value or a vector of length 2 specifying the error multiplication factor for GPS locations. If a single value is provided, it will be used for both x and y axes. Default is 0.1 (i.e. GPS errors are 10x more accurate than Argos \code{lc} 3.
#' @param emf.x A numeric vector of length 6 specifying the error multiplication factors for the x-axis for each location class (in order: 3, 2, 1, 0, A, B, where Z is assumed equal to B). Default values are based on the EMF values for Argos satellite telemetry data.
#' @param emf.y A numeric vector of length 6 specifying the error multiplication factors for the y-axis for each location class (in order: 3, 2, 1, 0, A, B, where Z is assumed equal to B). Default values are based on the EMF values for Argos satellite telemetry data.
#' @return A data frame with columns \code{lc}, \code{x.sd}, and \code{y.sd} containing the error multiplication factors for each location class. The location classes included are "G" for GPS and "3", "2", "1", "0", "A", "B", and "Z" for Argos satellite telemetry data.
#' @export
get_emf <- function (gps = 0.1, emf.x = c(1, 1.54, 3.72, 13.51, 23.9, 44.22),
                     emf.y = c(1, 1.29, 2.55, 14.99, 22, 32.53))
{
  if (!length(gps) %in% 1:2)
    stop("GPS emf must be a vector of length 1 or 2")
  if (length(emf.x) != 6)
    stop("Argos emf.x must be a vector of length 6")
  if (length(emf.y) != 6)
    stop("Argos emf.y must be a vector of length 6")
  if (length(gps) == 1)
    gps <- c(gps, gps)
  data.frame(emf.x = c(gps[1], emf.x, emf.x[6]), emf.y = c(gps[2],
                                                           emf.y, emf.y[6]), lc = as.character(c("G", "3", "2",
                                                                                                 "1", "0", "A", "B", "Z")))
}

checkErrorData <- function(data, coord=c("x","y")){
  if(any(is.na(data[,coord[1]]) & !is.na(data[,coord[2]])) | any(!is.na(data[,coord[1]]) & is.na(data[,coord[2]]))) stop("Missing values (NA) in coordinates must be in both x and y columns.")
  if(any(is.na(data[,coord[1]]) & (!is.na(data$smaj) | !is.na(data$smin) | !is.na(data$eor) | !is.na(data$x.sd) | !is.na(data$y.sd)))) stop("Measurement error terms must be NA when there are missing values in the coordinates.")
  if(any(is.na(data$smaj) & (!is.na(data$smin) & !is.na(data$eor)))) stop("When using the error ellipse model, smaj, smin, and eor must all be provided or all be NA.")
  if(any(is.na(data$x.sd) & !is.na(data$y.sd))) stop("When using the x- and y-axis error model, x.sd and y.sd must both be provided or both be NA.")
  if(any((!is.na(data$smaj) & !is.na(data$smin) & !is.na(data$eor)) & (!is.na(data$x.sd) | !is.na(data$y.sd)))) stop("Cannot provide both error ellipse and x- and y-axis error terms. If using the error ellipse, 'smaj', 'smin', and 'eor' must all be provided and 'x.sd' and 'y.sd' must both be NA. If using the x- and y-axis error model, 'x.sd' and 'y.sd' must both be provided and 'smaj', 'smin', and 'eor' must all be NA.")
  if(any((!is.na(data$smaj) | !is.na(data$smin) | !is.na(data$eor)) & (!is.na(data$x.sd) & !is.na(data$y.sd)))) stop("Cannot provide both error ellipse and x- and y-axis error terms. If using the error ellipse, 'smaj', 'smin', and 'eor' must all be provided and 'x.sd' and 'y.sd' must both be NA. If using the x- and y-axis error model, 'x.sd' and 'y.sd' must both be provided and 'smaj', 'smin', and 'eor' must all be NA.")
  if(isTRUE(any(data$eor<0 | data$eor > pi))) stop("Error ellipse orientation (eor) must be between 0 and pi radians.")
}
