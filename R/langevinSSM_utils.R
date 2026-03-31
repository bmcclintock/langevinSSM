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
#' A \code{dataLangevin} object containing an example simulated movement track.
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

## modified from aniMotum version 1.2-15 (https://github.com/ianjonsen/aniMotum)
#' @importFrom sf st_as_sf st_crs
#' @importFrom dplyr tibble
format_data <- function (x, id = "id", date = "date", lc = "lc", coord = c("x", "y"), epar = c("smaj", "smin", "eor"), sderr = c("x.sd", "y.sd"), tz = "UTC") {

  if (id %in% names(x))
    stopifnot(`id must be a character string` = is.character(id))
  else {
    stop("An 'id' variable must be included in the input data\n")
  }
  stopifnot(`date must be a character string` = is.character(date))
  stopifnot(`lc must be a character string` = is.character(lc))
  stopifnot(`coord must be a character vector with 1 or 2 elements` = all(is.character(coord)))
  stopifnot(`epar must be a character vector with 3 elements` = all(is.character(epar)))
  stopifnot(`sderr must be a character vector with 2 elements` = all(is.character(sderr)))
  if (all(!inherits(x, "sf"), "geometry" %in% names(x))) {
    if (inherits(x$geometry, "sfc")) {
      x <- sf::st_as_sf(x)
    }
  }
  if (inherits(x, "sf")) {
    coord <- "geometry"
  }
  stopifnot(`An id variable must be included in the input data; \n            see vignette('Overview', package = 'aniMotum')` = id %in%
              names(x))
  stopifnot(`A date/time variable must be included in the input data; \n            see vignette('Overview', package = 'aniMotum')` = date %in%
              names(x))
  stopifnot(`Coordinate variables must be included in the input data; \n            see vignette('Overview', package = 'aniMotum')` = all(coord %in%
                                                                                                                                            names(x)))
  if (inherits(x, "sf") & is.na(sf::st_crs(x))) {
    stop("\nCRS info is missing from input data sf object")
  }
  if (!lc %in% names(x)) {
    if (all(!epar %in% names(x)) & all(sderr %in% names(x))) {
      if (inherits(x, "data.frame", which = TRUE) == 1) {
        x <- data.frame(x, lc = rep("GL", nrow(x)))
      }
      else if (inherits(x, "tbl_df", which = TRUE) == 1) {
        x <- tibble(x, lc = "GL")
      }
      else if (inherits(x, "sf", which = TRUE) == 1) {
        x$lc <- rep("GL", nrow(x))
      }
      x <- x[, c(id, date, "lc", coord, sderr)]
    }
    else if (all(!epar %in% names(x)) & all(!sderr %in% names(x))) {
      message("Guessing that all observations are GPS locations.")
      if (inherits(x, "data.frame", which = TRUE) == 1) {
        x <- data.frame(x, lc = rep("G", nrow(x)))
      }
      else if (inherits(x, "tbl_df", which = TRUE) == 1) {
        x <- dplyr::tibble(x, lc = rep("G", nrow(x)))
      }
      else if (inherits(x, "sf", which = TRUE) == 1) {
        x$lc <- rep("G", nrow(x))
      }
      x <- x[, c(id, date, "lc", coord)]
    }
  }
  xx <- x
  xt.vars <- names(x)[!names(x) %in% c(id, date, lc, coord,
                                       epar, sderr)]
  if (all(!c("lon", "lat") %in% coord, coord != "geometry")) {
    pos1 <- grepl("lon", coord, ignore.case = TRUE)
    pos2 <- grepl("lat", coord, ignore.case = TRUE)
    if (!any(pos1)) {
      pos1 <- grepl("x", coord, ignore.case = TRUE)
      pos2 <- grepl("y", coord, ignore.case = TRUE)
    }
    coord <- coord[c(which(pos1), which(pos2))]
  }
  if (all(!epar %in% names(x), !sderr %in% names(x))) {
    x$smaj <- x$smin <- x$eor <- as.double(NA)
    x$x.sd <- x$y.sd <- as.double(NA)
    xx <- x[, c(id, date, lc, coord, epar, sderr, xt.vars)]
    if (all(!inherits(x, "sf"), all(coord %in% c("lon", "lat")))) {
      names(xx)[1:5] <- c("id", "date", "lc", coord)
      names(xx)[4:5] <- c("lon", "lat")
    }
    else if (all(!inherits(x, "sf"), any(!coord %in% c("lon",
                                                       "lat")))) {
      names(xx)[1:5] <- c("id", "date", "lc", "lon", "lat")
    }
    else if (inherits(x, "sf")) {
      names(xx)[1:4] <- c("id", "date", "lc", coord)
    }
  }
  if (all(epar %in% names(x), !sderr %in% names(x))) {
    x$x.sd <- x$y.sd <- as.double(NA)
    xx <- x[, c(id, date, lc, coord, epar, sderr, xt.vars)]
    if (all(!inherits(x, "sf"), coord != "geometry")) {
      names(xx)[1:8] <- c("id", "date", "lc", coord, "smaj",
                          "smin", "eor")
      names(xx)[4:5] <- c("lon", "lat")
    }
    else if (inherits(x, "sf")) {
      names(xx)[1:7] <- c("id", "date", "lc", coord, "smaj",
                          "smin", "eor")
    }
  }
  if (all(!epar %in% names(x), sderr %in% names(x))) {
    x$smaj <- x$smin <- x$eor <- as.double(NA)
    xx <- x[, c(id, date, lc, coord, epar, sderr, xt.vars)]
    if (!inherits(x, "sf")) {
      names(xx)[c(1:5, 9:10)] <- c("id", "date", "lc",
                                   coord, "x.sd", "y.sd")
    }
    else if (inherits(x, "sf")) {
      names(xx)[c(1:4, 8:9)] <- c("id", "date", "lc", coord,
                                  "x.sd", "y.sd")
    }
  }
  if (is.factor(xx$id))
    xx$id <- droplevels(xx$id)
  xx$id <- as.character(xx$id)
  if (!inherits(xx$date, "POSIXt")) {
    xx$date <- try(as.POSIXct(xx$date, tz = tz), silent = TRUE)
    if (inherits(xx$date, "try-error"))
      stop("dates must be in a standard format: YYYY-MM-DD HH:MM:SS")
  }
  xx <- xx[order(xx$date), ]
  return(xx)
}
