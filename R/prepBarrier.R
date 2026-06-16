#' Prepare a spatial barrier for Langevin models
#'
#' Converts a binary spatial mask into a signed distance field (SDF) required for barrier constraints in the \code{langevinSSM} package. The resulting raster is tagged with a special attribute allowing downstream functions (\code{\link{fitLangevin}}, \code{\link{simLangevin}}, \code{\link{getUD}}) to automatically detect and apply the barrier penalty.
#'
#' @param maskRast A \code{\link[terra]{SpatRaster-class}} object containing a binary mask. Values of \code{1} indicate allowed movement areas (e.g., water), and values of \code{0} indicate restricted areas (e.g., land).
#'
#' @return A \code{\link[terra]{SpatRaster-class}} object containing the calculated signed distance field. Positive values represent distance into the allowed area; negative values represent distance into the restricted area. The object is assigned a \code{"barLangevin" = TRUE} attribute.
#' @export
#' @importFrom stats na.omit
#' @importFrom terra distance unique nlyr crs "crs<-"
prepBarrier <- function(maskRast) {

  if (!inherits(maskRast, "SpatRaster")) {
    stop("'maskRast' must be a terra::SpatRaster object.")
  }

  if (terra::nlyr(maskRast) > 1) {
    stop("Barrier rasters must be static (single layer).")
  }

  unq_vals <- stats::na.omit(terra::unique(maskRast)[[1]])
  if (!all(unq_vals %in% c(0, 1))) {
    stop("Barrier raster must only contain binary values (0 and 1).")
  }

  orig_crs <- terra::crs(maskRast)

  # If no CRS is provided (e.g., simulated data), temporarily assign "local"
  # to prevent terra's internal C++ warning about unknown CRSs
  if (orig_crs == "") {
    terra::crs(maskRast) <- "local"
  }

  message("   Calculating signed distance field for barrier...")

  # Positive allowed, negative prohibited
  sdf <- terra::distance(maskRast, target = 1) - terra::distance(maskRast, target = 0)

  terra::crs(sdf) <- orig_crs

  # Tag the SpatRaster with our custom attribute
  attr(sdf, "barLangevin") <- TRUE

  return(sdf)
}
