#' Prepare a spatial barrier for Langevin models
#'
#' Converts a binary spatial mask into a signed distance field (SDF) required for barrier constraints in the \code{langevinSSM} package. The resulting raster is tagged with a special attribute allowing downstream functions (\code{\link{fitLangevin}}, \code{\link{simLangevin}}, \code{\link{getUD}}) to automatically detect and apply the barrier penalty.
#'
#' @param mask_rast A \code{\link[terra]{SpatRaster-class}} object containing a binary mask. Values of \code{1} indicate allowed movement areas (e.g., water), and values of \code{0} indicate restricted areas (e.g., land).
#'
#' @return A \code{\link[terra]{SpatRaster-class}} object containing the calculated signed distance field. Positive values represent distance into the allowed area; negative values represent distance into the restricted area. The object is assigned a \code{"barLangevin" = TRUE} attribute.
#' @export
#' @importFrom stats na.omit
#' @importFrom terra distance unique nlyr
prepBarrier <- function(mask_rast) {

  if (!inherits(mask_rast, "SpatRaster")) {
    stop("'mask_rast' must be a terra::SpatRaster object.")
  }

  if (terra::nlyr(mask_rast) > 1) {
    stop("Barrier rasters must be static (single layer).")
  }

  unq_vals <- stats::na.omit(terra::unique(mask_rast)[[1]])
  if (!all(unq_vals %in% c(0, 1))) {
    stop("Barrier raster must only contain binary values (0 and 1).")
  }

  message("   Calculating signed distance field for barrier...")

  # Positive allowed, negative prohibited
  sdf <- terra::distance(mask_rast, target = 1) - terra::distance(mask_rast, target = 0)

  # Tag the SpatRaster with our custom attribute
  attr(sdf, "barLangevin") <- TRUE

  return(sdf)
}
