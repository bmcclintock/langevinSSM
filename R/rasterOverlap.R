#' Calculate Bhattacharyya's Affinity between two \code{\link[terra]{SpatRaster-class}} objects
#'
#' This function calculates the Bhattacharyya's affinity (also known as the Bhattacharyya coefficient) between two rasters, which is a measure of similarity between two probability distributions. The rasters are first normalized to sum to 1 (to represent probability distributions), and then the affinity is calculated as the sum of the square root of the product of the two distributions across all cells.
#'
#' @param r1 A \code{\link[terra]{SpatRaster-class}} objects object. If any values are negative, \code{rasterOverlap} assumes \code{r1} is on the log scale, and the raster will be exponentiated before calculating the affinity (and a warning will be triggered).
#' @param r2 A \code{\link[terra]{SpatRaster-class}} objects object. If any values are negative, \code{rasterOverlap} assumes \code{r2} is on the log scale, and the raster will be exponentiated before calculating the affinity (and a warning will be triggered).
#' @return A numeric value between 0 (no overlap) and 1 (identical distributions).
#' @importFrom terra compareGeom global
#' @export
rasterOverlap <- function(r1, r2) {

  if (!inherits(r1, "SpatRaster") || !inherits(r2, "SpatRaster")) {
    stop("Both inputs must be terra::SpatRaster objects.")
  }

  if (!terra::compareGeom(r1, r2, stopOnError = FALSE)) {
    stop("Rasters do not have the same geometry (extent, resolution, or CRS). Please resample/project first.")
  }

  min_r1 <- terra::global(r1, "min", na.rm = TRUE)[[1]]
  if (!is.na(min_r1) && min_r1 < 0) {
    warning("Negative values found in r1. Assuming log-scale and exponentiating.")
    r1 <- exp(r1)
  }

  min_r2 <- terra::global(r2, "min", na.rm = TRUE)[[1]]
  if (!is.na(min_r2) && min_r2 < 0) {
    warning("Negative values found in r2. Assuming log-scale and exponentiating.")
    r2 <- exp(r2)
  }

  sum_r1 <- terra::global(r1, "sum", na.rm = TRUE)[[1]]
  sum_r2 <- terra::global(r2, "sum", na.rm = TRUE)[[1]]

  if (is.na(sum_r1) || is.na(sum_r2) || sum_r1 == 0 || sum_r2 == 0) {
    stop("One or both rasters sum to 0 or NA. Cannot normalize into a probability distribution.")
  }

  p <- r1 / sum_r1
  q <- r2 / sum_r2

  pq_sqrt <- sqrt(p * q)
  bc_val <- terra::global(pq_sqrt, "sum", na.rm = TRUE)[[1]]

  return(bc_val)
}
