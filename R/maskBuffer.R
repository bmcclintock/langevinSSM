#' Buffer a binary mask by a specified number of grid cells
#'
#' Expands the restricted areas (0s) into the allowed areas (1s) to create
#' a padding zone along the boundary.
#'
#' @param maskRast A SpatRaster binary mask (1 = allowed, 0 = restricted).
#' @param bufferCells Integer. The number of grid cells to expand the boundary by.
#' @return A SpatRaster binary mask with the padded coastline.
#' @importFrom terra distance
#' @export
maskBuffer <- function(maskRast, bufferCells) {

  if (!inherits(maskRast, "SpatRaster")) stop("'maskRast' must be a SpatRaster")
  if(!is.numeric(bufferCells) || length(bufferCells) != 1 || bufferCells < 0) stop("'bufferCells' must be a single non-negative numeric value.")

  if (bufferCells <= 0) return(maskRast)

  cell_res <- terra::res(maskRast)[1]
  buffer_dist <- bufferCells * cell_res

  message(sprintf("   Buffering mask by %d cells (%g spatial units)...",
                  bufferCells, buffer_dist))

  dist_to_land <- terra::distance(maskRast, target = 1)

  buffered_mask <- terra::ifel(dist_to_land > buffer_dist, 1, 0)
  names(buffered_mask) <- names(maskRast)

  return(buffered_mask)
}
