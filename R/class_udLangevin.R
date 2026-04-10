#' Constructor and validator for udLangevin class
#'
#' @param x A list containing the utilization distribution raster(s).
#' @return A validated \code{udLangevin} object.
#' @noRd
class_udLangevin <- function(x) {
  if (!is.list(x)) stop("Object must be a list to be a 'udLangevin' object.")

  if (is.null(x$UD)) stop("'udLangevin' object must contain a 'UD' element.")
  if (!inherits(x$UD, "SpatRaster")) stop("The 'UD' element must be a SpatRaster.")

  class(x) <- unique(c("udLangevin", class(x)))

  return(x)
}
