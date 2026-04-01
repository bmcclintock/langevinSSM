#' Plot Raster
#'
#' Plot a \code{\link[terra]{SpatRaster-class}} object using \pkg{ggplot2} and \pkg{viridis}.
#'
#' @param rast A \code{\link[terra]{SpatRaster-class}} object to plot.
#' @param legend.title Character string for the legend title. Default: \code{NULL}, in which case the name of the raster layer will be used as the legend title.
#' @return A ggplot object containing the raster plot.
#' @examples
#' # exampleCovs included in package; see ?exampleCovs for details
#' UD <- getUD(exampleCovs, beta = c(-4, 6, 5, -0.1) )
#' plotRaster(UD, legend.title = paste0("log(",expression(pi),")"))
#'
#' @importFrom terra crds values res
## #' @importFrom ggplot2 ggplot geom_raster coord_equal aes
## #' @importFrom viridis scale_fill_viridis
#' @export
plotRaster <- function (rast, legend.title = NULL)
{
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for plotting rasters. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package \"viridis\" needed for plotting rasters. Please install it.", call. = FALSE)
  }

  covmap <- data.frame(terra::crds(rast), val = as.numeric(terra::values(rast)))

  p <- ggplot2::ggplot(covmap, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = val)) +
    ggplot2::coord_equal()

  if(is.null(legend.title)) {
    legend.title <- names(rast)
  }
  p <- p + viridis::scale_fill_viridis(name = legend.title)

  return(p)
}
