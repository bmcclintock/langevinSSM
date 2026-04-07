#' Plot Raster
#'
#' Plot a \code{\link[terra]{SpatRaster-class}} object using \pkg{ggplot2} and \pkg{viridis}.
#' Supports multi-layer rasters by automatically creating a faceted multi-panel plot.
#'
#' @param rast A \code{\link[terra]{SpatRaster-class}} object to plot.
#' @param legend.title Character string for the legend title. Default: \code{NULL}, in which case the name of the first raster layer will be used as the legend title.
#' @return A ggplot object containing the raster plot.
#' @examples
#' # exampleCovs included in package; see ?exampleCovs for details
#' UD <- getUD(exampleCovs, beta = c(-4, 6, 5, -0.1) )
#' plotRaster(UD, legend.title = paste0("log(",expression(pi),")"))
#'
## #' @importFrom ggplot2 ggplot geom_raster coord_equal aes
## #' @importFrom viridis scale_fill_viridis
#' @importFrom terra as.data.frame nlyr time
#' @export
plotRaster <- function (rast, legend.title = NULL)
{
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for plotting rasters. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package \"viridis\" needed for plotting rasters. Please install it.", call. = FALSE)
  }

  if(is.null(legend.title)) {
    legend.title <- names(rast)[1]
  }

  n_layers <- terra::nlyr(rast)
  layer_times <- terra::time(rast)

  if (n_layers > 1) {
    if (!is.null(layer_times) && !all(is.na(layer_times))) {
      # Use time attributes for facet labels if available
      new_names <- paste0("Time: ", layer_times)
    } else {
      # Fallback if no time attributes exist
      new_names <- paste0("Layer ", seq_len(n_layers))
    }
    names(rast) <- new_names
  }

  covmap_wide <- terra::as.data.frame(rast, xy = TRUE, na.rm = TRUE)
  layer_names <- setdiff(names(covmap_wide), c("x", "y"))

  covmap_long <- do.call(rbind, lapply(layer_names, function(lyr) {
    data.frame(x = covmap_wide$x, y = covmap_wide$y, layer = lyr, val = covmap_wide[[lyr]])
  }))

  # Lock the factor levels to ensure chronological sorting in the plot
  covmap_long$layer <- factor(covmap_long$layer, levels = layer_names)

  p <- ggplot2::ggplot(covmap_long, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = val)) +
    ggplot2::coord_equal()

  if (length(layer_names) > 1) {
    p <- p + ggplot2::facet_wrap(~ layer)
  }

  p <- p + viridis::scale_fill_viridis(name = legend.title)

  return(p)
}
