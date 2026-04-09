#' Plot Raster
#'
#' Plot a \code{\link[terra]{SpatRaster-class}} object using \pkg{ggplot2}.
#' Supports multi-layer rasters by automatically creating a faceted multi-panel plot.
#'
#' @param rast A \code{\link[terra]{SpatRaster-class}} object to plot.
#' @param legend.title Character string for the legend title. Default: \code{NULL}, in which case the name of the first raster layer will be used as the legend title.
#' @param extent Optional. A numeric vector of length 4 \code{c(xmin, xmax, ymin, ymax)} or a \code{\link[terra]{SpatExtent}} object defining the bounding box. Default: \code{NULL}.
#' @param time Optional. Indicates which layer(s) of a dynamic raster to plot. Can be a numeric index, a layer name, or a value matching the raster's \code{\link[terra]{time}} attribute (e.g., a \code{POSIXct}/\code{Date} or numeric object). If \code{NULL} (default), all layers are plotted.
#' @return A \code{\link[ggplot2]{ggplot}} object containing the raster plot.
#' @examples
#' # exampleCovs included in package; see ?exampleCovs for details
#' UD <- getUD(exampleCovs, beta = c(-4, 6, 5, -0.1) )
#' plotRaster(UD, legend.title = paste0("log(",expression(pi),")"))
#'
#' @importFrom terra as.data.frame nlyr time ext
#' @export
plotRaster <- function(rast, legend.title = NULL, extent = NULL, time = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for plotting rasters. Please install it.", call. = FALSE)
  }

  if(is.null(legend.title)) {
    legend.title <- names(rast)[1]
  }

  # --- Time Subsetting ---
  if (!is.null(time)) {
    if (terra::nlyr(rast) > 1) {
      rast_times <- terra::time(rast)
      if (!is.null(rast_times) && any(rast_times %in% time)) {
        rast <- rast[[which(rast_times %in% time)]]
      } else {
        rast <- tryCatch(rast[[time]], error = function(e) stop("Invalid 'time' argument. Could not match to a layer index, name, or time attribute."))
      }
    }
  }

  # --- Extent Parsing ---
  xlim <- ylim <- NULL
  if (!is.null(extent)) {
    crop_ext <- tryCatch(terra::ext(extent), error = function(e) NULL)
    if (!is.null(crop_ext)) {
      ext_vec <- as.vector(crop_ext)
      xlim <- c(ext_vec["xmin"], ext_vec["xmax"])
      ylim <- c(ext_vec["ymin"], ext_vec["ymax"])
    } else {
      warning("'extent' must be a numeric vector of length 4 (xmin, xmax, ymin, ymax) or a SpatExtent object. Ignoring 'extent'.")
    }
  }

  n_layers <- terra::nlyr(rast)
  layer_times <- terra::time(rast)

  if (n_layers > 1) {
    if (!is.null(layer_times) && !all(is.na(layer_times))) {
      # Use time attributes for facet labels if available
      new_names <- paste0("Time: ", layer_times)
    } else {
      new_names <- paste0("Layer ", seq_len(n_layers))
    }
    names(rast) <- new_names
  }

  covmap_wide <- terra::as.data.frame(rast, xy = TRUE, na.rm = TRUE)
  layer_names <- setdiff(names(covmap_wide), c("x", "y"))

  covmap_long <- do.call(rbind, lapply(layer_names, function(lyr) {
    data.frame(x = covmap_wide$x, y = covmap_wide$y, layer = lyr, val = covmap_wide[[lyr]])
  }))

  # Lock the factor levels to ensure chronological sorting in the plot facets
  covmap_long$layer <- factor(covmap_long$layer, levels = layer_names)

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = covmap_long, ggplot2::aes(x = x, y = y, fill = val)) +
    ggplot2::scale_fill_viridis_c(name = legend.title, option = "viridis", na.value = "transparent")

  # Apply the coordinates once here to avoid ggplot2 overwriting warnings
  if (!is.null(xlim) && !is.null(ylim)) {
    p <- p + ggplot2::coord_equal(xlim = xlim, ylim = ylim)
  } else {
    p <- p + ggplot2::coord_equal()
  }

  if (length(layer_names) > 1) {
    p <- p + ggplot2::facet_wrap(~ layer)
  }

  return(p)
}
