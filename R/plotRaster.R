#' Plot Raster
#'
#' Plot a \code{\link[terra]{SpatRaster-class}} object using \pkg{ggplot2}.
#' Supports multi-layer rasters by automatically creating a faceted multi-panel plot.
#'
#' @param rast A \code{\link[terra]{SpatRaster-class}} object to plot.
#' @param legend.title Character string for the legend title. Default: \code{NULL}, in which case the name of the first raster layer (see \code{\link[terra]{names}}) will be used as the legend title.
#' @param extent Optional. A numeric vector of length 4 \code{c(xmin, xmax, ymin, ymax)} or a \code{\link[terra]{SpatExtent}} object defining the bounding box. Default: \code{NULL}.
#' @param time Optional. Indicates which layer(s) of a dynamic raster to plot. Can be a numeric index, a layer name, or a value matching the raster's \code{\link[terra]{time}} attribute (e.g., a \code{POSIXct}/\code{Date} or numeric object). If \code{NULL} (default), all layers are plotted.
#' @param maskRast \code{\link[terra]{SpatRaster-class}} object for areas to be masked out (set to \code{0}) before plotting the raster. Default: \code{NULL} (no mask).
#' @param ... Additional arguments passed to \code{\link[ggplot2]{scale_fill_viridis_c}} (e.g., \code{direction = -1}).
#' @return A \code{\link[ggplot2]{ggplot}} object containing the raster plot.
#' @examples
#' # exampleCovs included in package; see ?exampleCovs for details
#' plotRaster(exampleCovs$cov1)
#'
#' @importFrom terra as.data.frame nlyr time ext crop ifel compareGeom
#' @export
plotRaster <- function(rast, legend.title = NULL, extent = NULL, time = NULL, maskRast = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for plotting rasters. Please install it.", call. = FALSE)
  }

  orig_name <- names(rast)[1]
  if(is.null(legend.title)) legend.title <- orig_name

  # --- Time Subsetting ---
  if (!is.null(time)) {
    if (terra::nlyr(rast) > 1) {
      rast_times <- terra::time(rast)
      subset_idx <- if (!is.null(rast_times) && any(rast_times %in% time)) {
        which(rast_times %in% time)
      } else {
        tryCatch({ seq_len(terra::nlyr(rast))[time] },
                 error = function(e) stop("Invalid 'time' argument."))
      }
      if(any(is.na(subset_idx))) stop("Invalid 'time' argument.")

      rast <- rast[[subset_idx]]

      if (terra::nlyr(rast) == 1) {
        t_val <- terra::time(rast)
        names(rast) <- if (!is.null(t_val) && !is.na(t_val)) paste0("Time: ", t_val) else paste0("Layer ", subset_idx)
      }
    }
  }

  # --- Extent Parsing ---
  xlim <- ylim <- NULL
  crop_ext <- NULL
  if (!is.null(extent)) {
    crop_ext <- tryCatch(terra::ext(extent), error = function(e) NULL)
    if (!is.null(crop_ext)) {
      ext_vec <- as.vector(crop_ext)
      xlim <- c(ext_vec["xmin"], ext_vec["xmax"]); ylim <- c(ext_vec["ymin"], ext_vec["ymax"])
      rast <- tryCatch({ terra::crop(rast, crop_ext) },
                       error = function(e) { warning("terra::crop failed. Ignoring extent."); return(rast) })
    } else { warning("Invalid extent object. Ignoring.") }
  }

  # --- Masking ---
  if (!is.null(maskRast)) {

    if(!inherits(maskRast, "SpatRaster")) stop("'maskRast' must be a SpatRaster")
    if (!terra::compareGeom(rast, maskRast, stopOnError = FALSE)) stop("The 'maskRast' raster must share the same projection (CRS), extent, and resolution as the rasters in 'spatialCovs'.")

    maskRast_eval <- tryCatch({
      if (!is.null(crop_ext)) terra::crop(maskRast, crop_ext) else maskRast
    }, error = function(e) {
      warning("terra::crop failed on maskRast. Ignoring maskRast.")
      return(NULL)
    })
    if(!is.null(maskRast_eval)) {
      rast <- terra::ifel(maskRast_eval <= 0, NA, rast)
    }
  }

  # --- Renaming for Facets ---
  n_layers <- terra::nlyr(rast)
  layer_times <- terra::time(rast)
  if (n_layers > 1) {
    names(rast) <- if (!is.null(layer_times) && !all(is.na(layer_times))) paste0("Time: ", layer_times) else paste0("Layer ", seq_len(n_layers))
  }

  covmap_wide <- terra::as.data.frame(rast, xy = TRUE, na.rm = FALSE)

  layer_names <- setdiff(names(covmap_wide), c("x", "y"))
  covmap_long <- do.call(rbind, lapply(layer_names, function(lyr) {
    data.frame(x = covmap_wide$x, y = covmap_wide$y, layer = lyr, val = covmap_wide[[lyr]])
  }))
  covmap_long$layer <- factor(covmap_long$layer, levels = layer_names)

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = covmap_long, ggplot2::aes(x = x, y = y, fill = val)) +
    ggplot2::scale_fill_viridis_c(name = legend.title, option = "viridis", na.value = "transparent", ...) +
    ggplot2::theme_minimal()

  p <- p + (if (!is.null(xlim)) ggplot2::coord_equal(xlim = xlim, ylim = ylim) else ggplot2::coord_equal())
  if (length(layer_names) > 1) p <- p + ggplot2::facet_wrap(~ layer)

  return(p)
}
