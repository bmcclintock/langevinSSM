#' Plot a fitLangevin object
#'
#' Plots the estimated Utilization Distribution (UD) and overlays the estimated true locations (mu) using \code{ggplot2}. If the original data is provided, the observed locations are also plotted for comparison.
#'
#' @param x A \code{fitLangevin} object.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects used to fit the model.
#' @param log Logical indicating whether to plot the log UD (default) or the probability UD.
#' @param extent Optional. A numeric vector of length 4 \code{c(xmin, xmax, ymin, ymax)} or a \code{\link[terra]{SpatExtent}} object defining the bounding box to zoom the plot. If \code{NULL} (default), the extent is automatically calculated from the track data with a 10\% buffer.
#' @param data Optional \code{dataLangevin} object. If provided, the observed coordinates will be plotted beneath the estimated locations.
#' @param time Optional. Indicates which layer(s) of a dynamic UD to plot. Can be a numeric index, a layer name, or a \code{POSIXct}/\code{Date} object matching the raster's \code{\link[terra]{time}} attribute. If \code{NULL} (default), all UD layers are plotted.
#' @param compact Logical indicating whether to plot all tracks on a single panel (\code{TRUE}, default) or plot each track separately (\code{FALSE}).
#' @param ... Additional arguments (currently ignored, kept for S3 compatibility).
#'
#' @return A \code{\link[ggplot2]{ggplot}} object containing the plot (or a list of \code{\link[ggplot2]{ggplot}} objects if \code{compact = FALSE}), which can be further modified by the user.
#' @examples
#' \dontrun{
#' # exampleDat included in package; see ?exampleDat for details
#' # exampleCovs included in package; see ?exampleCovs for details
#' fit <- fitLangevin(exampleDat, spatialCovs = exampleCovs, silent=TRUE)
#' p <- plot(fit, spatialCovs = exampleCovs, data = exampleDat)
#' p
#' p + ggplot2::labs(title = "New Plot Title", subtitle = "Optional Subtitle")
#' }
#'
#' @method plot fitLangevin
#' @importFrom terra as.data.frame nlyr time crop ext
#' @export
plot.fitLangevin <- function(x, spatialCovs, log = TRUE, extent = NULL, data = NULL, time = NULL, compact = TRUE, ...) {

  if (missing(spatialCovs)) {
    stop("You must provide the 'spatialCovs' list to compute and plot the estimated utilization distribution.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for plotting rasters. Please install it.", call. = FALSE)
  }

  verify_signatures(x, data = data, spatialCovs = spatialCovs)

  track_id <- x$estimates$random$mu$est$id
  if (is.null(track_id)) track_id <- rep("1", nrow(x$estimates$random$mu$est))

  track_df <- data.frame(
    x = x$estimates$random$mu$est[, "mu.x"],
    y = x$estimates$random$mu$est[, "mu.y"],
    id = as.character(track_id),
    type = "Estimated"
  )

  if (!is.null(data)) {
    if (inherits(data, "dataLangevin")) {
      obs_df <- data.frame(
        x = data$x,
        y = data$y,
        id = as.character(data$id),
        type = "Observed"
      )
      track_df <- rbind(obs_df, track_df)
    } else {
      warning("'data' is not a dataLangevin object. Skipping observed locations.")
    }
  }

  track_df$type <- factor(track_df$type, levels = c("Observed", "Estimated"))

  rn <- rownames(x$estimates$natural)
  beta_idx <- which(grepl("^beta", rn))
  beta_est <- x$estimates$natural[beta_idx, "Estimate"]

  ud_full <- getUD(spatialCovs = spatialCovs, beta = beta_est, log = log)

  if (!is.null(time)) {
    if (terra::nlyr(ud_full) > 1) {
      ud_times <- terra::time(ud_full)
      if (!is.null(ud_times) && any(ud_times %in% time)) {
        ud_full <- ud_full[[which(ud_times %in% time)]]
      } else {
        ud_full <- tryCatch({
          ud_full[[time]]
        }, error = function(e) {
          stop("The provided 'time' argument could not be matched to a layer index, name, or time attribute in 'spatialCovs'.")
        })
      }
    } else {
      warning("The 'time' argument was provided, but the utilization distribution is static.")
    }
  }

  plot_ids <- if (compact) "all" else as.character(unique(track_df$id))
  plot_list <- list()

  for (pid in plot_ids) {
    if (compact) {
      trk_sub <- track_df
    } else {
      trk_sub <- track_df[track_df$id == pid, ]
    }

    current_extent <- extent
    if (is.null(current_extent)) {
      x_min <- min(trk_sub$x, na.rm = TRUE)
      x_max <- max(trk_sub$x, na.rm = TRUE)
      y_min <- min(trk_sub$y, na.rm = TRUE)
      y_max <- max(trk_sub$y, na.rm = TRUE)

      x_range <- x_max - x_min
      y_range <- y_max - y_min
      max_range <- max(x_range, y_range)

      if (max_range == 0) max_range <- 1

      x_mid <- (x_max + x_min) / 2
      y_mid <- (y_max + y_min) / 2

      current_extent <- c(
        x_mid - 0.6 * max_range,
        x_mid + 0.6 * max_range,
        y_mid - 0.6 * max_range,
        y_mid + 0.6 * max_range
      )
    }

    xlim <- NULL
    ylim <- NULL
    ud_sub <- ud_full

    if (!is.null(current_extent)) {
      crop_ext <- tryCatch(terra::ext(current_extent), error = function(e) NULL)

      if (!is.null(crop_ext)) {
        ext_vec <- as.vector(crop_ext)
        xlim <- c(ext_vec["xmin"], ext_vec["xmax"])
        ylim <- c(ext_vec["ymin"], ext_vec["ymax"])

        ud_sub <- tryCatch(terra::crop(ud_sub, crop_ext), error = function(e) {
          warning("Failed to crop raster to extent. The calculated/provided extent may not overlap with the raster.")
          return(ud_sub)
        })
      } else {
        warning("'extent' must be a numeric vector of length 4 (xmin, xmax, ymin, ymax) or a SpatExtent object (see ?terra::ext). Ignoring 'extent'.")
      }
    }

    ud_df <- terra::as.data.frame(ud_sub, xy = TRUE, na.rm = TRUE)
    layer_names <- setdiff(names(ud_df), c("x", "y"))

    ud_long <- do.call(rbind, lapply(layer_names, function(lyr) {
      data.frame("x" = ud_df$x, "y" = ud_df$y, "time_layer" = lyr, "UD" = ud_df[[lyr]])
    }))

    title_text <- if (compact) "Estimated utilization distribution and tracks" else paste("Estimated utilization distribution and tracks - ID:", pid)

    p <- ggplot2::ggplot() +
      ggplot2::geom_raster(data = ud_long, ggplot2::aes(x = x, y = y, fill = UD)) +
      ggplot2::scale_fill_viridis_c(name = ifelse(log, expression(log(pi)), expression(pi)), option = "viridis", na.value = "transparent") +
      ggplot2::geom_path(data = trk_sub, ggplot2::aes(x = x, y = y, color = type, linetype = type, group = interaction(id, type)), linewidth = 0.5, alpha = 0.8) +
      ggplot2::geom_point(data = trk_sub, ggplot2::aes(x = x, y = y, color = type, shape = type), size = 1, alpha = 0.8) +
      ggplot2::scale_color_manual(name = "Track", values = c("Observed" = "lightgrey", "Estimated" = "tomato")) +
      ggplot2::scale_linetype_manual(name = "Track", values = c("Observed" = "dashed", "Estimated" = "solid")) +
      ggplot2::scale_shape_manual(name = "Track", values = c("Observed" = 16, "Estimated" = 16)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "easting (x)", y = "northing (y)", title = title_text)

    if (!is.null(xlim) && !is.null(ylim)) {
      p <- p + ggplot2::coord_equal(xlim = xlim, ylim = ylim)
    } else {
      p <- p + ggplot2::coord_equal()
    }

    if (length(layer_names) > 1) {
      p <- p + ggplot2::facet_wrap(~ time_layer)
    }

    plot_list[[pid]] <- p
  }

  if (compact || length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }
}
