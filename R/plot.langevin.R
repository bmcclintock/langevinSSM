#' Plot Langevin model fits and data
#'
#' These methods plot objects associated with the \code{langevinSSM} package using \code{ggplot2}.
#' \itemize{
#'   \item \strong{\code{plot.fitLangevin}:} Plots the estimated Utilization Distribution (UD) and overlays the estimated true locations (mu). If the original data is provided, observed locations are also plotted.
#'   \item \strong{\code{plot.dataLangevin}:} Plots the spatial covariates and overlays the observed locations. If the data is a simulated \code{simLangevin} object, both true (latent) and observed locations are plotted.
#'   \item \strong{\code{plot.simLangevin}:} If \code{beta} is provided, plots the theoretical UD based on those coefficients and overlays true (latent) and observed locations. If \code{beta} is omitted, defaults to \code{plot.dataLangevin} behavior (plotting individual covariates).
#' }
#'
#' @param x A \code{fitLangevin}, \code{dataLangevin}, or \code{simLangevin} object.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects. Used to compute the UD or plotted as the background.
#' @param beta Optional numeric vector of habitat selection coefficients (for \code{simLangevin} only). Must match the length of \code{spatialCovs}. If provided, plots the UD instead of individual covariates.
#' @param log Logical (for \code{fitLangevin} and \code{simLangevin}). Indicates whether to plot the log UD (\code{TRUE}, default) or the probability UD.
#' @param extent Optional. A numeric vector of length 4 \code{c(xmin, xmax, ymin, ymax)} or a \code{\link[terra]{SpatExtent}} object defining the bounding box. If \code{NULL} (default), the extent is automatically calculated from the track data.
#' @param data Optional \code{dataLangevin} object (for \code{fitLangevin} only). If provided, observed coordinates will be plotted beneath the estimated locations.
#' @param time Optional (for \code{fitLangevin} only). Indicates which layer(s) of a dynamic UD to plot. Can be a numeric index, a layer name, or a \code{POSIXct}/\code{Date} object. If \code{NULL} (default), all UD layers are plotted.
#' @param compact Logical indicating whether to plot all tracks on a single panel (\code{TRUE}, default) or plot each track separately (\code{FALSE}).
#' @param ... Additional arguments (currently ignored, kept for S3 compatibility).
#'
#' @return A \code{\link[ggplot2]{ggplot}} object, or a list of \code{ggplot} objects depending on the input type and \code{compact} argument.
#'
#' @name plot.langevin
#' @importFrom terra as.data.frame nlyr time crop ext
NULL

#' @rdname plot.langevin
#' @method plot fitLangevin
#' @export
plot.fitLangevin <- function(x, spatialCovs, log = TRUE, extent = NULL, data = NULL, time = NULL, compact = TRUE, ...) {
  if (missing(spatialCovs)) stop("You must provide the 'spatialCovs' list.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package \"ggplot2\" needed for plotting. Please install it.", call. = FALSE)

  verify_signatures(x, data = data, spatialCovs = spatialCovs)

  # --- Prepare Tracks ---
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
      obs_df <- data.frame(x = data$x, y = data$y, id = as.character(data$id), type = "Observed")
      track_df <- rbind(obs_df, track_df)
    } else {
      warning("'data' is not a dataLangevin object. Skipping observed locations.")
    }
  }
  track_df$type <- factor(track_df$type, levels = c("Observed", "Estimated"))

  # --- Prepare UD Raster ---
  rn <- rownames(x$estimates$natural)
  beta_est <- x$estimates$natural[which(grepl("^beta", rn)), "Estimate"]
  ud_full <- getUD(spatialCovs = spatialCovs, beta = beta_est, log = log)

  if (!is.null(time)) {
    if (terra::nlyr(ud_full) > 1) {
      ud_times <- terra::time(ud_full)
      if (!is.null(ud_times) && any(ud_times %in% time)) {
        ud_full <- ud_full[[which(ud_times %in% time)]]
      } else {
        ud_full <- tryCatch(ud_full[[time]], error = function(e) stop("Invalid 'time' argument."))
      }
    } else {
      warning("'time' argument provided, but UD is static.")
    }
  }

  plot_ids <- if (compact) "all" else as.character(unique(track_df$id))
  plot_list <- list()

  for (pid in plot_ids) {
    title_text <- if (compact) "Estimated utilization distribution and tracks" else paste("Estimated UD and tracks - ID:", pid)
    fill_label <- ifelse(log, expression(log(pi)), expression(pi))

    plot_list[[pid]] <- .build_langevin_plot(
      track_df = track_df, pid = pid, raster_obj = ud_full, user_extent = extent,
      compact = compact, title_text = title_text, fill_label = fill_label,
      track_colors = c("Observed" = "lightgrey", "Estimated" = "tomato"),
      track_lines = c("Observed" = "dashed", "Estimated" = "solid")
    )
  }

  return(if (compact || length(plot_list) == 1) plot_list[[1]] else plot_list)
}

#' @rdname plot.langevin
#' @method plot simLangevin
#' @export
plot.simLangevin <- function(x, spatialCovs, beta = NULL, log = TRUE, extent = NULL, compact = TRUE, ...) {
  if (missing(spatialCovs)) stop("You must provide the 'spatialCovs' list.")

  # If no beta is provided, fall back to dataLangevin plotting behavior
  if (is.null(beta)) {
    return(plot.dataLangevin(x, spatialCovs, extent = extent, compact = compact, ...))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package \"ggplot2\" needed for plotting. Please install it.", call. = FALSE)

  if (length(beta) != length(spatialCovs)) {
    stop("The length of 'beta' must match the number of covariates in 'spatialCovs'.")
  }

  # --- Prepare Tracks ---
  obs_df <- data.frame(x = x$x, y = x$y, id = as.character(x$id), type = "Observed")
  true_df <- data.frame(x = x$mu.x, y = x$mu.y, id = as.character(x$id), type = "True (mu)")

  track_df <- rbind(obs_df, true_df)
  track_df$type <- factor(track_df$type, levels = c("Observed", "True (mu)"))

  track_colors <- c("Observed" = "lightgrey", "True (mu)" = "tomato")
  track_lines <- c("Observed" = "dashed", "True (mu)" = "solid")

  # --- Prepare Theoretical UD Raster ---
  ud_full <- getUD(spatialCovs = spatialCovs, beta = beta, log = log)

  plot_ids <- if (compact) "all" else as.character(unique(track_df$id))
  plot_list <- list()

  for (pid in plot_ids) {
    title_text <- if (compact) "Theoretical utilization distribution and tracks" else paste("Theoretical UD and tracks - ID:", pid)
    fill_label <- ifelse(log, expression(log(pi)), expression(pi))

    plot_list[[pid]] <- .build_langevin_plot(
      track_df = track_df, pid = pid, raster_obj = ud_full, user_extent = extent,
      compact = compact, title_text = title_text, fill_label = fill_label,
      track_colors = track_colors, track_lines = track_lines
    )
  }

  return(if (compact || length(plot_list) == 1) plot_list[[1]] else plot_list)
}

#' @rdname plot.langevin
#' @method plot dataLangevin
#' @export
plot.dataLangevin <- function(x, spatialCovs, extent = NULL, compact = TRUE, ...) {

  if (missing(spatialCovs)) stop("You must provide the 'spatialCovs' list.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package \"ggplot2\" needed for plotting. Please install it.", call. = FALSE)

  is_sim <- inherits(x, "simLangevin")

  # --- Prepare Tracks ---
  obs_df <- data.frame(x = x$x, y = x$y, id = as.character(x$id), type = "Observed")

  if (is_sim) {
    true_df <- data.frame(x = x$mu.x, y = x$mu.y, id = as.character(x$id), type = "True (mu)")
    track_df <- rbind(obs_df, true_df)
    track_df$type <- factor(track_df$type, levels = c("Observed", "True (mu)"))
    track_colors <- c("Observed" = "lightgrey", "True (mu)" = "tomato")
    track_lines <- c("Observed" = "dashed", "True (mu)" = "solid")
  } else {
    track_df <- obs_df
    track_df$type <- factor(track_df$type, levels = c("Observed"))
    track_colors <- c("Observed" = "lightgrey")
    track_lines <- c("Observed" = "dashed")
  }

  cov_names <- names(spatialCovs)
  if (is.null(cov_names)) cov_names <- paste0("Covariate_", seq_along(spatialCovs))

  plot_ids <- if (compact) "all" else as.character(unique(track_df$id))
  master_plot_list <- list()

  for (c_idx in seq_along(spatialCovs)) {
    cov_name <- cov_names[c_idx]
    cov_rast <- spatialCovs[[c_idx]] # We pass the FULL raster, no matter how many layers!
    cov_plot_list <- list()

    for (pid in plot_ids) {
      title_text <- if (compact) paste("Covariate:", cov_name, "- All Tracks") else paste("Covariate:", cov_name, "- Track ID:", pid)

      cov_plot_list[[pid]] <- .build_langevin_plot(
        track_df = track_df, pid = pid, raster_obj = cov_rast, user_extent = extent,
        compact = compact, title_text = title_text, fill_label = cov_name,
        track_colors = track_colors, track_lines = track_lines
      )
    }
    master_plot_list[[cov_name]] <- if (compact) cov_plot_list[[1]] else cov_plot_list
  }

  return(master_plot_list)
}

# --- Internal Helper for Plotting Engine ---
.build_langevin_plot <- function(track_df, pid, raster_obj, user_extent, compact, title_text, fill_label, track_colors, track_lines) {

  trk_sub <- if (compact) track_df else track_df[track_df$id == pid, ]

  # 1. Calculate or use extent
  current_extent <- user_extent
  if (is.null(current_extent)) {
    x_range <- max(trk_sub$x, na.rm = TRUE) - min(trk_sub$x, na.rm = TRUE)
    y_range <- max(trk_sub$y, na.rm = TRUE) - min(trk_sub$y, na.rm = TRUE)
    max_range <- max(x_range, y_range, 1) # Fallback to 1 if point is singular
    x_mid <- (max(trk_sub$x, na.rm = TRUE) + min(trk_sub$x, na.rm = TRUE)) / 2
    y_mid <- (max(trk_sub$y, na.rm = TRUE) + min(trk_sub$y, na.rm = TRUE)) / 2
    current_extent <- c(x_mid - 0.6 * max_range, x_mid + 0.6 * max_range, y_mid - 0.6 * max_range, y_mid + 0.6 * max_range)
  }

  xlim <- ylim <- NULL
  if (!is.null(current_extent)) {
    crop_ext <- tryCatch(terra::ext(current_extent), error = function(e) NULL)
    if (!is.null(crop_ext)) {
      ext_vec <- as.vector(crop_ext)
      xlim <- c(ext_vec["xmin"], ext_vec["xmax"])
      ylim <- c(ext_vec["ymin"], ext_vec["ymax"])
      raster_obj <- tryCatch(terra::crop(raster_obj, crop_ext), error = function(e) raster_obj)
    }
  }

  # 2. Rename raster layers to ensure distinct facet names (borrowing from plotRaster logic)
  n_layers <- terra::nlyr(raster_obj)
  layer_times <- terra::time(raster_obj)

  if (n_layers > 1) {
    if (!is.null(layer_times) && !all(is.na(layer_times))) {
      names(raster_obj) <- paste0("Time: ", layer_times)
    } else {
      names(raster_obj) <- paste0("Layer ", seq_len(n_layers))
    }
  }

  # 3. Raster to Dataframe
  rast_df <- terra::as.data.frame(raster_obj, xy = TRUE, na.rm = TRUE)
  layer_names <- setdiff(names(rast_df), c("x", "y"))

  rast_long <- do.call(rbind, lapply(layer_names, function(lyr) {
    data.frame("x" = rast_df$x, "y" = rast_df$y, "time_layer" = lyr, "Value" = rast_df[[lyr]])
  }))

  # Lock factor levels to maintain chronological/layer order in facets
  rast_long$time_layer <- factor(rast_long$time_layer, levels = layer_names)

  # 4. Build ggplot
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = rast_long, ggplot2::aes(x = x, y = y, fill = Value)) +
    ggplot2::scale_fill_viridis_c(name = fill_label, option = "viridis", na.value = "transparent") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "easting (x)", y = "northing (y)", title = title_text)

  # Track handling
  if(length(track_colors) > 1) {
    p <- p +
      ggplot2::geom_path(data = trk_sub, ggplot2::aes(x = x, y = y, color = type, linetype = type, group = interaction(id, type)), linewidth = 0.5, alpha = 0.8) +
      ggplot2::geom_point(data = trk_sub, ggplot2::aes(x = x, y = y, color = type, shape = type), size = 1, alpha = 0.8) +
      ggplot2::scale_color_manual(name = "Data Type", values = track_colors) +
      ggplot2::scale_linetype_manual(name = "Data Type", values = track_lines) +
      ggplot2::scale_shape_manual(name = "Data Type", values = c(16, 16))
  } else {
    p <- p +
      ggplot2::geom_path(data = trk_sub, ggplot2::aes(x = x, y = y, group = id), color = track_colors[1], linetype = track_lines[1], linewidth = 0.5, alpha = 0.8) +
      ggplot2::geom_point(data = trk_sub, ggplot2::aes(x = x, y = y), color = track_colors[1], shape = 16, size = 1, alpha = 0.8)
  }

  p <- p + if (!is.null(xlim) && !is.null(ylim)) ggplot2::coord_equal(xlim = xlim, ylim = ylim) else ggplot2::coord_equal()

  # THIS triggers facet_wrap for multi-layer rasters!
  if (length(layer_names) > 1) p <- p + ggplot2::facet_wrap(~ time_layer)

  return(p)
}
