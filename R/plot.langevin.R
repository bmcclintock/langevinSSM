#' Plot Langevin model fits and data
#'
#' These methods plot objects associated with the \code{langevinSSM} package using \code{ggplot2}.
#' \itemize{
#'   \item \strong{\code{plot.fitLangevin}:} Plots the estimated Utilization Distribution (UD) and overlays the estimated true locations (mu). If the original data is provided, observed locations are also plotted.
#'   \item \strong{\code{plot.dataLangevin}:} Plots the spatial covariates and overlays the observed locations. If the data is a simulated \code{simLangevin} object, both true (latent) and observed locations are plotted.
#'   \item \strong{\code{plot.simLangevin}:} If \code{beta} is provided, plots the theoretical UD based on those coefficients and overlays true (latent) and observed locations. If \code{beta} is omitted, defaults to \code{plot.dataLangevin} behavior (plotting individual covariates).
#'   \item \strong{\code{plotUD}:} Plots the estimated Utilization Distribution (UD). If the SpatRaster stack contains uncertainty metrics (SE and CV), these are also plotted. A \code{log} argument allows plotting the log of the SE.
#' }
#'
#' @param x A \code{fitLangevin}, \code{dataLangevin}, \code{simLangevin}, or \code{udLangevin} object.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects. Used to compute the UD or plotted as the background.
#' @param beta Optional numeric vector of habitat selection coefficients (for \code{simLangevin} only). Must match the length of \code{spatialCovs}. If provided, plots the UD instead of individual covariates.
#' @param log Logical. For \code{fitLangevin} and \code{simLangevin}, indicates whether to plot the log UD (default: \code{TRUE}) or the probability UD. For \code{udLangevin}, indicates whether to plot the log of the standard error (\code{TRUE}) or natural standard error (default: \code{TRUE}).
#' @param extent Optional. A numeric vector of length 4 \code{c(xmin, xmax, ymin, ymax)} or a \code{\link[terra]{SpatExtent}} object defining the bounding box. If \code{NULL} (default), the extent is automatically calculated from the track data.
#' @param data Optional \code{dataLangevin} object (for \code{fitLangevin} only). If provided, observed coordinates will be plotted beneath the estimated locations.
#' @param time Optional. Indicates which layer(s) of a dynamic UD or covariate to plot. Can be a numeric index, a layer name, or a \code{POSIXct}/\code{Date} object. If \code{NULL} (default), all layers are plotted.
#' @param compact Logical indicating whether to plot all tracks on a single panel (\code{TRUE}, default) or plot each track separately (\code{FALSE}).
#' @param ... Additional arguments passed to \code{\link{plotRaster}}.
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

  ud_full <- getUD(spatialCovs = spatialCovs, beta = beta_est, log = log, plot = FALSE)
  ud_layer_name <- if (log) "log_UD" else "UD"

  # extract all layers that match the target name
  target_idx <- which(names(ud_full) == ud_layer_name)
  ud_raster <- ud_full[[target_idx]]

  if (!is.null(time) && terra::nlyr(ud_raster) == 1) {
    warning("'time' argument provided, but the resulting UD is static. Ignoring 'time'.")
  }

  plot_ids <- if (compact) "all" else as.character(unique(track_df$id))
  plot_list <- list()

  for (pid in plot_ids) {
    title_text <- if (compact) "Estimated utilization distribution and tracks" else paste("Estimated UD and tracks - ID:", pid)
    fill_label <- ifelse(log, expression(log(pi)), expression(pi))

    plot_list[[pid]] <- .build_langevin_plot(
      track_df = track_df, pid = pid, raster_obj = ud_raster, user_extent = extent, time = time,
      compact = compact, title_text = title_text, fill_label = fill_label,
      track_colors = c("Observed" = "lightgrey", "Estimated" = "tomato"),
      track_lines = c("Observed" = "dashed", "Estimated" = "solid"), ...
    )
  }

  return(if (compact || length(plot_list) == 1) plot_list[[1]] else plot_list)
}

#' @rdname plot.langevin
#' @method plot simLangevin
#' @export
plot.simLangevin <- function(x, spatialCovs, beta = NULL, log = TRUE, extent = NULL, time = NULL, compact = TRUE, ...) {
  if (missing(spatialCovs)) stop("You must provide the 'spatialCovs' list.")

  # If no beta is provided, fall back to dataLangevin plotting behavior
  if (is.null(beta)) {
    return(plot.dataLangevin(x, spatialCovs, extent = extent, time = time, compact = compact, ...))
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

  ud_full <- getUD(spatialCovs = spatialCovs, beta = beta, log = log, plot = FALSE)
  ud_layer_name <- if (log) "log_UD" else "UD"

  # extract all layers that match the target name
  target_idx <- which(names(ud_full) == ud_layer_name)
  ud_raster <- ud_full[[target_idx]]

  if (!is.null(time) && terra::nlyr(ud_raster) == 1) {
    warning("'time' argument provided, but the resulting UD is static. Ignoring 'time'.")
  }

  plot_ids <- if (compact) "all" else as.character(unique(track_df$id))
  plot_list <- list()

  for (pid in plot_ids) {
    title_text <- if (compact) "Theoretical utilization distribution and tracks" else paste("Theoretical UD and tracks - ID:", pid)
    fill_label <- ifelse(log, expression(log(pi)), expression(pi))

    plot_list[[pid]] <- .build_langevin_plot(
      track_df = track_df, pid = pid, raster_obj = ud_raster, user_extent = extent, time = time,
      compact = compact, title_text = title_text, fill_label = fill_label,
      track_colors = track_colors, track_lines = track_lines, ...
    )
  }

  return(if (compact || length(plot_list) == 1) plot_list[[1]] else plot_list)
}

#' @rdname plot.langevin
#' @method plot dataLangevin
#' @export
plot.dataLangevin <- function(x, spatialCovs, extent = NULL, time = NULL, compact = TRUE, ...) {

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

  if (!is.null(time) && all(vapply(spatialCovs, terra::nlyr, numeric(1)) == 1)) {
    warning("'time' argument provided, but all spatial covariates are static. Ignoring 'time'.")
  }

  cov_names <- names(spatialCovs)
  if (is.null(cov_names)) cov_names <- paste0("Covariate_", seq_along(spatialCovs))

  plot_ids <- if (compact) "all" else as.character(unique(track_df$id))
  main_plot_list <- list()

  for (c_idx in seq_along(spatialCovs)) {
    cov_name <- cov_names[c_idx]
    cov_rast <- spatialCovs[[c_idx]]

    names(cov_rast) <- rep(cov_name, terra::nlyr(cov_rast))

    cov_plot_list <- list()

    for (pid in plot_ids) {
      title_text <- if (compact) paste("Covariate:", cov_name, "- All Tracks") else paste("Covariate:", cov_name, "- Track ID:", pid)

      cov_plot_list[[pid]] <- .build_langevin_plot(
        track_df = track_df, pid = pid, raster_obj = cov_rast, user_extent = extent, time = time,
        compact = compact, title_text = title_text, fill_label = cov_name,
        track_colors = track_colors, track_lines = track_lines, ...
      )
    }
    main_plot_list[[cov_name]] <- if (compact) cov_plot_list[[1]] else cov_plot_list
  }

  return(main_plot_list)
}

#' @details Because \code{getUD} returns a standard \code{\link[terra]{SpatRaster}} object, users are free to bypass \code{plotUD} and visualize the rasters using base \code{plot()}, \code{ggplot2}, or \code{tidyterra} to suit their specific needs.
#' @rdname plot.langevin
#' @export
plotUD <- function(x, log = TRUE, extent = NULL, time = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package \"ggplot2\" needed for plotting. Please install it.", call. = FALSE)

  plot_list <- list()
  layer_names <- names(x)

  # 1. Plot Base UD
  ud_name <- if (any(grepl("^log_UD", layer_names))) "log_UD" else "UD"

  # Safely extract all layers that match the target name
  target_idx <- which(layer_names == ud_name)
  ud_rast <- x[[target_idx]]

  is_log_ud <- (ud_name == "log_UD")

  ud_title <- ifelse(is_log_ud, "Utilization distribution (log scale)", "Utilization distribution")
  ud_legend <- expression(pi(x)) # if (is_log_ud) expression(log(pi(x))) else expression(pi(x))

  plot_list[["UD"]] <- plotRaster(ud_rast, legend.title = ud_legend, extent = extent, time = time, ...) +
    ggplot2::labs(title = ud_title)

  # Define common legend expressions for uncertainty plots
  se_legend <- "SE" # if (log) expression(log(SE(pi(x)))) else expression(SE(pi(x)))
  cv_legend <- "CV" #expression(CV(pi(x)))

  # 2. Plot Delta Method Uncertainty
  if ("UD_SE_delta" %in% layer_names && "UD_CV_delta" %in% layer_names) {
    se_idx <- which(layer_names == "UD_SE_delta")
    cv_idx <- which(layer_names == "UD_CV_delta")

    se_rast_delta <- x[[se_idx]]

    if (log) {
      se_rast_delta <- log(se_rast_delta)
      names(se_rast_delta) <- rep("log_UD_SE_delta", terra::nlyr(se_rast_delta))
    }

    plot_list[["SE_delta"]] <- plotRaster(se_rast_delta, legend.title = se_legend, extent = extent, time = time, ...) +
      ggplot2::labs(title = ifelse(log, "UD log standard error (Delta method)", "UD standard error (Delta method)"))

    plot_list[["CV_delta"]] <- plotRaster(x[[cv_idx]], legend.title = cv_legend, extent = extent, time = time, ...) +
      ggplot2::labs(title = "UD coefficient of variation (Delta method)")
  }

  # 3. Plot Simulated Uncertainty (Welford's Algorithm)
  if ("UD_SE_sim" %in% layer_names && "UD_CV_sim" %in% layer_names) {
    se_sim_idx <- which(layer_names == "UD_SE_sim")
    cv_sim_idx <- which(layer_names == "UD_CV_sim")

    se_rast_sim <- x[[se_sim_idx]]

    if (log) {
      se_rast_sim <- log(se_rast_sim)
      names(se_rast_sim) <- rep("log_UD_SE_sim", terra::nlyr(se_rast_sim))
    }

    plot_list[["SE_sim"]] <- plotRaster(se_rast_sim, legend.title = se_legend, extent = extent, time = time, ...) +
      ggplot2::labs(title = ifelse(log, "UD log standard error (simulated)", "UD standard error (simulated)"))

    plot_list[["CV_sim"]] <- plotRaster(x[[cv_sim_idx]], legend.title = cv_legend, extent = extent, time = time, ...) +
      ggplot2::labs(title = "UD coefficient of variation (simulated)")
  }

  if (length(plot_list) == 1) return(plot_list[[1]])

  return(plot_list)
}

# --- Internal Helper for Plotting Engine ---
.build_langevin_plot <- function(track_df, pid, raster_obj, user_extent, time, compact, title_text, fill_label, track_colors, track_lines, ...) {

  trk_sub <- if (compact) track_df else track_df[track_df$id == pid, ]

  current_extent <- user_extent
  if (is.null(current_extent)) {
    x_range <- max(trk_sub$x, na.rm = TRUE) - min(trk_sub$x, na.rm = TRUE)
    y_range <- max(trk_sub$y, na.rm = TRUE) - min(trk_sub$y, na.rm = TRUE)
    max_range <- max(x_range, y_range, 1) # Fallback to 1 if point is singular
    x_mid <- (max(trk_sub$x, na.rm = TRUE) + min(trk_sub$x, na.rm = TRUE)) / 2
    y_mid <- (max(trk_sub$y, na.rm = TRUE) + min(trk_sub$y, na.rm = TRUE)) / 2
    current_extent <- c(x_mid - 0.6 * max_range, x_mid + 0.6 * max_range, y_mid - 0.6 * max_range, y_mid + 0.6 * max_range)
  }

  crop_ext <- NULL
  if (!is.null(current_extent)) {
    crop_ext <- tryCatch(terra::ext(current_extent), error = function(e) NULL)
    if (!is.null(crop_ext)) {
      raster_obj <- tryCatch({
        terra::crop(raster_obj, crop_ext)
      }, error = function(e) {
        warning("terra::crop failed. Ignoring extent.")
        return(raster_obj)
      })
    }
  }

  p <- plotRaster(raster_obj, legend.title = fill_label, extent = crop_ext, time = time, ...)

  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "easting (x)", y = "northing (y)", title = title_text)

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

  return(p)
}
