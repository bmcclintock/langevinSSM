#' Plot Langevin model fits and data
#'
#' These methods plot objects associated with the \code{langevinSSM} package using \code{ggplot2}.
#' \itemize{
#'   \item \strong{\code{plot.fitLangevin}:} Plots the estimated Utilization Distribution (UD) and overlays the estimated true locations (mu). If the original data is provided, observed locations are also plotted.
#'   \item \strong{\code{plot.dataLangevin}:} Plots the spatial covariates and overlays the observed locations. If the data is a simulated \code{simLangevin} object, both true (latent) and observed locations are plotted.
#'   \item \strong{\code{plot.simLangevin}:} If \code{beta} is provided, plots the theoretical UD based on those coefficients and overlays true (latent) and observed locations. If \code{beta} is omitted, defaults to \code{plot.dataLangevin} behavior (plotting tracks on spatial covariates).
#'   \item \strong{\code{plot.regLangevin}:} Plots the region of interest defined by the mask and the relative probability of presence within that region.
#'   \item \strong{\code{plot.resLangevin}:} Generates Q-Q and ACF diagnostic plots for One-Step-Ahead (OSA) residuals.
#'   \item \strong{\code{plotUD}:} Plots the estimated utilization distribution (UD). If the SpatRaster stack contains uncertainty metrics (SE and CV), these are also plotted. A \code{log} argument allows plotting the log of the SE.
#' }
#'
#' @param x A \code{fitLangevin}, \code{dataLangevin}, \code{simLangevin}, \code{regLangevin}, \code{udLangevin}, or \code{resLangevin} object.
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects. Used to compute the UD or plotted as the background.
#' @param beta Optional numeric vector of habitat selection coefficients (for \code{simLangevin} only). Must match the length of \code{spatialCovs} minus the barrier (if present). If provided, plots the UD instead of individual covariates.
#' @param log Logical. Indicates whether to plot the Utilization Distribution (UD) on the log scale (\code{TRUE}) or the probability scale (\code{FALSE}). For \code{plot.regLangevin}, the default is \code{FALSE}. For all other UD plotting methods, the default is \code{TRUE}. When plotting a \code{SpatRaster} that contains uncertainty metrics via \code{plotUD}, this argument also toggles the standard error layers between log-scale SE and natural-scale SE.
#' @param extent Optional. A numeric vector of length 4 \code{c(xmin, xmax, ymin, ymax)} or a \code{\link[terra]{SpatExtent}} object defining the bounding box. If \code{NULL} (default), the extent is automatically calculated from the track data.
#' @param data Optional \code{dataLangevin} object (for \code{fitLangevin} only). If provided, observed coordinates will be plotted beneath the estimated locations.
#' @param time Optional. Indicates which layer(s) of a dynamic UD or covariate to plot. Can be a numeric index, a layer name, or a \code{POSIXct}/\code{Date} object. If \code{NULL} (default), all layers are plotted.
#' @param compact Logical indicating whether to plot all tracks on a single panel (\code{TRUE}, default) or plot each track separately (\code{FALSE}).
#' @param tracks Optional. Vector of track IDs to plot separately, or \code{"all"} to plot each track individually. If \code{NULL} (default), residuals for all tracks are aggregated into a single set of plots (used only for \code{plot.resLangevin}).
#' @param maskRast \code{\link[terra]{SpatRaster-class}} object for areas to be masked out (set to \code{0}) before plotting the UD. Default: \code{NULL} (no mask).
#' @param ... Additional arguments passed to internal plotting methods.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object, or a list of \code{ggplot} objects depending on the input type and \code{compact} or \code{tracks} argument.
#'
#' @name plot.langevin
#' @importFrom terra as.data.frame nlyr time crop ext ifel mask trim compareGeom
NULL

#' @rdname plot.langevin
#' @method plot fitLangevin
#' @export
plot.fitLangevin <- function(x, spatialCovs, log = TRUE, extent = NULL, data = NULL, time = NULL, compact = TRUE, maskRast = NULL, ...) {
  if (missing(spatialCovs)) stop("'spatialCovs' must be provided.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package \"ggplot2\" needed for plotting. Please install it.", call. = FALSE)

  verify_signatures(x, data = data, spatialCovs = spatialCovs)
  boundsWarning(x)

  # --- Prepare Tracks ---
  has_est <- !is.null(x$estimates$random$mu$est)
  track_df <- data.frame()

  if (has_est) {
    track_id <- x$estimates$random$mu$est$id
    if (is.null(track_id)) track_id <- rep("1", nrow(x$estimates$random$mu$est))

    track_df <- data.frame(
      x = x$estimates$random$mu$est[, "mu.x"],
      y = x$estimates$random$mu$est[, "mu.y"],
      id = as.character(track_id),
      type = "Estimated"
    )
  }

  if (!is.null(data)) {
    if (inherits(data, "dataLangevin")) {
      obs_df <- data.frame(x = data$x, y = data$y, id = as.character(data$id), type = "Observed")
      track_df <- rbind(obs_df, track_df)
    } else {
      warning("'data' is not a dataLangevin object. Skipping observed locations.")
    }
  }

  # style the tracks based on what is available
  if (nrow(track_df) == 0) {
    if (!has_est) message("Note: Model was fit without measurement error. Provide 'data' to overlay the tracks.")
  } else {
    if (has_est && !is.null(data)) {
      track_df$type <- factor(track_df$type, levels = c("Observed", "Estimated"))
      track_colors <- c("Observed" = "lightgrey", "Estimated" = "#E69F00")
      track_lines <- c("Observed" = "dashed", "Estimated" = "solid")
    } else if (has_est && is.null(data)) {
      track_df$type <- factor(track_df$type, levels = c("Estimated"))
      track_colors <- c("Estimated" = "#E69F00")
      track_lines <- c("Estimated" = "solid")
    } else if (!has_est && !is.null(data)) {
      track_df$type <- factor(track_df$type, levels = c("Observed"))
      track_colors <- c("Observed" = "black")
      track_lines <- c("Observed" = "solid")
    }
  }

  barrier <- x$conditions$barrier
  lambda <- x$conditions$lambda
  scaleFactor <- x$conditions$scaleFactor
  rn <- rownames(x$estimates$natural)
  beta_est <- x$estimates$natural[which(grepl("^beta", rn)), "Estimate"]

  ud_full <- getUD(spatialCovs = spatialCovs, beta = beta_est, barrier=barrier, lambda = lambda, log = log, plot = FALSE, maskRast = maskRast, scaleFactor = scaleFactor)
  ud_layer_name <- if (log) "log_UD" else "UD"

  # extract all layers that match the target name
  target_idx <- which(names(ud_full) == ud_layer_name)
  ud_raster <- ud_full[[target_idx]]

  if (!is.null(time) && terra::nlyr(ud_raster) == 1) {
    warning("'time' argument provided, but the resulting UD is static. Ignoring 'time'.")
  }

  plot_ids <- if (compact || nrow(track_df) == 0) "all" else as.character(unique(track_df$id))
  plot_list <- list()

  for (pid in plot_ids) {
    title_text <- if (compact || nrow(track_df) == 0) "Estimated utilization distribution and tracks" else paste("Estimated UD and tracks - ID:", pid)
    fill_label <- if (log) expression(log(pi(x))) else expression(pi(x))

    plot_list[[pid]] <- .build_langevin_plot(
      track_df = track_df, pid = pid, raster_obj = ud_raster, user_extent = extent, time = time,
      compact = compact, title_text = title_text, fill_label = fill_label,
      track_colors = track_colors,
      track_lines = track_lines,
      maskRast = maskRast, ...
    )
  }

  return(if (compact || length(plot_list) == 1) plot_list[[1]] else plot_list)
}

#' @rdname plot.langevin
#' @method plot simLangevin
#' @export
plot.simLangevin <- function(x, spatialCovs, beta = NULL, log = TRUE, extent = NULL, time = NULL, compact = TRUE, maskRast = NULL, ...) {
  if (missing(spatialCovs)) stop("You must provide the 'spatialCovs' list.")

  if (is.null(beta)){
    return(plot.dataLangevin(x, spatialCovs, extent = extent, time = time, compact = compact, maskRast = maskRast, ...))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package \"ggplot2\" needed for plotting. Please install it.", call. = FALSE)

  barrier <- attr(x, "barrier")
  lambda <- attr(x, "lambda")
  sf_attr <- attr(x, "scaleFactor")
  scaleFactor <- if (!is.null(sf_attr)) sf_attr else 1

  # Check that lengths line up after handling the barrier
  expected_length <- length(spatialCovs)
  if (!is.null(barrier) && barrier %in% names(spatialCovs)) {
    expected_length <- expected_length - 1
  }

  if (length(beta) != expected_length) {
    stop(sprintf("The length of 'beta' must match the number of habitat covariates in 'spatialCovs' (%d).", expected_length))
  }

  # --- Prepare Tracks ---
  err_cols <- c("smaj", "smin", "eor", "x.err", "y.err")
  existing_err_cols <- intersect(err_cols, names(x))
  has_error <- FALSE
  if (length(existing_err_cols) > 0) {
    has_error <- any(!is.na(x[existing_err_cols]))
  }

  true_df <- data.frame(x = x$mu.x, y = x$mu.y, id = as.character(x$id), type = "True")

  if (has_error) {
    obs_df <- data.frame(x = x$x, y = x$y, id = as.character(x$id), type = "Observed")
    track_df <- rbind(obs_df, true_df)
    track_df$type <- factor(track_df$type, levels = c("Observed", "True"))

    track_colors <- c("Observed" = "lightgrey", "True" = "#E69F00")
    track_lines <- c("Observed" = "dashed", "True" = "solid")
  } else {
    track_df <- true_df
    track_df$type <- factor(track_df$type, levels = c("True"))

    track_colors <- c("True" = "#E69F00")
    track_lines <- c("True" = "solid")
  }

  ud_full <- getUD(spatialCovs = spatialCovs, beta = beta, barrier = barrier, lambda = lambda, log = log, plot = FALSE, maskRast = maskRast, scaleFactor = scaleFactor)
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
    fill_label <- if (log) expression(log(pi(x))) else expression(pi(x))

    plot_list[[pid]] <- .build_langevin_plot(
      track_df = track_df, pid = pid, raster_obj = ud_raster, user_extent = extent, time = time,
      compact = compact, title_text = title_text, fill_label = fill_label,
      track_colors = track_colors, track_lines = track_lines,
      maskRast = maskRast, ...
    )
  }

  return(if (compact || length(plot_list) == 1) plot_list[[1]] else plot_list)
}

#' @rdname plot.langevin
#' @method plot dataLangevin
#' @export
plot.dataLangevin <- function(x, spatialCovs, extent = NULL, time = NULL, compact = TRUE, maskRast = NULL, ...) {

  if (missing(spatialCovs)) stop("You must provide the 'spatialCovs' list.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package \"ggplot2\" needed for plotting. Please install it.", call. = FALSE)

  is_sim <- inherits(x, "simLangevin")

  # --- Prepare Tracks ---
  if (is_sim) {
    err_cols <- c("smaj", "smin", "eor", "x.err", "y.err")
    existing_err_cols <- intersect(err_cols, names(x))
    has_error <- FALSE
    if (length(existing_err_cols) > 0) {
      has_error <- any(!is.na(x[existing_err_cols]))
    }

    true_df <- data.frame(x = x$mu.x, y = x$mu.y, id = as.character(x$id), type = "True")

    if (has_error) {
      obs_df <- data.frame(x = x$x, y = x$y, id = as.character(x$id), type = "Observed")
      track_df <- rbind(obs_df, true_df)
      track_df$type <- factor(track_df$type, levels = c("Observed", "True"))
      track_colors <- c("Observed" = "lightgrey", "True" = "#E69F00")
      track_lines <- c("Observed" = "dashed", "True" = "solid")
    } else {
      track_df <- true_df
      track_df$type <- factor(track_df$type, levels = c("True"))
      track_colors <- c("True" = "#E69F00")
      track_lines <- c("True" = "solid")
    }
  } else {
    obs_df <- data.frame(x = x$x, y = x$y, id = as.character(x$id), type = "Observed")
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
        track_colors = track_colors, track_lines = track_lines,
        maskRast = maskRast, ...
      )
    }
    main_plot_list[[cov_name]] <- if (compact) cov_plot_list[[1]] else cov_plot_list
  }

  return(main_plot_list)
}

#' @rdname plot.langevin
#' @method plot resLangevin
#' @importFrom stats acf qnorm qchisq ppoints quantile dnorm dchisq na.omit
#' @export
plot.resLangevin <- function(x, tracks = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {

    stop("Package 'ggplot2' is required for plotting residuals. Please install it.")
  }

  resids <- x

  if (is.null(tracks)) {
    return(.generate_res_plots(resids))
  } else {
    if (length(tracks) == 1 && tracks[1] == "all") {
      tracks_to_plot <- unique(resids$id)
    } else {
      tracks_to_plot <- tracks
    }

    out_list <- list()
    for (trk in tracks_to_plot) {
      sub_data <- resids[resids$id == trk, ]
      if (nrow(sub_data) > 0) {
        out_list[[as.character(trk)]] <- .generate_res_plots(sub_data, paste0(" - ID: ", trk))
      } else {
        warning("No data found for track ID: ", trk, ". Skipping this track.")
      }
    }
    return(out_list)
  }
}

#' @rdname plot.langevin
#' @method plot regLangevin
#' @export
plot.regLangevin <- function(x, extent = NULL, log = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package \"ggplot2\" needed for plotting. Please install it.", call. = FALSE)

  region_prob <- x$prob_raster

  if (log) {
    region_prob <- log(region_prob)
    fill_label <- expression(log(pi(x)))
  } else {
    fill_label <- expression(pi(x))
  }

  # Filter out everything outside the region of interest
  m_na <- terra::ifel(x$mask == 1, 1, NA)

  # terra::mask automatically recycles a 1-layer mask across a multi-layer raster
  region_prob <- terra::mask(region_prob, m_na)

  # Find the true bounding box of the active region
  active_ext <- tryCatch(terra::ext(terra::trim(region_prob)), error = function(e) NULL)

  crop_ext <- extent
  if (is.null(crop_ext)) {
    crop_ext <- active_ext
  } else if (!is.null(active_ext)) {
    # Check if the user's extent crops out any part of the active region
    user_ext <- tryCatch(terra::ext(crop_ext), error = function(e) NULL)
    if (!is.null(user_ext)) {
      # Vectors from terra::ext() format as: xmin, xmax, ymin, ymax
      if (active_ext[1] < user_ext[1] || active_ext[2] > user_ext[2] ||
          active_ext[3] < user_ext[3] || active_ext[4] > user_ext[4]) {
        warning("The provided 'extent' crops out parts of the region of interest. The visible pixels will not sum to the total regional probability.")
      }
    }
  }

  n_layers <- terra::nlyr(region_prob)
  prob_strings <- sprintf("%.2f%%", x$Point_Estimate * 100)

  if (n_layers == 1) {
    title_str <- paste("Regional Probability:", prob_strings[1])
  } else {
    title_str <- "Regional Probability"
  }

  p <- .build_langevin_plot(
    track_df = NULL, pid = "all", raster_obj = region_prob, user_extent = crop_ext, time = NULL,
    compact = TRUE, title_text = title_str, fill_label = fill_label,
    ...
  )

  # Intercept and safely overwrite the facet labels created by plotRaster
  if (n_layers > 1) {
    layer_times <- tryCatch(terra::time(region_prob), error = function(e) NULL)

    if (!is.null(layer_times) && !all(is.na(layer_times))) {
      old_names <- paste0("Time: ", layer_times)
    } else {
      old_names <- paste0("Layer ", seq_len(n_layers))
    }

    # Map the old names to the new strings containing the probabilities
    new_names <- paste0(old_names, "\n(Prob: ", prob_strings, ")")
    names(new_names) <- old_names

    p <- p + ggplot2::facet_wrap(~ layer, labeller = ggplot2::as_labeller(new_names))
  }

  return(p)
}

#' @details Because \code{getUD} returns a standard \code{\link[terra]{SpatRaster}} object, users are free to bypass \code{plotUD} and visualize the rasters using base \code{plot()}, \code{ggplot2}, or \code{tidyterra} to suit their specific needs.
#' @rdname plot.langevin
#' @export
plotUD <- function(x, log = TRUE, extent = NULL, time = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package \"ggplot2\" needed for plotting. Please install it.", call. = FALSE)

  if(!inherits(x, "SpatRaster")) stop("'x' must be a SpatRaster object containing the UD and optionally its uncertainty layers.")

  plot_list <- list()
  layer_names <- names(x)

  has_log_ud <- any(grepl("^log_UD$", layer_names))
  has_ud <- any(grepl("^UD$", layer_names))

  ud_name <- if (has_log_ud) "log_UD" else "UD"
  target_idx <- which(layer_names == ud_name)
  if(length(target_idx)) {
    ud_rast <- x[[target_idx]]
  } else {
    stop("No UD layer found in the provided SpatRaster. Expected a layer named 'UD' or 'log_UD'.")
  }

  # Actively transform base UD based on log argument regardless of raster origin
  if (log) {
    if (!has_log_ud) {
      ud_rast <- log(ud_rast)
    }
    ud_title <- "Utilization distribution (log scale)"
    ud_legend <- expression(log(pi(x)))
    se_legend <- "log(SE)"
  } else {
    if (has_log_ud) {
      # Safely normalize the linear predictor back to a probability distribution
      max_log <- terra::global(ud_rast, "max", na.rm = TRUE)$max
      for(k in seq_len(terra::nlyr(ud_rast))) {
        ud_rast[[k]] <- exp(ud_rast[[k]] - max_log[k])
      }
      layer_sums <- terra::global(ud_rast, "sum", na.rm = TRUE)$sum
      for(k in seq_len(terra::nlyr(ud_rast))) {
        ud_rast[[k]] <- ud_rast[[k]] / layer_sums[k]
      }
    }
    ud_title <- "Utilization distribution"
    ud_legend <- expression(pi(x))
    se_legend <- "SE"
  }

  plot_list[["UD"]] <- plotRaster(ud_rast, legend.title = ud_legend, extent = extent, time = time, ...) +
    ggplot2::labs(title = ud_title, fill = ud_legend)

  cv_legend <- "CV"

  # plot Delta method uncertainty
  if ("UD_SE_delta" %in% layer_names && "UD_CV_delta" %in% layer_names) {
    se_idx <- which(layer_names == "UD_SE_delta")
    cv_idx <- which(layer_names == "UD_CV_delta")

    se_rast_delta <- x[[se_idx]]

    if (log) {
      se_rast_delta <- log(se_rast_delta)
      names(se_rast_delta) <- rep("log_UD_SE_delta", terra::nlyr(se_rast_delta))
    }

    plot_list[["SE_delta"]] <- plotRaster(se_rast_delta, legend.title = se_legend, extent = extent, time = time, ...) +
      ggplot2::labs(title = ifelse(log, "UD log standard error (Delta method)", "UD standard error (Delta method)"), fill = se_legend)

    plot_list[["CV_delta"]] <- plotRaster(x[[cv_idx]], legend.title = cv_legend, extent = extent, time = time, ...) +
      ggplot2::labs(title = "UD coefficient of variation (Delta method)", fill = cv_legend)
  }

  # plot Monte Carlo uncertainty
  if ("UD_SE_sim" %in% layer_names && "UD_CV_sim" %in% layer_names) {
    se_sim_idx <- which(layer_names == "UD_SE_sim")
    cv_sim_idx <- which(layer_names == "UD_CV_sim")

    se_rast_sim <- x[[se_sim_idx]]
    if (log) {
      se_rast_sim <- log(se_rast_sim)
      names(se_rast_sim) <- rep("log_UD_SE_sim", terra::nlyr(se_rast_sim))
    }
    plot_list[["SE_sim"]] <- plotRaster(se_rast_sim, legend.title = se_legend, extent = extent, time = time, ...) +
      ggplot2::labs(title = ifelse(log, "UD log standard error (simulated)", "UD standard error (simulated)"), fill = se_legend)

    plot_list[["CV_sim"]] <- plotRaster(x[[cv_sim_idx]], legend.title = cv_legend, extent = extent, time = time, ...) +
      ggplot2::labs(title = "UD coefficient of variation (simulated)", fill = cv_legend)
  }

  # return single plot if no uncertainty, otherwise return the full list
  if (length(plot_list) == 1) return(plot_list[[1]])

  return(plot_list)
}

# --- Internal Helpers for plot.resLangevin ---

.make_acf_plot <- function(x, title) {
  x <- stats::na.omit(x)
  if(length(x) < 2) return(NULL)

  acf_obj <- stats::acf(x, plot = FALSE)
  acf_df <- data.frame("lag" = acf_obj$lag[, 1, 1], "acf" = acf_obj$acf[, 1, 1])
  clim <- stats::qnorm(0.975) / sqrt(acf_obj$n.used)

  ggplot2::ggplot(acf_df, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(yintercept = c(0, -clim, clim),
                        color = c("black", "blue", "blue"),
                        linetype = c("solid", "dashed", "dashed")) +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0)) +
    ggplot2::labs(title = title, x = "Lag", y = "ACF") +
    ggplot2::theme_minimal()
}

.make_qq_norm <- function(x, title) {
  x <- sort(stats::na.omit(x))
  n <- length(x)
  if(n < 2) return(NULL)

  p <- stats::ppoints(n)
  z <- stats::qnorm(p)

  Qx <- stats::quantile(x, c(0.25, 0.75), names = FALSE)
  Qz <- stats::qnorm(c(0.25, 0.75))
  slope <- diff(Qx) / diff(Qz)
  int <- Qx[1] - slope * Qz[1]

  se <- slope / stats::dnorm(z) * sqrt(p * (1 - p) / n)
  fit_line <- int + slope * z
  upper <- fit_line + 1.96 * se
  lower <- fit_line - 1.96 * se

  outlier <- x < lower | x > upper

  df <- data.frame("sample" = x, "theoretical" = z, "lower" = lower, "upper" = upper, "outlier" = outlier)

  ggplot2::ggplot(df, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
    ggplot2::geom_abline(intercept = int, slope = slope, color = "red", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(color = outlier), alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("FALSE" = "royalblue", "TRUE" = "red")) +
    ggplot2::labs(title = title, x = "Theoretical Quantiles", y = "Sample Quantiles") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

.make_qq_chisq <- function(x, title) {
  x <- sort(stats::na.omit(x))
  n <- length(x)
  if(n < 2) return(NULL)

  p <- stats::ppoints(n)
  z <- stats::qchisq(p, df = 2)

  slope <- 1
  int <- 0

  se <- slope / stats::dchisq(z, df = 2) * sqrt(p * (1 - p) / n)
  fit_line <- int + slope * z
  upper <- fit_line + 1.96 * se
  lower <- fit_line - 1.96 * se

  outlier <- x < lower | x > upper

  df <- data.frame("sample" = x, "theoretical" = z, "lower" = lower, "upper" = upper, "outlier" = outlier)

  ggplot2::ggplot(df, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
    ggplot2::geom_abline(intercept = int, slope = slope, color = "red", linewidth = 1, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(color = outlier), alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("FALSE" = "royalblue", "TRUE" = "red")) +
    ggplot2::labs(title = title, x = "Theoretical Quantiles (Chi-sq, df=2)", y = "Sample Squared Mahalan Distance") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

.generate_res_plots <- function(data_subset, title_suffix = "") {
  res_x_full <- data_subset$residual.x[!is.na(data_subset$residual.x)]
  res_y_full <- data_subset$residual.y[!is.na(data_subset$residual.y)]

  valid_idx <- which(!is.na(data_subset$residual.x) & !is.na(data_subset$residual.y))
  res_x_pair <- data_subset$residual.x[valid_idx]
  res_y_pair <- data_subset$residual.y[valid_idx]
  D2 <- res_x_pair^2 + res_y_pair^2

  p_qq_x <- .make_qq_norm(res_x_full, paste0("Normal Q-Q Plot: x residuals", title_suffix))
  p_qq_y <- .make_qq_norm(res_y_full, paste0("Normal Q-Q Plot: y residuals", title_suffix))
  p_acf_x <- .make_acf_plot(res_x_full, paste0("ACF: x residuals", title_suffix))
  p_acf_y <- .make_acf_plot(res_y_full, paste0("ACF: y residuals", title_suffix))

  p_qq_mah <- .make_qq_chisq(D2, paste0("Mahalanobis Q-Q Plot", title_suffix))
  p_acf_mah <- .make_acf_plot(D2, paste0("ACF: Mahalanobis Distance", title_suffix))

  return(list(
    qq_x = p_qq_x, qq_y = p_qq_y,
    acf_x = p_acf_x, acf_y = p_acf_y,
    qq_mah = p_qq_mah, acf_mah = p_acf_mah
  ))
}

# --- Internal Helper for Plotting Engine ---
.build_langevin_plot <- function(track_df, pid, raster_obj, user_extent = NULL, time = NULL, compact = TRUE, title_text = "", fill_label = "", track_colors = NULL, track_lines = NULL, maskRast = NULL, ...) {

  trk_sub <- NULL
  if (!is.null(track_df)) {
    trk_sub <- if (compact) track_df else track_df[track_df$id == pid, ]
  }

  current_extent <- user_extent
  if (is.null(current_extent) && !is.null(trk_sub) && nrow(trk_sub) > 0) {
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
      raster_obj <- tryCatch(
        terra::crop(raster_obj, crop_ext),
        error = function(e) {
          warning("terra::crop failed. Ignoring extent.")
          return(raster_obj)
        }
      )
    }
  }

  if (!is.null(maskRast)) {
    if(!inherits(maskRast, "SpatRaster")) stop("'maskRast' must be a SpatRaster")

    if (!is.null(crop_ext)) {
      maskRast <- tryCatch(
        terra::crop(maskRast, crop_ext),
        error = function(e) {
          warning("terra::crop failed on maskRast. Ignoring maskRast.")
          return(NULL)
        }
      )
    }

    if(!is.null(maskRast)) {
      if (!terra::compareGeom(raster_obj, maskRast, stopOnError = FALSE)) {
        stop("The 'maskRast' raster must share the same projection (CRS), extent, and resolution as the rasters in 'spatialCovs'.")
      }
      raster_obj <- terra::ifel(maskRast <= 0, NA, raster_obj)
    }
  }

  p <- plotRaster(raster_obj, legend.title = fill_label, extent = crop_ext, time = time, ...)

  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "easting (x)", y = "northing (y)", title = title_text, fill = fill_label)

  if (!is.null(trk_sub) && nrow(trk_sub) > 0) {
    if(length(track_colors) > 1) {
      p <- p +
        ggplot2::geom_path(data = trk_sub, ggplot2::aes(x = x, y = y, color = type, linetype = type, group = interaction(id, type)), linewidth = 0.5, alpha = 0.8, na.rm = TRUE) +
        ggplot2::geom_point(data = trk_sub, ggplot2::aes(x = x, y = y, color = type, shape = type), size = 1, alpha = 0.8, na.rm = TRUE) +
        ggplot2::scale_color_manual(name = "Tracks", values = track_colors) +
        ggplot2::scale_linetype_manual(name = "Tracks", values = track_lines) +
        ggplot2::scale_shape_manual(name = "Tracks", values = c(16, 16)) +
        # FORCE LEGEND ORDER: Pin the discrete Tracks legend to the top (order = 1)
        ggplot2::guides(
          color = ggplot2::guide_legend(order = 1),
          linetype = ggplot2::guide_legend(order = 1),
          shape = ggplot2::guide_legend(order = 1)
        )
    } else {
      p <- p +
        ggplot2::geom_path(data = trk_sub, ggplot2::aes(x = x, y = y, group = id), color = track_colors[1], linetype = track_lines[1], linewidth = 0.5, alpha = 0.8, na.rm = TRUE) +
        ggplot2::geom_point(data = trk_sub, ggplot2::aes(x = x, y = y), color = track_colors[1], shape = 16, size = 1, alpha = 0.8, na.rm = TRUE)
    }
  }

  return(p)
}
