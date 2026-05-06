#' Assess mask leakage for a fitted Langevin model or tracking data
#'
#' Evaluates the proportion and severity of track locations that cross into restricted spatial boundaries,
#' by integrating the exact 2D probability density of the observation error or estimated uncertainty
#' over the restricted raster cells.
#'
#' @param x A \code{fitLangevin}, \code{simLangevin}, or \code{dataLangevin} object.
#' @param maskRast A SpatRaster binary mask (1 = allowed, 0 = restricted).
#' @param level Numeric between 0 and 1. The proportion of the observation's spatial density that must fall within the restricted area to be considered a leak. Default is 0.999 (meaning 99.9\% of the density must be restricted).
#' @param tolerance Numeric. A buffer distance (in the spatial units of the raster). Density falling within this distance of allowed cells is not flagged as a leakage. Default is 0.
#' @param coord Character vector of length 2 specifying the names of the coordinate columns in \code{x} (e.g., c("mu.x", "mu.y"), c("x", "y")). Default is c("mu.x", "mu.y").
#' @return A list containing two data frames:
#' \describe{
#'   \item{summary}{A single-row data frame containing overall leakage metrics:
#'     \itemize{
#'       \item \code{Total_Locations}: The total number of valid locations evaluated.
#'       \item \code{Leaked_Locations}: The number of locations meeting the threshold.
#'       \item \code{Percent_Leaked}: The percentage of evaluated locations classified as leaked.
#'       \item \code{Max_Restricted_Density} & \code{Mean_Restricted_Density}: Restricted area density metrics.
#'       \item \code{Max_Leak_Depth} & \code{Mean_Leak_Depth}: Penetration distance metrics.
#'     }
#'   }
#'   \item{leaked_data}{A data frame containing only the offending locations, with columns:
#'     \itemize{
#'       \item \code{obs_index}: The row index of the location.
#'       \item \code{restricted_density}: The proportion of the error distribution in the restricted area.
#'       \item \code{leak_depth}: The distance the point penetrated into the restricted area.
#'     }
#'   }
#' }
#'
#' @importFrom terra as.matrix xFromCol yFromRow res extract ifel
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
maskLeakage <- function(x, maskRast, level = 0.999, tolerance = 0, coord = c("mu.x", "mu.y")) {

  if (!inherits(x, c("fitLangevin", "simLangevin", "dataLangevin"))) {
    stop("'x' must be a fitLangevin, simLangevin, or dataLangevin object.")
  }

  if (missing(maskRast)) stop("'maskRast' must be provided.")
  if(!inherits(maskRast, "SpatRaster")) stop("'maskRast' must be a SpatRaster")

  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("'level' must be a numeric value between 0 and 1 (e.g., 0.999).")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance < 0) {
    stop("'tolerance' must be a single non-negative numeric value.")
  }

  # Extract the coordinates dynamically based on the object class
  if (inherits(x, "fitLangevin")) {
    if(all(coord==c("x", "y"))) stop("'coord' must be the estimated true locations ('mu.x', 'mu.y') for 'fitLangevin' objects.")
    pts_df <- x$estimates$random$mu$est
    se_df <- x$estimates$random$mu$se
  } else {
    pts_df <- as.data.frame(x)
  }
  if(!all(coord %in% names(pts_df))) stop(sprintf("'coord' columns '%s' not found.", paste(coord, collapse = ", ")))

  pts_x <- pts_df[[coord[1]]]
  pts_y <- pts_df[[coord[2]]]
  n_pts <- nrow(pts_df)
  p_restricted <- numeric(n_pts)

  # Calculate leakage depth using SDF
  sdf_rast <- suppressMessages(prepBarrier(maskRast))
  pts <- cbind(pts_x, pts_y)

  # Dynamically extract the last column to avoid ID argument errors
  dist_ext <- terra::extract(sdf_rast, pts)
  dist_vals <- dist_ext[, ncol(dist_ext)]
  leak_depth <- ifelse(!is.na(dist_vals) & dist_vals < 0, abs(dist_vals), 0)

  # Adjust mask based on tolerance buffer
  if (tolerance > 0) {
    maskRast_eval <- terra::ifel(sdf_rast >= -tolerance, 1, 0)
  } else {
    maskRast_eval <- maskRast
  }

  mask_mat <- terra::as.matrix(maskRast_eval, wide = TRUE)
  x_coords <- terra::xFromCol(maskRast_eval, 1:ncol(maskRast_eval))
  y_coords <- terra::yFromRow(maskRast_eval, 1:nrow(maskRast_eval))

  cell_res <- terra::res(maskRast_eval)
  res_x <- cell_res[1]
  res_y <- cell_res[2]
  cell_area <- res_x * res_y

  box_sigma <- 4.0 # 4-sigma captures >99.9% of the 2D error distribution

  fast_bvn <- function(x, y, mux, muy, sigx, sigy, covxy) {
    rho <- covxy / (sigx * sigy)
    rho <- pmax(pmin(rho, 0.999999), -0.999999)

    zx <- (x - mux) / sigx
    zy <- (y - muy) / sigy
    z <- zx^2 - 2 * rho * zx * zy + zy^2

    denom <- 2 * pi * sigx * sigy * sqrt(1 - rho^2)
    return(exp(-z / (2 * (1 - rho^2))) / denom)
  }

  pb <- utils::txtProgressBar(min = 0, max = n_pts, style = 3, width = 50)

  for (i in 1:n_pts) {
    if (is.na(pts_x[i]) || is.na(pts_y[i])) {
      p_restricted[i] <- NA
      utils::setTxtProgressBar(pb, i)
      next
    }

    var_x <- NA_real_
    var_y <- NA_real_
    cov_xy <- 0.0

    if (inherits(x, "fitLangevin")) {
      if (!is.null(se_df)) {
        var_x <- se_df[[coord[1]]][i]^2
        var_y <- se_df[[coord[2]]][i]^2
      }
    } else {
      if ("x.err" %in% names(pts_df) && !is.na(pts_df$x.err[i]) && !is.na(pts_df$y.err[i])) {
        var_x <- pts_df$x.err[i]^2
        var_y <- pts_df$y.err[i]^2
      } else if ("smaj" %in% names(pts_df) && !is.na(pts_df$smaj[i])) {
        M2 <- (pts_df$smaj[i] / sqrt(2.0))^2
        m2 <- (pts_df$smin[i] / sqrt(2.0))^2
        c_rad <- pts_df$eor[i]
        s2c <- sin(c_rad)^2
        c2c <- cos(c_rad)^2

        var_x <- M2 * s2c + m2 * c2c
        var_y <- M2 * c2c + m2 * s2c
        cov_xy <- (M2 - m2) * cos(c_rad) * sin(c_rad)
      }
    }

    # handle zero-error (or missing error) locations as discrete points
    if (is.na(var_x) || is.na(var_y) || var_x == 0 || var_y == 0) {
      pt_ext <- terra::extract(maskRast_eval, matrix(c(pts_x[i], pts_y[i]), ncol = 2))
      dist_val <- pt_ext[, ncol(pt_ext)]
      p_restricted[i] <- ifelse(!is.na(dist_val) && dist_val == 0, 1.0, 0.0)
      utils::setTxtProgressBar(pb, i)
      next
    }

    sig_x <- sqrt(var_x)
    sig_y <- sqrt(var_y)

    # Pad bounding box by half a raster cell to guarantee centroid capture for high-precision GPS
    x_min <- pts_x[i] - box_sigma * sig_x - res_x / 2
    x_max <- pts_x[i] + box_sigma * sig_x + res_x / 2
    y_min <- pts_y[i] - box_sigma * sig_y - res_y / 2
    y_max <- pts_y[i] + box_sigma * sig_y + res_y / 2

    col_idx <- which(x_coords >= x_min & x_coords <= x_max)
    row_idx <- which(y_coords >= y_min & y_coords <= y_max)

    if (length(col_idx) == 0 || length(row_idx) == 0) {
      p_restricted[i] <- NA_real_ # Point is completely off the spatial grid
    } else {
      local_mask <- mask_mat[row_idx, col_idx, drop = FALSE]

      if (all(local_mask == 1, na.rm = TRUE)) {
        p_restricted[i] <- 0.0
      } else if (all(local_mask == 0, na.rm = TRUE)) {
        p_restricted[i] <- 1.0
      } else {
        # Find cells that represent restricted terrain (0)
        rest_rc <- which(local_mask == 0, arr.ind = TRUE)

        if (nrow(rest_rc) > 0) {
          # Map local indices back to geographic coordinates
          rest_x <- x_coords[col_idx[rest_rc[, 2]]]
          rest_y <- y_coords[row_idx[rest_rc[, 1]]]

          dens <- fast_bvn(rest_x, rest_y, pts_x[i], pts_y[i], sig_x, sig_y, cov_xy)
          raw_sum <- sum(dens) * cell_area

          # Clamp to avoid floating point math exceeding 1.0
          p_restricted[i] <- min(raw_sum, 1.0)
        } else {
          p_restricted[i] <- 0.0
        }
      }
    }
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  # ==============================================================================
  # SUMMARIZE RESULTS
  # ==============================================================================
  is_leaked <- !is.na(p_restricted) & (p_restricted >= level)

  total_pts <- sum(!is.na(p_restricted))
  n_leaked <- sum(is_leaked, na.rm = TRUE)
  pct_leaked <- if (total_pts > 0) (n_leaked / total_pts) * 100 else 0

  max_leak <- if (n_leaked > 0) max(p_restricted[is_leaked], na.rm = TRUE) else 0
  mean_leak <- if (n_leaked > 0) mean(p_restricted[is_leaked], na.rm = TRUE) else 0

  max_depth <- if (n_leaked > 0) max(leak_depth[is_leaked], na.rm = TRUE) else 0
  mean_depth <- if (n_leaked > 0) mean(leak_depth[is_leaked], na.rm = TRUE) else 0

  summary_df <- data.frame(
    Total_Locations = total_pts,
    Leaked_Locations = n_leaked,
    Percent_Leaked = round(pct_leaked, 3),
    Max_Restricted_Density = max_leak,
    Mean_Restricted_Density = mean_leak,
    Max_Leak_Depth = max_depth,
    Mean_Leak_Depth = mean_depth
  )

  # Isolate and sort the offending locations
  leaked_df <- pts_df[is_leaked, , drop = FALSE]
  if (n_leaked > 0) {
    leaked_df$obs_index <- which(is_leaked)
    leaked_df$restricted_density <- p_restricted[is_leaked]
    leaked_df$leak_depth <- leak_depth[is_leaked]

    # Reorder columns to place obs_index first
    leaked_df <- leaked_df[, c("obs_index", setdiff(names(leaked_df), "obs_index"))]

    leaked_df <- leaked_df[order(-leaked_df$restricted_density), , drop = FALSE]
    rownames(leaked_df) <- NULL
  }

  message("\n   --- Barrier Leakage Assessment ---")
  if (n_leaked == 0) {
    message(sprintf("      0 locations leaked beyond the %.1f%% density threshold (tolerance: %g).", level * 100, tolerance))
  } else {
    message(sprintf("      %d locations (%.2f%%) leaked beyond the %.1f%% density threshold (tolerance: %g).",
                    n_leaked, pct_leaked, level * 100, tolerance))
    message(sprintf("      Maximum restricted density: %.4f", max_leak))
    message(sprintf("      Maximum leakage depth: %g", max_depth))
  }
  message("   ----------------------------------\n")

  return(list(
    summary = summary_df,
    leaked_data = leaked_df
  ))
}
