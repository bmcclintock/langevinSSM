# Extract and convert raster times
getRasterTimes <- function(r, time.unit) {
  t_vals <- terra::time(r)

  if (terra::nlyr(r) == 1 && (is.null(t_vals) || all(is.na(t_vals)))) {
    return(0)
  } else {
    if (inherits(t_vals, "POSIXt") || inherits(t_vals, "Date")) {
      return(as.numeric(difftime(t_vals, as.POSIXct("1970-01-01 00:00:00", tz = "UTC"), units = time.unit)))
    } else {
      return(as.numeric(t_vals))
    }
  }
}

#' @importFrom terra rast values nlyr time as.array crds res ext compareGeom
prepareRaster <- function(spatialCovs, scaleFactor=1, time.unit="hours", data = NULL, coord = NULL) {

  if(!is.list(spatialCovs)) stop('spatialCovs must be a list')
  spatialcovnames <- names(spatialCovs)
  if(is.null(spatialcovnames)) stop('spatialCovs must be a named list')
  nbSpatialCovs <- length(spatialcovnames)

  for(j in 1:nbSpatialCovs) {

    if(!inherits(spatialCovs[[j]], "SpatRaster")) {
      stop("spatialCovs$", spatialcovnames[j], " must be of class 'SpatRaster'")
    }

    if (j > 1) {
      if (!terra::compareGeom(spatialCovs[[1]], spatialCovs[[j]], stopOnError = FALSE)) {
        stop("All rasters in the 'spatialCovs' list must share the exact same projection (CRS), extent, and resolution. Mismatch detected in: ", spatialcovnames[j])
      }
    }

    if(any(is.na(terra::values(spatialCovs[[j]])))) {
      stop("missing values are not permitted in spatialCovs$", spatialcovnames[j])
    }

    if(terra::nlyr(spatialCovs[[j]]) > 1) {
      t_vals <- terra::time(spatialCovs[[j]])

      if(is.null(t_vals) || all(is.na(t_vals))) {
        stop("spatialCovs$", spatialcovnames[j], " is a multi-layer raster that must have time values set (see ?terra::time)")
      } else if(!is.null(data)) {
        if(!("date" %in% names(data))) {
          stop("spatialCovs$", spatialcovnames[j], " requires a 'date' column in 'data' to match the raster's dynamic layers")
        }

        d_vals <- data$date

        is_t_time <- inherits(t_vals, c("POSIXt", "Date"))
        is_d_time <- inherits(d_vals, c("POSIXt", "Date"))

        # --- Type Consistency Check ---
        if (is_t_time && !is_d_time) {
          stop("Type mismatch: spatialCovs$", spatialcovnames[j], " has POSIXt/Date time values, but data$date is numeric.")
        } else if (!is_t_time && is_d_time) {
          stop("Type mismatch: spatialCovs$", spatialcovnames[j], " has numeric time values, but data$date is POSIXt/Date.")
        } else if (!is_t_time && !is_d_time) {
          if (!is.numeric(t_vals) || !is.numeric(d_vals)) {
            stop("Time values must be either numeric or POSIXt/Date.")
          }
        }

        # --- Temporal Bounding Check ---
        min_t <- min(t_vals, na.rm = TRUE)
        max_t <- max(t_vals, na.rm = TRUE)
        min_d <- min(d_vals, na.rm = TRUE)
        max_d <- max(d_vals, na.rm = TRUE)

        if (min_d < min_t || max_d > max_t) {
          stop("The tracking data times fall outside the temporal boundaries of 'spatialCovs$", spatialcovnames[j], "'. Ensure min(data$date) and max(data$date) are strictly within the raster's time range.")
        }
      }
    }
  }

  rasterStack <- terra::rast(spatialCovs)

  if (!is.null(data) && all(coord %in% names(data))) {
    cov_ext <- as.vector(terra::ext(rasterStack))

    data_xmin <- min(data[[coord[1]]], na.rm = TRUE)
    data_xmax <- max(data[[coord[1]]], na.rm = TRUE)
    data_ymin <- min(data[[coord[2]]], na.rm = TRUE)
    data_ymax <- max(data[[coord[2]]], na.rm = TRUE)

    if (data_xmin > cov_ext["xmax"] || data_xmax < cov_ext["xmin"] ||
        data_ymin > cov_ext["ymax"] || data_ymax < cov_ext["ymin"]) {
      stop("The tracking data do not overlap with 'spatialCovs'. Please ensure they share the same metric projection and cover the same area.")
    }

    if (data_xmin < cov_ext["xmin"] || data_xmax > cov_ext["xmax"] ||
        data_ymin < cov_ext["ymin"] || data_ymax > cov_ext["ymax"]) {
      stop("Some tracking locations fall outside the boundaries of 'spatialCovs'. Expand the extent of the rasters.")
    }

    # 2. point-specific 3-sigma (99.7%) probabilistic bounding box (for edge warnings)
    err_x <- rep(0, nrow(data))
    err_y <- rep(0, nrow(data))

    # calculate exact marginal standard deviations from the observation error covariance matrix
    # (assuming neutral scaling parameters psi=1, tau=c(1,1) prior to model fitting)
    if (all(c("smaj", "smin", "eor") %in% names(data))) {
      kf_idx <- which(!is.na(data$smaj) & !is.na(data$smin) & !is.na(data$eor))
      if (length(kf_idx) > 0) {
        M2 <- (data$smaj[kf_idx]^2) / 2
        m2 <- (data$smin[kf_idx]^2) / 2
        s2c <- sin(data$eor[kf_idx])^2
        c2c <- cos(data$eor[kf_idx])^2

        err_x[kf_idx] <- sqrt(M2 * s2c + m2 * c2c)
        err_y[kf_idx] <- sqrt(M2 * c2c + m2 * s2c)
      }
    }

    if (all(c("x.err", "y.err") %in% names(data))) {
      ls_idx <- which(!is.na(data$x.err) & !is.na(data$y.err))
      if (length(ls_idx) > 0) {
        err_x[ls_idx] <- data$x.err[ls_idx]
        err_y[ls_idx] <- data$y.err[ls_idx]
      }
    }

    prob_xmin <- min(data[[coord[1]]] - 3 * err_x, na.rm = TRUE)
    prob_xmax <- max(data[[coord[1]]] + 3 * err_x, na.rm = TRUE)
    prob_ymin <- min(data[[coord[2]]] - 3 * err_y, na.rm = TRUE)
    prob_ymax <- max(data[[coord[2]]] + 3 * err_y, na.rm = TRUE)

    if (prob_xmin < cov_ext["xmin"] || prob_xmax > cov_ext["xmax"] ||
        prob_ymin < cov_ext["ymin"] || prob_ymax > cov_ext["ymax"]) {
      stop("Some tracking locations are dangerously close to the edge of 'spatialCovs' relative to their measurement error. Because the Langevin model estimates true locations (mu) that can deviate from observed coordinates, the optimizer will likely push these locations outside the raster extent during model fitting. The spatial extent of the rasters must be extended before proceeding.")
    }
  }

  if(!is.null(data) && any(spatialcovnames %in% names(data))) stop("spatialCovs cannot have same names as data")
  if(anyDuplicated(spatialcovnames)) stop("spatialCovs must have unique names")

  vals_array <- terra::as.array(rasterStack)

  # Permute to [ncol, nrow, nlayer] so it matches the C++ idx formula
  vals_array <- aperm(vals_array, c(2, 1, 3))

  times_list <- lapply(spatialCovs, getRasterTimes, time.unit=time.unit)

  n_zvals_cov <- sapply(spatialCovs, terra::nlyr)
  cov_offset_R <- c(0, cumsum(n_zvals_cov)[-length(n_zvals_cov)])

  all_z_values_R <- unlist(times_list)

  rasterData <- list(
    raster_vals = vals_array,
    raster_coords = terra::crds(rasterStack)/scaleFactor,
    raster_resolution = terra::res(rasterStack)/scaleFactor,
    raster_extent = as.vector(terra::ext(rasterStack)/scaleFactor),
    n_covs = length(n_zvals_cov),
    all_z_values = as.numeric(all_z_values_R), # The flattened raster slice times
    n_zvals_cov = as.integer(n_zvals_cov),
    cov_offset = as.integer(cov_offset_R)
  )
  return(rasterData)
}
