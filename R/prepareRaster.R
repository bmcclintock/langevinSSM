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
      } else if(!is.null(data) && !("date" %in% names(data))) {
        stop("spatialCovs$", spatialcovnames[j], " requires a 'date' column in 'data' to match the raster's dynamic layers")
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

    err_x_vals <- c(data$x.sd, data$smaj)
    err_y_vals <- c(data$y.sd, data$smaj)

    max_err_x <- if (all(is.na(err_x_vals))) 0 else max(err_x_vals, na.rm = TRUE)
    max_err_y <- if (all(is.na(err_y_vals))) 0 else max(err_y_vals, na.rm = TRUE)

    # Use a 3-sigma (approx 99% CI) buffer
    buffer_x <- 3 * max_err_x
    buffer_y <- 3 * max_err_y

    if (max_err_x > 0 || max_err_y > 0) {
      if ((data_xmin - buffer_x) < cov_ext["xmin"] || (data_xmax + buffer_x) > cov_ext["xmax"] ||
          (data_ymin - buffer_y) < cov_ext["ymin"] || (data_ymax + buffer_y) > cov_ext["ymax"]) {
        warning("Some tracking locations are close to the edge of 'spatialCovs' relative to their measurement error. Because the Langevin model estimates true locations (mu) that can deviate from observed coordinates, the model may attempt to push locations outside the raster extent during fitting, causing convergence failures or crashes. Consider expanding the spatial extent of your rasters.")
      }
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
