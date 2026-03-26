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

#' @importFrom terra rast values nlyr time as.array crds res
prepareRaster <- function(spatialCovs, scaleFactor=1, time.unit="hours", data = NULL) {

  if(!is.list(spatialCovs)) stop('spatialCovs must be a list')
  spatialcovnames <- names(spatialCovs)
  if(is.null(spatialcovnames)) stop('spatialCovs must be a named list')
  nbSpatialCovs <- length(spatialcovnames)

  #if (!requireNamespace("terra", quietly = TRUE)) {
  #  stop("Package \"terra\" needed for spatial covariates. Please install it.", call. = FALSE)
  #}

  rasterStack <- terra::rast(spatialCovs)

  for(j in 1:nbSpatialCovs) {
    if(!inherits(spatialCovs[[j]], "SpatRaster")) {
      stop("spatialCovs$", spatialcovnames[j], " must be of class 'SpatRaster'")
    }

    if(any(is.na(terra::values(spatialCovs[[j]])))) {
      stop("missing values are not permitted in spatialCovs$", spatialcovnames[j])
    }

    if(terra::nlyr(spatialCovs[[j]]) > 1) {

      t_vals <- terra::time(spatialCovs[[j]])

      # Check if time/Z values are set
      if(is.null(t_vals) || all(is.na(t_vals))) {
        stop("spatialCovs$", spatialcovnames[j], " is a multi-layer raster that must have time values set (see ?terra::time)")
      }

      else if(!is.null(data) && !("date" %in% names(data))) {
        stop("spatialCovs$", spatialcovnames[j], " requires a 'date' column in 'data' to match the raster's dynamic layers")
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
