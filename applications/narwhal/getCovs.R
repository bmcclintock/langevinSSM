#' @param buffer_meters distance (in meters) to expand the raster boundaries.
#' @param target_res cell resolution of the returned raster.
getCovs <- function(buffer_meters = 0, target_res = 821.408) {

  # 1. Base extent from original (IBCAO?) raster
  base_xmin <- 469606.2
  base_xmax <- 628137.9
  base_ymin <- 7968124
  base_ymax <- 8212082

  target_crs <- "EPSG:32617" # WGS 84 / UTM zone 17N

  buf_xmin <- base_xmin - buffer_meters
  buf_xmax <- base_xmax + buffer_meters
  buf_ymin <- base_ymin - buffer_meters
  buf_ymax <- base_ymax + buffer_meters

  template <- terra::rast(
    xmin = buf_xmin, xmax = buf_xmax,
    ymin = buf_ymin, ymax = buf_ymax,
    resolution = c(target_res, target_res),
    crs = target_crs
  )

  template_ll <- terra::project(template, "EPSG:4326")

  ext_ll <- as.vector(terra::ext(template_ll))

  # add a 0.5-degree buffer to the web query to prevent edge clipping during reprojection
  lon1 <- ext_ll["xmin"] - 0.5
  lon2 <- ext_ll["xmax"] + 0.5
  lat1 <- ext_ll["ymin"] - 0.5
  lat2 <- ext_ll["ymax"] + 0.5

  message(sprintf("Querying NOAA ETOPO Database for bounding box:\n Lon: [%.2f, %.2f], Lat: [%.2f, %.2f]", lon1, lon2, lat1, lat2))

  noaa_raw <- marmap::getNOAA.bathy(lon1 = lon1, lon2 = lon2,
                                    lat1 = lat1, lat2 = lat2,
                                    resolution = 1, keep = FALSE)

  bathy_ll <- terra::rast(marmap::as.raster(noaa_raw))

  message("Projecting NOAA data to target UTM grid...")

  bathy_utm <- terra::project(bathy_ll, template, method = "bilinear")
  names(bathy_utm) <- "bathy"

  mask_utm <- terra::ifel(bathy_utm > 0, 1, 0)
  names(mask_utm) <- "land_mask"

  return(list(
    bathy = bathy_utm,
    mask = mask_utm
  ))
}
