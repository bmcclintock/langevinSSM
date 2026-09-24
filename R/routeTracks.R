#' Reroute tracks around barriers using \code{pathroutr}
#'
#' This function converts a binary barrier mask into an \code{sf} polygon and uses the \code{pathroutr} package to route observed coordinates around the barrier. It returns the original \code{dataLangevin} object with routed coordinates (\code{mu.x_pr}, \code{mu.y_pr}) appended, which can then be passed to \code{\link{fitLangevin}} as initial values for the latent locations (\code{mu}).
#'
#' @param data A formatted \code{dataLangevin} object returned by \code{\link{formatData}}.
#' @param maskRast A \code{\link[terra]{SpatRaster-class}} object containing a binary mask. Values of \code{1} indicate allowed movement areas (e.g., water), and values of \code{0} indicate restricted areas (e.g., land).
#'
#' @details
#' \code{\link{fitLangevin}} can sometimes struggle to push tracks around complex coastlines. \code{routeTracks} builds a visibility network graph inside the allowed areas (\code{maskRast = 1}) and bends any invalid track segments around the barrier. The resulting coordinates are appended to the dataset as \code{mu.x_pr} and \code{mu.y_pr}.
#'
#' @return The original \code{dataLangevin} data frame, appended with two new columns: \code{mu.x_pr} and \code{mu.y_pr}. These represent the re-routed tracks.
#'
#' @examples
#' if (requireNamespace("pathroutr", quietly = TRUE)) {
#'
#'   library(ggplot2)
#'
#'   # format the tracking data
#'   formatDat <- formatData(unformatDat, time.unit = "hours")
#'
#'   # create a dummy mask (e.g., 1 = water, 0 = land)
#'   mask_rast <- exampleCovs[[1]]
#'   coords <- terra::crds(mask_rast)
#'   terra::values(mask_rast) <- ifelse(coords[, "x"] >= 1010, 1, 0)
#'   names(mask_rast) <- "coast_barrier"
#'
#'   # route the tracks around barrier
#'   formatDat_routed <- routeTracks(formatDat, maskRast = mask_rast)
#'
#'   barrier <- prepBarrier(mask_rast)
#'   plot(formatDat, spatialCovs = list(coast_barrier=barrier),
#'        maskRast = mask_rast)$coast_barrier +
#'   ggplot2::geom_path(data = formatDat_routed,
#'                      aes(x = mu.x_pr, y = mu.y_pr, group = id, color = as.factor(id)),
#'                      linewidth = 0.8) +
#'   ggplot2::geom_point(data = formatDat_routed,
#'                       aes(x = mu.x_pr, y = mu.y_pr, color = as.factor(id)),
#'                       shape = 16, size = 1.5) +
#'   ggplot2::scale_color_manual(
#'     name = "Track ID",
#'     values = c(
#'       "1" = "yellow",
#'       "2" = "lightblue",
#'       "3" = "orange",
#'       "Observed" = "lightgrey"
#'     )
#'   )
#' }
#'
#' @export
routeTracks <- function(data, maskRast) {

  if (!requireNamespace("pathroutr", quietly = TRUE)) {
    stop("The 'pathroutr' package is required to use this function. You can install it from GitHub using remotes::install_github('jmlondon/pathroutr').")
  }

  if (!"package:sf" %in% search()) {
    do.call("require", list(package = "sf", quietly = TRUE, character.only = TRUE))
  }
  if (!"package:dplyr" %in% search()) {
    do.call("require", list(package = "dplyr", quietly = TRUE, character.only = TRUE))
  }

  if (!inherits(data, "dataLangevin")) {
    stop("'data' must be a 'dataLangevin' object. Please format your data using formatData() first.")
  }
  if (!inherits(maskRast, "SpatRaster")) {
    stop("'maskRast' must be a terra::SpatRaster object.")
  }

  coord_cols <- attr(data, "coord")
  if (is.null(coord_cols)) coord_cols <- c("x", "y")

  x_col <- coord_cols[1]
  y_col <- coord_cols[2]

  message("   Extracting restricted areas from raster mask...")

  # Isolate the restricted areas (0) by setting allowed areas (1) to NA.
  nogo_rast <- terra::ifel(maskRast == 0, 1, NA)

  # Convert the restricted raster cells to an sf polygon and dissolve internal borders
  nogo_poly <- terra::as.polygons(nogo_rast) %>%
    sf::st_as_sf() %>%
    sf::st_union()

  # DENSIFY POLYGON: Add vertices along perfectly straight raster edges so pathroutr
  # has valid nodes to snap to, rather than firing points into the extreme corners
  poly_res <- max(terra::res(maskRast))
  nogo_poly <- sf::st_segmentize(nogo_poly, dfMaxLength = poly_res)

  has_missing <- any(is.na(data[[x_col]])) | any(is.na(data[[y_col]]))

  if (has_missing) {
    message("   Pre-interpolating missing locations to route the entire track...")

    # Interpolate NA coordinates first so they can be pushed through pathroutr alongside observations
    full_dat <- data %>%
      dplyr::arrange(id, date) %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(
        temp_x = if (sum(!is.na(.data[[x_col]])) >= 2) stats::approx(x = date, y = .data[[x_col]], xout = date, rule = 2)$y else .data[[x_col]],
        temp_y = if (sum(!is.na(.data[[y_col]])) >= 2) stats::approx(x = date, y = .data[[y_col]], xout = date, rule = 2)$y else .data[[y_col]]
      ) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  } else {
    full_dat <- data
    full_dat$temp_x <- full_dat[[x_col]]
    full_dat$temp_y <- full_dat[[y_col]]
  }

  full_dat <- full_dat[!is.na(full_dat$temp_x) & !is.na(full_dat$temp_y), ]
  full_dat <- full_dat[order(full_dat$id, full_dat$date), ]

  # Identify which points fall inside the barrier
  full_sf <- sf::st_as_sf(full_dat, coords = c("temp_x", "temp_y"), crs = sf::st_crs(nogo_poly))
  in_barrier <- sf::st_intersects(full_sf, nogo_poly, sparse = FALSE)[, 1]

  # Snap ALL points (observed or interpolated) that fall on land to the coastline
  if (any(in_barrier)) {
    message("   Snapping land-bound points to the barrier boundary before routing...")

    nogo_boundary <- sf::st_cast(nogo_poly, "MULTILINESTRING")
    bad_pts <- full_sf[in_barrier, ]
    nearest_lines <- sf::st_nearest_points(bad_pts, nogo_boundary)
    nearest_pts <- sf::st_cast(nearest_lines, "POINT")
    bound_pts <- nearest_pts[seq(2, length(nearest_pts), by = 2)]

    safe_coords_pre <- sf::st_coordinates(bound_pts)
    full_dat$temp_x[in_barrier] <- safe_coords_pre[, 1]
    full_dat$temp_y[in_barrier] <- safe_coords_pre[, 2]
  }

  # Build the visibility graph in the water
  message("   Building visibility graph network around the barrier (this may take a moment)...")
  vis_graph <- pathroutr::prt_visgraph(barrier = nogo_poly)

  message("   Rerouting barrier-crossing segments track-by-track...")

  track_ids <- unique(full_dat$id)

  routed_list <- lapply(track_ids, function(trk_id) {
    message(sprintf("      Routing track: %s", trk_id))

    # Extract the specific track and explicitly cast to sf to guarantee class retention
    # We map the newly interpolated and snapped coordinates (temp_x, temp_y)
    trk_df <- full_dat[full_dat$id == trk_id, ]
    trk_sf <- sf::st_as_sf(trk_df, coords = c("temp_x", "temp_y"), crs = sf::st_crs(nogo_poly))

    # TRIM: pathroutr cannot route a track if it begins or ends inside the barrier.
    # We wrap this in try() because prt_trim will crash if the ENTIRE track is on land.
    trimmed_sf <- try(pathroutr::prt_trim(trkpts = trk_sf, barrier = nogo_poly), silent = TRUE)

    # If trimming crashes or leaves < 2 points, skip routing for this track entirely
    if (inherits(trimmed_sf, "try-error") || nrow(trimmed_sf) < 2) {
      warning(sprintf("Track %s is entirely on land or lacks sufficient valid points. Skipping routing for this track.", trk_id), call. = FALSE)
    } else {
      trk_sf <- trimmed_sf

      # Calculate the shortest path routing
      routes <- pathroutr::prt_reroute(trkpts = trk_sf, barrier = nogo_poly, vis_graph = vis_graph)

      # Overwrite the points with the updated safe geometries
      if (nrow(routes) > 0) {
        trk_sf <- pathroutr::prt_update_points(rrt_pts = routes, trkpts = trk_sf)
      }
    }

    return(trk_sf)
  })

  # Combine back into a single sf object (bind_rows safely handles mismatched columns if pathroutr added any)
  track_routed_sf <- dplyr::bind_rows(routed_list)
  if (!inherits(track_routed_sf, "sf")) {
    track_routed_sf <- sf::st_as_sf(track_routed_sf)
  }

  # Extract the safely routed coordinates from the spatial object.
  safe_coords <- sf::st_coordinates(track_routed_sf)
  obs_routed <- data.frame(
    id = track_routed_sf$id,
    date = track_routed_sf$date,
    mu.x_pr = safe_coords[, 1],
    mu.y_pr = safe_coords[, 2],
    pr_idx = 1:nrow(track_routed_sf) # Keep sequence to track inserted waypoints
  )

  message("   Mapping routed points to dataset...")

  # Prevent suffix duplication if the user is re-routing already routed data
  if ("mu.x_pr" %in% names(data)) data$mu.x_pr <- NULL
  if ("mu.y_pr" %in% names(data)) data$mu.y_pr <- NULL

  # Join original data to obs_routed using full_join to retain pathroutr's inserted rows
  data_routed <- data %>%
    dplyr::full_join(obs_routed, by = c("id", "date")) %>%
    dplyr::arrange(id, date, pr_idx)

  # Flag inserted waypoints. pathroutr duplicates the attributes of the segment's starting node.
  # For any given id & date block, the first row is the original node; subsequent rows are inserted.
  data_routed <- data_routed %>%
    dplyr::group_by(id, date) %>%
    dplyr::mutate(is_inserted = dplyr::row_number() > 1) %>%
    dplyr::ungroup()

  # Wipe out observation data for inserted waypoints (they are latent padding locations, not observations)
  cols_to_na <- setdiff(names(data), c("id", "date", "dt"))
  for(col in cols_to_na) {
    if(col %in% names(data_routed)) {
      data_routed[data_routed$is_inserted, col] <- NA
    }
  }

  # Bring in temp_x and temp_y as fallback for mu.x_pr / mu.y_pr
  data_routed <- data_routed %>%
    dplyr::left_join(dplyr::select(full_dat, id, date, temp_x, temp_y), by = c("id", "date")) %>%
    dplyr::mutate(
      mu.x_pr = ifelse(is.na(mu.x_pr), temp_x, mu.x_pr),
      mu.y_pr = ifelse(is.na(mu.y_pr), temp_y, mu.y_pr)
    )

  # Fallback to the original raw coordinates for endpoints that lacked valid points to route
  missing_pr <- is.na(data_routed$mu.x_pr)
  if (any(missing_pr)) {
    data_routed$mu.x_pr[missing_pr] <- data_routed[[x_col]][missing_pr]
    data_routed$mu.y_pr[missing_pr] <- data_routed[[y_col]][missing_pr]
  }

  # Clean up temporary columns
  data_routed$temp_x <- NULL
  data_routed$temp_y <- NULL
  data_routed$pr_idx <- NULL

  # Interpolate timestamps for inserted waypoints
  data_routed$date_num <- as.numeric(data_routed$date)
  data_routed$date_num[data_routed$is_inserted] <- NA

  data_routed <- data_routed %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      date_num = stats::approx(x = 1:dplyr::n(), y = date_num, xout = 1:dplyr::n(), rule = 2)$y
    ) %>%
    dplyr::ungroup()

  # Restore date format
  if (inherits(data$date, "POSIXt")) {
    tz_attr <- attr(data$date, "tzone")
    if (is.null(tz_attr)) tz_attr <- "UTC"
    data_routed$date <- as.POSIXct(data_routed$date_num, origin = "1970-01-01", tz = tz_attr)
  } else if (inherits(data$date, "Date")) {
    data_routed$date <- as.Date(data_routed$date_num, origin = "1970-01-01")
  } else {
    data_routed$date <- data_routed$date_num
  }

  data_routed$date_num <- NULL
  data_routed$is_inserted <- NULL

  # Recalculate dt
  time.unit <- attr(data, "time.unit")
  is_numeric_date <- is.numeric(data_routed$date)

  data_routed <- data_routed %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      dt = if (is_numeric_date) {
        c(0, as.numeric(diff(date)))
      } else {
        c(0, as.numeric(difftime(date[-1], date[-dplyr::n()], units = time.unit)))
      }
    ) %>%
    dplyr::ungroup()

  # Re-apply the dataLangevin class to ensure downstream compatibility
  data_routed <- as.data.frame(data_routed)
  class(data_routed) <- c("dataLangevin", "data.frame")
  attr(data_routed, "time.unit") <- time.unit
  attr(data_routed, "coord") <- coord_cols

  message("   Done. Appended 'mu.x_pr' and 'mu.y_pr' to the dataset.")

  return(data_routed)
}
