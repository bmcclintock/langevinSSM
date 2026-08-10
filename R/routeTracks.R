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

  message("Extracting restricted areas from raster mask...")

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

  message("Extracting and formatting observed points...")

  # Extract only the actual observations to route using base R to prevent dplyr from adding tbl_df classes
  obs_dat <- as.data.frame(data[!is.na(data[[x_col]]) & !is.na(data[[y_col]]), ])
  obs_dat <- obs_dat[order(obs_dat$id, obs_dat$date), ]

  # Build the visibility graph in the water
  message("Building visibility graph network around the barrier (this may take a moment)...")
  vis_graph <- pathroutr::prt_visgraph(barrier = nogo_poly)

  message("Rerouting barrier-crossing segments track-by-track...")

  track_ids <- unique(obs_dat$id)

  routed_list <- lapply(track_ids, function(trk_id) {
    message(sprintf("   Routing track: %s", trk_id))

    # Extract the specific track and explicitly cast to sf to guarantee class retention
    trk_df <- obs_dat[obs_dat$id == trk_id, ]
    trk_sf <- sf::st_as_sf(trk_df, coords = c(x_col, y_col), crs = sf::st_crs(nogo_poly))

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

      # Overwrite the points with the updated safe geometries (explicitly named arguments!)
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
  # Doing this here prevents column-mismatch crashes since prt_trim dropped rows.
  safe_coords <- sf::st_coordinates(track_routed_sf)
  obs_routed <- data.frame(
    id = track_routed_sf$id,
    date = track_routed_sf$date,
    mu.x_pr = safe_coords[, 1],
    mu.y_pr = safe_coords[, 2]
  )

  message("Interpolating latent state matrices for unobserved (NA) times...")

  # Prevent ".x" and ".y" suffix duplication if the user is re-routing already routed data
  if ("mu.x_pr" %in% names(data)) data$mu.x_pr <- NULL
  if ("mu.y_pr" %in% names(data)) data$mu.y_pr <- NULL

  # Join the safe points back to the FULL original dataset (including NA prediction times)
  # and interpolate the gaps along the new safe path.
  data_routed <- data %>%
    dplyr::left_join(obs_routed, by = c("id", "date")) %>%
    dplyr::arrange(id, date) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      # We use rule = 2 so any NAs at the absolute beginning or end of a track are carried forward/backward
      # stats::approx requires at least 2 non-NA points to interpolate, so we wrap it in a safety check.
      mu.x_pr = if (sum(!is.na(mu.x_pr)) >= 2) stats::approx(x = date, y = mu.x_pr, xout = date, rule = 2)$y else mu.x_pr,
      mu.y_pr = if (sum(!is.na(mu.y_pr)) >= 2) stats::approx(x = date, y = mu.y_pr, xout = date, rule = 2)$y else mu.y_pr
    ) %>%
    dplyr::ungroup()

  # Fallback to the original raw coordinates for any tracks/endpoints that lacked enough valid points to route
  missing_pr <- is.na(data_routed$mu.x_pr)
  if (any(missing_pr)) {
    data_routed$mu.x_pr[missing_pr] <- data_routed[[x_col]][missing_pr]
    data_routed$mu.y_pr[missing_pr] <- data_routed[[y_col]][missing_pr]
  }

  # Re-apply the dataLangevin class to ensure downstream compatibility
  data_routed <- as.data.frame(data_routed)
  class(data_routed) <- c("dataLangevin", "data.frame")
  attr(data_routed, "time.unit") <- attr(data, "time.unit")
  attr(data_routed, "coord") <- coord_cols

  message("Done! Appended 'mu.x_pr' and 'mu.y_pr' to the dataset.")

  return(data_routed)
}
