#' Prepare a spatial barrier for Langevin models
#'
#' Converts a barrier (a binary raster mask, or polygons) into the signed distance field (SDF) required for barrier constraints in the \code{langevinSSM} package. The resulting raster is tagged with a special attribute allowing downstream functions (\code{\link{fitLangevin}}, \code{\link{simLangevin}}, \code{\link{getUD}}) to automatically detect and apply the barrier penalty.
#'
#' @param x Either a \code{\link[terra]{SpatRaster-class}} object containing a binary mask (values of \code{1} indicate allowed movement areas, e.g. water, and values of \code{0} restricted areas, e.g. land), or polygons (an \code{sf}, \code{sfc} or \code{\link[terra]{SpatVector-class}} object) delimiting the restricted (or, with \code{restricted = FALSE}, the allowed) area.
#' @param template A \code{\link[terra]{SpatRaster-class}} defining the grid (extent, resolution and CRS) on which the SDF of polygons is computed. Required when \code{x} is polygons and ignored otherwise. Typically one of the habitat covariate rasters, so that the barrier shares their grid.
#' @param restricted Logical. When \code{x} is polygons, do they delimit the allowed area (\code{FALSE}, default) or the restricted area (\code{TRUE}; e.g. land polygons for a marine animal)?
#'
#' @details
#' The SDF \eqn{d} is the Euclidean distance to the boundary between the allowed and restricted areas, positive in the allowed area and negative in the restricted area. Its zero contour is the boundary itself, and \eqn{|\nabla d| = 1}. With \eqn{d} exact, the barrier force \eqn{-\lambda\min(d,0)\nabla d} is the penalization of Delporte et al. (2025), \eqn{-(x - \pi(x))/\lambda_D} with \eqn{\pi(x)} the projection onto the allowed area and \eqn{\lambda_D = 1/\lambda}.
#'
#' \strong{Polygons} (recommended when available, e.g. a vector coastline): the distance from every cell centre of \code{template} to the polygon boundary is computed exactly (\code{\link[sf]{st_distance}}), and the sign from polygon membership of the cell centre. The coastline is therefore represented as given, not as a raster staircase. \code{x} and \code{template} must share a projected CRS (distances in map units). Computation time grows with the number of cells and polygon vertices (roughly a minute per million cells for a detailed coastline).
#'
#' \strong{Raster mask}: the boundary is taken to be the cell edges between allowed and restricted cells. The distance between cell centres of opposite type is computed by \code{\link[terra]{distance}} and reduced by half a cell (\eqn{d = \mathrm{sign}(d_c)(|d_c| - r/2)}, with \eqn{r} the cell size), so that \eqn{d} is (approximately) the distance to that cell edge: \eqn{|d| = r/2} at the centres of the cells adjacent to the boundary and \eqn{|\nabla d| \approx 1} across it. Without the correction the field jumps by two cells across the boundary, which is a kink in the barrier force and an artificial steepening of it near the boundary; the interpolated zero contour (the cell edges) is the same either way. Along diagonal or irregular boundaries the cell-edge staircase is an approximation of the true boundary; supply polygons when they are available.
#'
#' @return A \code{\link[terra]{SpatRaster-class}} object (one layer named \code{"barrier"}) containing the signed distance field. Positive values represent distance into the allowed area; negative values represent distance into the restricted area. The object is assigned a \code{"barLangevin" = TRUE} attribute, and a \code{"barrierSource"} attribute (\code{"raster"} or \code{"polygons"}).
#'
#' @references
#' Delporte, A., Ditlevsen, S. and Samson, A., 2025. Spatial constraints improve filtering of measurement noise from animal tracks. arXiv preprint arXiv:2511.22430.
#'
#' @export
#' @importFrom stats na.omit
#' @importFrom terra distance unique nlyr crs "crs<-" classify res ncell xyFromCell setValues rast is.lonlat
prepBarrier <- function(x, template = NULL, restricted = FALSE) {

  if (inherits(x, c("sf", "sfc", "SpatVector"))) {
    sdf <- .sdf_polygons(x, template, restricted)
    src <- "polygons"
  } else if (inherits(x, "SpatRaster")) {
    sdf <- .sdf_mask(x)
    src <- "raster"
  } else {
    stop("'x' must be a terra::SpatRaster mask, or polygons (sf, sfc or terra::SpatVector).")
  }
  names(sdf) <- "barrier"

  # Tag the SpatRaster with our custom attribute
  attr(sdf, "barLangevin") <- TRUE
  attr(sdf, "barrierSource") <- src

  return(sdf)
}

# Signed distance to the cell-edge boundary of a binary mask (1 = allowed, 0 = restricted)
.sdf_mask <- function(maskRast) {
  if (terra::nlyr(maskRast) > 1) {
    stop("Barrier rasters must be static (single layer).")
  }

  unq_vals <- stats::na.omit(terra::unique(maskRast)[[1]])
  if (!all(unq_vals %in% c(0, 1))) {
    stop("Barrier raster must only contain binary values (0 and 1).")
  }
  if (!all(c(0, 1) %in% unq_vals)) {
    stop("Barrier raster must contain both allowed (1) and restricted (0) cells.")
  }
  r <- terra::res(maskRast)
  if (abs(r[1] - r[2]) > 1e-6 * max(r)) {
    stop("Barrier raster must have square cells (equal x and y resolution).")
  }

  if (.is_lonlat(maskRast)) {
    warning("The barrier mask has a longitude/latitude CRS (terra's default for rasters created without a CRS whose extent lies within +-180/+-90): distances are computed in metres, not in the units of the coordinates. Use a projected CRS (or crs = \"local\" for unitless coordinates).", call. = FALSE)
  }
  orig_crs <- terra::crs(maskRast)

  # If no CRS is provided (e.g., simulated data), temporarily assign "local"
  # to prevent terra's internal C++ warning about unknown CRSs
  if (orig_crs == "") {
    terra::crs(maskRast) <- "local"
  }

  message("   Calculating signed distance field for barrier...")

  # distance() of a raster with NA cells gives, for each NA cell, the distance to the nearest non-NA
  # cell centre (and 0 for non-NA cells)
  d_allowed <- terra::distance(terra::classify(maskRast, cbind(1, NA)))     # allowed cells -> nearest restricted
  d_restricted <- terra::distance(terra::classify(maskRast, cbind(0, NA)))  # restricted cells -> nearest allowed
  sdf <- d_allowed - d_restricted                                           # positive allowed, negative prohibited
  # centre-to-centre -> centre-to-cell-edge
  sdf <- sign(sdf) * (abs(sdf) - r[1] / 2)

  terra::crs(sdf) <- orig_crs
  sdf
}

# Exact signed distance from the cell centres of 'template' to the boundary of polygons
.sdf_polygons <- function(x, template, restricted) {
  if (is.null(template) || !inherits(template, "SpatRaster")) {
    stop("'template' (a terra::SpatRaster defining the grid) is required when 'x' is polygons.")
  }
  if (!is.logical(restricted) || length(restricted) != 1 || is.na(restricted)) stop("'restricted' must be TRUE or FALSE.")
  poly <- if (inherits(x, "SpatVector")) sf::st_as_sf(x) else sf::st_as_sf(sf::st_geometry(x))
  gtype <- unique(as.character(sf::st_geometry_type(poly)))
  if (!all(gtype %in% c("POLYGON", "MULTIPOLYGON"))) stop("'x' must contain (multi)polygons.")
  tcrs <- terra::crs(template)
  if (tcrs == "" || is.na(sf::st_crs(poly))) stop("'x' and 'template' must both have a (projected) CRS.")
  if (.is_lonlat(template)) stop("'template' must have a projected CRS (distances are computed in map units).")
  if (sf::st_crs(poly) != sf::st_crs(tcrs)) poly <- sf::st_transform(poly, sf::st_crs(tcrs))
  if (isTRUE(sf::st_is_longlat(poly))) stop("'x' must have a projected CRS.")

  message("   Calculating signed distance field for barrier from polygons...")
  g <- sf::st_union(sf::st_make_valid(sf::st_geometry(poly)))
  boundary <- sf::st_boundary(g)

  tmpl <- terra::rast(template[[1]])
  xy <- terra::xyFromCell(tmpl, seq_len(terra::ncell(tmpl)))
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = sf::st_crs(g))
  inside <- lengths(sf::st_intersects(pts, g)) > 0
  dist <- numeric(nrow(xy))
  for (ch in split(seq_len(nrow(xy)), ceiling(seq_len(nrow(xy)) / 50000))) {
    dist[ch] <- as.numeric(sf::st_distance(pts[ch, ], boundary))
  }
  allowed <- if (restricted) !inside else inside
  if (all(allowed) || !any(allowed)) stop("The polygons do not split the template grid into allowed and restricted cells.")
  terra::setValues(tmpl, ifelse(allowed, dist, -dist))
}

# TRUE only when the raster has a known lon/lat CRS. A missing CRS is not
# treated as lon/lat, and terra's "unknown crs" warning is not passed on.
.is_lonlat <- function(r) {
  cr <- terra::crs(r)
  if (is.na(cr) || !nzchar(cr)) return(FALSE)
  isTRUE(suppressWarnings(terra::is.lonlat(r)))
}
