# --- prepBarrier (raster half-cell correction, polygons) and the smooth barrier onset ---

test_that("prepBarrier on a raster mask gives the distance to the cell-edge boundary", {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 1000, 0, 1000), crs = "local")
  terra::values(r) <- ifelse(terra::xFromCell(r, seq_len(terra::ncell(r))) > 400, 1, 0)   # boundary at x = 400
  sdf <- suppressMessages(prepBarrier(r))
  expect_true(isTRUE(attr(sdf, "barLangevin")))
  expect_identical(attr(sdf, "barrierSource"), "raster")
  xc <- terra::xFromCell(sdf, seq_len(terra::ncell(sdf)))
  expect_equal(as.numeric(terra::values(sdf)), xc - 400, tolerance = 1e-8)
  expect_false(anyNA(terra::values(sdf)))
})

test_that("prepBarrier validates the mask", {
  r <- terra::rast(nrows = 5, ncols = 5, ext = c(0, 500, 0, 500), crs = "local")
  terra::values(r) <- 1
  expect_error(suppressMessages(prepBarrier(r)), "both allowed")
  terra::values(r) <- 2
  expect_error(suppressMessages(prepBarrier(r)), "binary")
  expect_error(prepBarrier(1), "SpatRaster")
})

test_that("prepBarrier on polygons gives the exact signed distance", {
  skip_if_not_installed("sf")
  tmpl <- terra::rast(nrows = 20, ncols = 20, ext = c(0, 20, 0, 20), crs = "EPSG:32616")
  land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(c(-50, -50), c(7.3, -50), c(7.3, 70), c(-50, 70), c(-50, -50)))),
                                          crs = 32616))
  sdf <- suppressMessages(prepBarrier(land, template = tmpl, restricted = TRUE))
  expect_identical(attr(sdf, "barrierSource"), "polygons")
  xc <- terra::xFromCell(sdf, seq_len(terra::ncell(sdf)))
  expect_equal(as.numeric(terra::values(sdf)), xc - 7.3, tolerance = 1e-8)
  # restricted = FALSE: the polygon is the allowed area
  sdf2 <- suppressMessages(prepBarrier(land, template = tmpl, restricted = FALSE))
  expect_equal(as.numeric(terra::values(sdf2)), 7.3 - xc, tolerance = 1e-8)
  expect_error(prepBarrier(land), "template")
})

test_that(".barrier_potential is C4 at the boundary and at the onset length, and 0 recovers the quadratic", {
  B <- langevinSSM:::.barrier_potential
  d <- seq(-5, 3, by = 0.01)
  expect_equal(B(d, 0), ifelse(d < 0, 0.5 * d^2, 0))
  l <- 2
  f <- function(u) B(-u, l)   # u = depth into the restricted zone
  h <- 1e-3
  # value continuity at the onset length (inner polynomial vs outer quadratic)
  expect_equal(f(l), l^2 / 7)
  expect_equal(l^2 * (1/7 - 1/2 + 1/2), 0.5 * l^2 - 0.5 * l * l + l^2 / 7)
  # one-sided first and second derivatives agree at u = 0 and u = l
  d1 <- function(u, s) s * (f(u + s * h) - f(u)) / h
  d2 <- function(u, s) (f(u + 2 * s * h) - 2 * f(u + s * h) + f(u)) / h^2
  expect_equal(d1(0, 1), 0, tolerance = 1e-6); expect_equal(d2(0, 1), 0, tolerance = 1e-4)
  expect_equal(d1(l, -1), d1(l, 1), tolerance = 1e-2)
  expect_equal(d2(l, -1), d2(l, 1), tolerance = 1e-2)
  expect_equal(d2(l, 1), 1, tolerance = 1e-2)
  # beyond the onset: force = u - l/2
  u <- c(3, 4, 6)
  expect_equal((f(u + h) - f(u - h)) / (2 * h), u - l / 2, tolerance = 1e-6)
  # allowed side unchanged
  expect_true(all(B(d[d >= 0], l) == 0))
})

test_that(".resolve_onset handles 'res', numbers and bad input", {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 50, 0, 50))
  expect_equal(langevinSSM:::.resolve_onset("res", r), 5)
  expect_equal(langevinSSM:::.resolve_onset(NULL, r), 0)
  expect_equal(langevinSSM:::.resolve_onset(12.5, r), 12.5)
  expect_error(langevinSSM:::.resolve_onset(-1, r), "non-negative")
  expect_error(langevinSSM:::.resolve_onset("cell", r), "res")
})

make_onset_covs <- function() {
  r <- terra::rast(nrows = 40, ncols = 40, ext = c(-1000, 1000, -1000, 1000))
  terra::values(r) <- ifelse(terra::xFromCell(r, seq_len(terra::ncell(r))) >= -600, 1, 0)
  bar <- suppressMessages(prepBarrier(r))
  cov1 <- r
  terra::values(cov1) <- terra::yFromCell(r, seq_len(terra::ncell(r))) / 1000
  names(cov1) <- "cov1"
  list(cov1 = cov1, coast_barrier = bar)
}

test_that("simLangevin stores barrierOnset and fitLangevin picks it up and stores it", {
  covs <- make_onset_covs()
  set.seed(1)
  sim <- suppressMessages(simLangevin(model = "overdamped", spatialCovs = covs, barrier = "coast_barrier",
                                      par = list(beta = 0.5, sigma = 20), lambda = 0.01, barrierOnset = "res",
                                      timeStep = 1, obsPerAnimal = 60,
                                      initialPosition = matrix(c(-500, 0), 1)))
  expect_equal(attr(sim, "barrierOnset"), 50)
  fit <- suppressWarnings(suppressMessages(fitLangevin(sim, model = "overdamped", spatialCovs = covs,
                                                       par = list(beta = 0.5, sigma = 20))))
  expect_equal(fit$conditions$barrierOnset, 50)
})

test_that("the onset changes the likelihood only through restricted-side locations", {
  covs <- make_onset_covs()
  set.seed(2)
  sim <- suppressMessages(simLangevin(model = "overdamped", spatialCovs = covs, barrier = "coast_barrier",
                                      par = list(beta = 0.5, sigma = 20), lambda = 0.01,
                                      timeStep = 1, obsPerAnimal = 40,
                                      initialPosition = matrix(c(200, 0), 1)))
  ok <- all(sim$mu.x > -600 + 2 * 50 + 1)   # never within onset + 1 cell of land
  skip_if_not(ok, "simulated track approached the barrier")
  f0 <- suppressWarnings(suppressMessages(fitLangevin(sim, model = "overdamped", spatialCovs = covs, lambda = 0.01,
                                                      barrierOnset = 0, par = list(beta = 0.5, sigma = 20))))
  f1 <- suppressWarnings(suppressMessages(fitLangevin(sim, model = "overdamped", spatialCovs = covs, lambda = 0.01,
                                                      barrierOnset = 50, par = list(beta = 0.5, sigma = 20))))
  expect_equal(f0$objective, f1$objective, tolerance = 1e-6)
  expect_equal(f1$conditions$barrierOnset, 50)
})
