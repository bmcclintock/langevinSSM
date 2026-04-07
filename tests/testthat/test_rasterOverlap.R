# tests/testthat/test_rasterOverlap.R

# --- Helper functions for mock data ---
get_mock_raster <- function(vals = 1, extent = c(0, 10, 0, 10), res = 1) {
  r <- terra::rast(ext = extent, res = res)
  terra::values(r) <- rep(vals, length.out = terra::ncell(r))
  return(r)
}

test_that("rasterOverlap catches non-SpatRaster inputs", {
  r1 <- get_mock_raster()
  bad_input <- matrix(1, nrow = 10, ncol = 10)

  expect_error(rasterOverlap(r1, bad_input), "Both inputs must be terra::SpatRaster objects")
  expect_error(rasterOverlap(bad_input, r1), "Both inputs must be terra::SpatRaster objects")
})

test_that("rasterOverlap catches mismatched geometries", {
  r1 <- get_mock_raster(extent = c(0, 10, 0, 10))

  r2_bad_ext <- get_mock_raster(extent = c(5, 15, 5, 15))
  expect_error(rasterOverlap(r1, r2_bad_ext), "Rasters do not have the same geometry")

  r2_bad_res <- get_mock_raster(res = 0.5)
  expect_error(rasterOverlap(r1, r2_bad_res), "Rasters do not have the same geometry")
})

test_that("rasterOverlap handles negative values (log scale) with a warning", {
  r_pos <- get_mock_raster(vals = c(1, 2, 3, 4))
  r_neg <- get_mock_raster(vals = c(-1, -2, -3, -4))

  expect_warning(res1 <- rasterOverlap(r_neg, r_pos), "Negative values found in r1. Assuming log-scale")
  expect_warning(res2 <- rasterOverlap(r_pos, r_neg), "Negative values found in r2. Assuming log-scale")

  expect_true(is.numeric(res1))
  expect_true(is.numeric(res2))
})

test_that("rasterOverlap catches rasters that cannot be normalized", {
  r1 <- get_mock_raster(vals = 1)
  r_zero <- get_mock_raster(vals = 0)

  expect_error(rasterOverlap(r1, r_zero), "One or more layers in r1 or r2 sum to 0 or NA")

  r_na <- get_mock_raster(vals = NA)

  expect_error(rasterOverlap(r1, r_na), "One or more layers in r1 or r2 sum to 0 or NA")
})

test_that("rasterOverlap correctly calculates affinity for identical distributions", {
  r1 <- get_mock_raster(vals = runif(100, 0.1, 1))

  affinity <- rasterOverlap(r1, r1)

  expect_equal(as.numeric(affinity), 1, tolerance = 1e-6)
})

test_that("rasterOverlap correctly calculates affinity for entirely disjoint distributions", {
  vals1 <- c(rep(1, 50), rep(0, 50))
  vals2 <- c(rep(0, 50), rep(1, 50))

  r1 <- get_mock_raster(vals = vals1)
  r2 <- get_mock_raster(vals = vals2)

  affinity <- rasterOverlap(r1, r2)

  expect_equal(as.numeric(affinity), 0, tolerance = 1e-6)
})
