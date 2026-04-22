# tests/testthat/test_simCov.R

test_that("simCov generates a SpatRaster with correct dimensions and extent", {
  # Use a small scale to speed up the test
  sca_val <- 10
  r <- simCov(sca = sca_val)

  # Check class
  expect_s4_class(r, "SpatRaster")

  # Check dimensions: should be (2*sca) x (2*sca)
  expected_dim <- 2 * sca_val
  expect_equal(terra::nrow(r), expected_dim)
  expect_equal(terra::ncol(r), expected_dim)
  expect_equal(terra::nlyr(r), 1)

  # Check extent: [-sca - 0.5, sca + 0.5]
  ext_r <- terra::ext(r)

  # FIX: Use as.numeric() to strip the names before checking equality
  expect_equal(as.numeric(ext_r$xmin), -sca_val)
  expect_equal(as.numeric(ext_r$xmax), sca_val)
  expect_equal(as.numeric(ext_r$ymin), -sca_val)
  expect_equal(as.numeric(ext_r$ymax), sca_val)
})

test_that("simCov generates sensible, finite numeric values", {
  r <- simCov(sca = 5, sigma2 = 0.5)
  vals <- terra::values(r)

  # Ensure valid data extraction
  expect_true(is.numeric(vals))
  expect_false(any(is.na(vals)))
  expect_true(all(is.finite(vals)))

  # Because it's a random field, all values should not be identical
  expect_gt(length(unique(as.vector(vals))), 1)
})

test_that("simCov is completely reproducible given a set.seed()", {
  set.seed(12345)
  r1 <- simCov(sca = 5, irange = 0.2, sigma2 = 0.1)

  set.seed(12345)
  r2 <- simCov(sca = 5, irange = 0.2, sigma2 = 0.1)

  # The rasters should be mathematically identical
  expect_equal(terra::values(r1), terra::values(r2))
})

test_that("simCov successfully handles custom FFT padding (M and N)", {
  # Force specific padding sizes
  r <- simCov(sca = 5, M = 64, N = 64)

  # The output raster should STILL be the original requested grid size (11x11),
  # as the padding is only used internally by fields::matern.image.cov
  expect_equal(terra::nrow(r), 10)
  expect_equal(terra::ncol(r), 10)
})

test_that("simCov fails gracefully on invalid inputs", {
  # Non-numeric scale
  expect_error(simCov(sca = "100"), "non-numeric")

  # FIX: Update expected error string and suppress the NaN warning
  # caused by negative log2 padding calculations
  expect_error(suppressWarnings(simCov(sca = -5)), "must be a non-negative number")

  # Zero or negative variance (sigma2) passed to sqrt() generates NaNs or errors
  expect_warning(simCov(sca = 5, sigma2 = -0.1), "NaNs produced")
})
