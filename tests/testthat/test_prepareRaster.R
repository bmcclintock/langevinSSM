# tests/testthat/test_prepareRaster.R

# --- Helper functions to generate mock objects for testing ---
get_valid_raster <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 100, 0, 100))
  terra::values(r) <- 1:100
  return(r)
}

get_valid_data <- function() {
  data.frame(
    id = rep("A", 3),
    date = as.POSIXct(c("2023-01-01 10:00:00", "2023-01-01 12:00:00", "2023-01-01 14:00:00"), tz = "UTC"),
    x = c(20, 50, 80),
    y = c(20, 50, 80)
  )
}

# ------------------------------------------------------------------

test_that("Input structure validation catches bad lists and names", {
  r <- get_valid_raster()

  # Not a list
  expect_error(prepareRaster(spatialCovs = r), "must be a list")

  # Not a named list
  expect_error(prepareRaster(spatialCovs = list(r)), "must be a named list")

  # Duplicate names
  expect_error(prepareRaster(spatialCovs = list(cov1 = r, cov1 = r)), "must have unique names")

  # Contains non-SpatRaster object
  expect_error(prepareRaster(spatialCovs = list(cov1 = r, cov2 = "not_a_raster")), "must be of class 'SpatRaster'")
})

# ------------------------------------------------------------------

test_that("Raster geometry matching catches misaligned rasters", {
  r1 <- get_valid_raster()

  # Create a raster with a different extent
  r2 <- terra::rast(nrows = 10, ncols = 10, ext = c(10, 110, 10, 110))
  terra::values(r2) <- 1:100

  expect_error(prepareRaster(spatialCovs = list(cov1 = r1, cov2 = r2)),
               "share the exact same projection \\(CRS\\), extent, and resolution")
})

# ------------------------------------------------------------------

test_that("Missing values in rasters are caught", {
  r <- get_valid_raster()
  terra::values(r)[5] <- NA # Inject an NA value

  expect_error(prepareRaster(spatialCovs = list(cov1 = r)),
               "missing values are not permitted")
})

# ------------------------------------------------------------------

test_that("Spatial overlap checks catch out-of-bounds data", {
  r <- get_valid_raster()
  covs <- list(cov1 = r)

  bad_data_total <- data.frame(x = c(150, 200), y = c(150, 200))

  expect_error(prepareRaster(spatialCovs = covs, data = bad_data_total, coord = c("x", "y")),
               "do not overlap with 'spatialCovs'")

  bad_data_partial <- data.frame(x = c(50, 150), y = c(50, 50))
  expect_error(prepareRaster(spatialCovs = covs, data = bad_data_partial, coord = c("x", "y")),
               "fall outside the boundaries of 'spatialCovs'")
})

# ------------------------------------------------------------------

test_that("Dynamic/Multi-layer raster time requirements are enforced", {
  r <- get_valid_raster()
  # Create a 2-layer raster
  r_multi <- c(r, r)
  covs <- list(cov1 = r_multi)
  valid_data <- get_valid_data()

  # Fails because the raster has no time values set
  expect_error(prepareRaster(spatialCovs = covs),
               "is a multi-layer raster that must have time values set")

  # Set time values on the raster
  terra::time(r_multi) <- as.POSIXct(c("2023-01-01", "2023-01-02"), tz = "UTC")
  covs_with_time <- list(cov1 = r_multi)

  # Fails because data doesn't have a 'date' column
  data_no_date <- valid_data
  data_no_date$date <- NULL
  expect_error(prepareRaster(spatialCovs = covs_with_time, data = data_no_date, coord = c("x", "y")),
               "requires a 'date' column in 'data'")

  # Succeeds when both raster times and data dates are present
  expect_type(prepareRaster(spatialCovs = covs_with_time, data = valid_data, coord = c("x", "y")), "list")
})

# ------------------------------------------------------------------

test_that("Column name conflicts between data and rasters are caught", {
  r <- get_valid_raster()
  valid_data <- get_valid_data()

  # Add a column to data that matches the raster name
  valid_data$habitat <- "forest"

  expect_error(prepareRaster(spatialCovs = list(habitat = r), data = valid_data, coord = c("x", "y")),
               "cannot have same names as data")
})

# ------------------------------------------------------------------

test_that("Valid inputs successfully build the C++ ready list structure", {
  r1 <- get_valid_raster()
  r2 <- get_valid_raster()
  covs <- list(cov1 = r1, cov2 = r2)
  valid_data <- get_valid_data()

  scale_factor <- 1000

  res <- prepareRaster(spatialCovs = covs, scaleFactor = scale_factor, data = valid_data, coord = c("x", "y"))

  # Check structure
  expect_type(res, "list")
  expect_true(all(c("raster_vals", "raster_coords", "raster_resolution",
                    "raster_extent", "n_covs", "all_z_values",
                    "n_zvals_cov", "cov_offset") %in% names(res)))

  # Check dimensions of permuted array: [ncol, nrow, nlayer]
  expect_equal(dim(res$raster_vals), c(10, 10, 2))

  # Check scale factor application on extent
  expected_extent <- as.vector(terra::ext(r1) / scale_factor)
  expect_equal(res$raster_extent, expected_extent)

  # Check single-layer times evaluate to 0
  expect_equal(res$all_z_values, c(0, 0))

  # Check covariate offsets
  expect_equal(res$n_covs, 2)
  expect_equal(res$n_zvals_cov, c(1, 1))
  expect_equal(res$cov_offset, c(0, 1))
})

# ------------------------------------------------------------------

test_that("Spatial overlap buffer warning catches data close to the edge", {
  r <- get_valid_raster() # Extent is 0 to 100
  covs <- list(cov1 = r)

  # Point is at x=90, y=50. It is physically inside the raster.
  # But x.sd is 5. The 3-sigma buffer is 15.
  # 90 + 15 = 105, which exceeds the max raster bound of 100.
  warn_data <- data.frame(
    id = "A",
    date = as.POSIXct("2023-01-01", tz="UTC"),
    x = 90,
    y = 50,
    x.sd = 5,
    y.sd = 1,
    smaj = NA
  )

  expect_warning(prepareRaster(spatialCovs = covs, data = warn_data, coord = c("x", "y")),
                 "close to the edge of 'spatialCovs' relative to their measurement error")
})
