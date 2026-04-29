# tests/testthat/test_checkBarrier.R

# --- Helper functions for mock data ---
get_mock_covs <- function() {
  # Large raster to safely contain 3-sigma error buffers: -100 to 100
  r <- terra::rast(nrows = 20, ncols = 20, ext = terra::ext(-100, 100, -100, 100))
  terra::values(r) <- runif(400)
  names(r) <- "habitat"

  # Barrier: Left half (x < 0) is allowed (1), Right half (x >= 0) is restricted (0)
  barrier <- r
  terra::values(barrier) <- ifelse(terra::crds(barrier)[, "x"] < 0, 1, 0)
  names(barrier) <- "coast_barrier"

  list(habitat = r, coast_barrier = barrier)
}

get_mock_track <- function() {
  # 20 points, all safely in the allowed zone (x = -50 to -10)
  df <- data.frame(
    id = "A",
    date = as.POSIXct("2023-01-01 00:00:00", tz = "UTC") + (1:20) * 3600,
    x = seq(-50, -10, length.out = 20),
    y = rep(0, 20),
    smaj = 1,
    smin = 0.5,
    eor = 45,
    lc = "G"
  )
  # Suppress projection warnings for mock data
  suppressWarnings(formatData(df, time.unit = "hours"))
}

# --- Tests ---

test_that("checkBarrier validates inputs correctly", {
  covs <- get_mock_covs()
  data <- get_mock_track()

  # 1. Not a dataLangevin object
  expect_error(
    suppressWarnings(suppressMessages(capture.output(checkBarrier(data = as.data.frame(data), spatialCovs = covs, barrier = "coast_barrier")))),
    "'data' must be a dataLangevin object"
  )

  # 2. Barrier name doesn't exist
  expect_error(
    suppressWarnings(suppressMessages(capture.output(checkBarrier(data = data, spatialCovs = covs, barrier = "non_existent_barrier")))),
    "The 'barrier' name must exist in 'spatialCovs'"
  )
})

test_that("checkBarrier works cleanly when no locations are in restricted zone", {
  covs <- get_mock_covs()
  data <- get_mock_track() # All points x <= -10 (safe)

  captured <- capture.output({
    res <- suppressWarnings(suppressMessages(checkBarrier(data, spatialCovs = covs, barrier = "coast_barrier")))
  })

  expect_type(res, "list")
  expect_null(res$implausible_locations)
  expect_equal(res$lambda_min, 0)
  expect_equal(res$Y_max, 0)

  # Data should remain completely unmodified
  expect_identical(res$filtered_data$x, data$x)
  expect_identical(res$filtered_data$smaj, data$smaj)
})

test_that("checkBarrier handles plausible restricted locations correctly", {
  covs <- get_mock_covs()
  data <- get_mock_track()

  # Push row 10 slightly into the restricted zone (x = 2, so 2 units inland)
  # Give it a massive error (smaj = 10) so it's statistically highly plausible to be in the water
  data$x[10] <- 2
  data$smaj[10] <- 10
  data$smin[10] <- 10

  captured <- capture.output({
    res <- suppressWarnings(suppressMessages(checkBarrier(data, spatialCovs = covs, barrier = "coast_barrier")))
  })

  # It should realize the point is recoverable
  expect_null(res$implausible_locations)
  expect_true(res$lambda_min > 0)
  expect_true(res$Y_max > 0)

  # It should NOT set the point to NA
  expect_false(is.na(res$filtered_data$x[10]))
})

test_that("checkBarrier identifies mathematically implausible locations and sets them to NA", {
  covs <- get_mock_covs()
  data <- get_mock_track()

  # Push row 15 deep into the restricted zone (x = 10, so 10 units inland)
  # Give it a tiny error (smaj = 0.01). Z-score will be huge.
  data$x[15] <- 10
  data$smaj[15] <- 0.01
  data$smin[15] <- 0.01

  captured <- capture.output({
    res <- suppressWarnings(suppressMessages(checkBarrier(data, spatialCovs = covs, barrier = "coast_barrier")))
  })

  # It should flag the point
  expect_false(is.null(res$implausible_locations))
  expect_equal(nrow(res$implausible_locations), 1)
  expect_equal(res$implausible_locations$row_index, 15)

  # It should set the coordinates and errors to NA
  expect_true(is.na(res$filtered_data$x[15]))
  expect_true(is.na(res$filtered_data$y[15]))
  expect_true(is.na(res$filtered_data$smaj[15]))
  expect_true(is.na(res$filtered_data$smin[15]))
  expect_true(is.na(res$filtered_data$eor[15]))

  # But the surrounding points should remain untouched
  expect_false(is.na(res$filtered_data$x[14]))
  expect_false(is.na(res$filtered_data$x[16]))
})

test_that("checkBarrier throws error if restricted point lacks measurement error", {
  covs <- get_mock_covs()
  data <- get_mock_track()

  # Push row 5 into the restricted zone, but remove its measurement error
  data$x[5] <- 5
  data$smaj[5] <- NA
  data$smin[5] <- NA
  data$eor[5] <- NA

  expect_error(
    suppressWarnings(suppressMessages(capture.output(checkBarrier(data, spatialCovs = covs, barrier = "coast_barrier")))),
    "Observation at row 5 is in restricted zone but lacks measurement error"
  )
})

test_that("checkBarrier works identically with x.err and y.err (GPS) models", {
  covs <- get_mock_covs()
  data <- get_mock_track()

  # Swap Argos error format for GPS error format
  data$x.err <- data$smaj
  data$y.err <- data$smin
  data$smaj <- NA
  data$smin <- NA
  data$eor <- NA

  # Push row 8 deep into the restricted zone with tiny GPS error
  data$x[8] <- 8
  data$x.err[8] <- 0.01
  data$y.err[8] <- 0.01

  captured <- capture.output({
    res <- suppressWarnings(suppressMessages(checkBarrier(data, spatialCovs = covs, barrier = "coast_barrier")))
  })

  # It should flag the point and set the specific GPS columns to NA
  expect_equal(res$implausible_locations$row_index, 8)
  expect_equal(res$implausible_locations$error_type, "LS")
  expect_true(is.na(res$filtered_data$x[8]))
  expect_true(is.na(res$filtered_data$x.err[8]))
  expect_true(is.na(res$filtered_data$y.err[8]))
})

