# tests/testthat/test_subSampleData.R

# --- Helper to create mock dataLangevin objects ---
get_mock_dataLangevin <- function(n_per_id = 100) {

  # Create hourly time steps
  dates_A <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC") + (0:(n_per_id - 1)) * 3600
  dates_B <- as.POSIXct("2024-01-02 00:00:00", tz = "UTC") + (0:(n_per_id - 1)) * 3600

  df <- data.frame(
    id = rep(c("A", "B"), each = n_per_id),
    date = c(dates_A, dates_B),
    x = rnorm(n_per_id * 2),
    y = rnorm(n_per_id * 2),
    smaj = 1, smin = 0.5, eor = 0,
    x.sd = 1, y.sd = 1,
    stringsAsFactors = FALSE
  )

  # Initialize standard dt vector (0 for first obs, 1 hour otherwise)
  df$dt <- rep(1, nrow(df))
  df$dt[c(1, n_per_id + 1)] <- 0

  # Assign proper class and attributes
  class(df) <- c("dataLangevin", "data.frame")
  attr(df, "time.unit") <- "hours"

  return(df)
}

test_that("subSampleData catches invalid inputs", {
  mock_data <- get_mock_dataLangevin(10)

  # Non-dataLangevin class
  bad_data <- unclass(mock_data)
  expect_error(subSampleData(bad_data), "formatted as a 'dataLangevin' object")

  # Invalid samplingRate
  expect_error(subSampleData(mock_data, samplingRate = 0), "samplingRate should be an integer >= 1")
  expect_error(subSampleData(mock_data, samplingRate = "5"), "samplingRate should be an integer >= 1")

  # Invalid propMissing
  expect_error(subSampleData(mock_data, propMissing = -0.1), "numeric value >=0 and <1")
  expect_error(subSampleData(mock_data, propMissing = 1), "numeric value >=0 and <1")

  # Missing columns for NA introduction trigger a warning
  expect_warning(
    subSampleData(mock_data, propMissing = 0.5, col_to_na = c("fake_col")),
    "No valid columns specified in 'col_to_na'"
  )
})

test_that("subSampleData returns correct data dimensions and attributes", {
  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  mock_data <- get_mock_dataLangevin(100) # 200 rows total

  # Subsampling rate = 2 (should return ceiling(200 / 2) = 100 rows)
  sub_rate2 <- subSampleData(mock_data, samplingRate = 2)
  expect_equal(nrow(sub_rate2), 100)
  expect_s3_class(sub_rate2, "dataLangevin")
  expect_equal(attr(sub_rate2, "time.unit"), "hours")

  # Subsampling rate = 10 (should return ceiling(200 / 10) = 20 rows)
  sub_rate10 <- subSampleData(mock_data, samplingRate = 10)
  expect_equal(nrow(sub_rate10), 20)
})

test_that("subSampleData always preserves the first observation of every track", {
  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  mock_data <- get_mock_dataLangevin(50)

  sub_data <- subSampleData(mock_data, samplingRate = 5)

  # Extract the very first timestamps of each track in the original dataset
  first_dates_orig <- tapply(mock_data$date, mock_data$id, min)

  # Extract the very first timestamps of each track in the subsampled dataset
  first_dates_sub <- tapply(sub_data$date, sub_data$id, min)

  expect_equal(as.numeric(first_dates_orig), as.numeric(first_dates_sub))
})

test_that("subSampleData properly recalculates the dt (time step) vector", {
  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  mock_data <- get_mock_dataLangevin(50)

  # Subsample
  sub_data <- subSampleData(mock_data, samplingRate = 5)

  # The first observation of each track must have a dt of 0
  first_obs_idx <- !duplicated(sub_data$id)
  expect_true(all(sub_data$dt[first_obs_idx] == 0))

  # For subsequent observations, dt should accurately reflect date differences in "hours"
  track_A <- sub_data[sub_data$id == "A", ]
  calculated_dt_A <- as.numeric(difftime(track_A$date[-1], track_A$date[-nrow(track_A)], units = "hours"))
  expect_equal(track_A$dt[-1], calculated_dt_A)
})

test_that("subSampleData randomly introduces missing values to specified columns", {
  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  mock_data <- get_mock_dataLangevin(100)

  # Prop missing = 0.5, rate = 1
  sub_data <- subSampleData(mock_data, samplingRate = 1, propMissing = 0.5, col_to_na = c("x", "y"))

  # Total missing should be roughly 50% of the non-first observations
  pool_size <- sum(duplicated(mock_data$id)) # 198 rows eligible for NA
  na_count_x <- sum(is.na(sub_data$x))

  # Verify NAs were successfully injected (binomial draw means it won't be exactly 99, but > 0)
  expect_gt(na_count_x, 0)

  # Verify NAs were only injected into the specified columns
  expect_true(any(is.na(sub_data$y)))
  expect_false(any(is.na(sub_data$smaj)))
  expect_false(any(is.na(sub_data$x.sd)))
})

test_that("subSampleData NEVER introduces NAs to the first observation of a track", {
  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  mock_data <- get_mock_dataLangevin(50)

  # Force a massive missing proportion (99%) to heavily test the safeguard
  sub_data <- subSampleData(mock_data, propMissing = 0.99, col_to_na = c("x", "y"))

  first_obs_idx <- !duplicated(sub_data$id)

  # Ensure the first observations remain perfectly intact
  expect_false(any(is.na(sub_data$x[first_obs_idx])))
  expect_false(any(is.na(sub_data$y[first_obs_idx])))
})
