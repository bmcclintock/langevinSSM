# tests/testthat/test_plot.fitLangevin.R

# --- Helper functions for mock data ---
get_mock_covs <- function(dynamic = FALSE) {
  r1 <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 10, 0, 10), vals = runif(100))
  if (dynamic) {
    r1 <- c(r1, r1 * 0.5)
    terra::time(r1) <- as.POSIXct(c("2023-01-01", "2023-01-02"), tz = "UTC")
  }
  list(cov1 = r1)
}

get_mock_fit <- function() {
  fit <- list(
    estimates = list(
      natural = matrix(c(1.5), nrow = 1, dimnames = list(c("beta"), "Estimate")),
      random = list(
        mu = list(
          est = data.frame(mu.x = c(4, 6), mu.y = c(4, 6), id = c("A", "B"))
        )
      )
    )
  )
  class(fit) <- "fitLangevin"
  return(fit)
}

get_mock_data <- function() {
  dat <- data.frame(x = c(3.9, 6.1), y = c(4.1, 5.9), id = c("A", "B"))
  class(dat) <- c("dataLangevin", "data.frame")
  return(dat)
}

test_that("plot.fitLangevin catches missing spatialCovs", {
  fit <- get_mock_fit()
  expect_error(plot(fit), "You must provide the 'spatialCovs'")
})

test_that("plot.fitLangevin works for static rasters with auto-extent", {
  fit <- get_mock_fit()
  covs <- get_mock_covs(dynamic = FALSE)

  # Auto-extent calculation triggers when extent = NULL
  p <- plot(fit, spatialCovs = covs)
  expect_s3_class(p, "ggplot")
})

test_that("plot.fitLangevin works with observed data overlay", {
  fit <- get_mock_fit()
  covs <- get_mock_covs(dynamic = FALSE)
  dat <- get_mock_data()

  p <- plot(fit, spatialCovs = covs, data = dat)
  expect_s3_class(p, "ggplot")
})

test_that("plot.fitLangevin handles manual extent arguments", {
  fit <- get_mock_fit()
  covs <- get_mock_covs(dynamic = FALSE)

  # Using numeric vector
  p1 <- plot(fit, spatialCovs = covs, extent = c(2, 8, 2, 8))
  expect_s3_class(p1, "ggplot")

  # Using SpatExtent
  p2 <- plot(fit, spatialCovs = covs, extent = terra::ext(2, 8, 2, 8))
  expect_s3_class(p2, "ggplot")
})

test_that("plot.fitLangevin handles dynamic rasters with and without time subsetting", {
  fit <- get_mock_fit()
  covs <- get_mock_covs(dynamic = TRUE)

  # All layers (faceting should trigger)
  p_all <- plot(fit, spatialCovs = covs)
  expect_s3_class(p_all, "ggplot")

  # Subsetting by time index
  p_sub <- plot(fit, spatialCovs = covs, time = 1)
  expect_s3_class(p_sub, "ggplot")
})

test_that("plot.fitLangevin compact = FALSE returns a list of individual plots", {
  fit <- get_mock_fit()
  covs <- get_mock_covs(dynamic = FALSE)

  p_list <- plot(fit, spatialCovs = covs, compact = FALSE)

  # Because the mock data has 2 IDs ("A" and "B"), it should return a list of 2 ggplots
  expect_true(is.list(p_list))
  expect_equal(length(p_list), 2)
  expect_s3_class(p_list[[1]], "ggplot")
  expect_s3_class(p_list[[2]], "ggplot")
})
