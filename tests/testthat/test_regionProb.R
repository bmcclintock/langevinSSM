# tests/testthat/test_regionProb.R

# --- Helpers ---
get_mock_covs_and_mask <- function(dynamic = FALSE) {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 10, 0, 10))
  terra::values(r) <- runif(100)
  names(r) <- "habitat"

  mask <- r
  terra::values(mask) <- ifelse(1:100 > 50, 1, 0)

  if (dynamic) {
    r2 <- r * 0.5
    r_dyn <- c(r, r2)
    terra::time(r_dyn) <- as.POSIXct(c("2023-01-01", "2023-01-02"), tz = "UTC")
    return(list(covs = list(habitat = r_dyn), mask = mask))
  }

  return(list(covs = list(habitat = r), mask = mask))
}

get_mock_fit_for_rp <- function() {
  fit <- list(
    estimates = list(
      natural = data.frame(Estimate = c(0.5, 2.0), row.names = c("beta_habitat", "sigma"))
    ),
    covariance = list(
      natural = matrix(c(0.1, 0, 0, 0.05), nrow = 2, dimnames = list(c("beta_habitat", "sigma"), c("beta_habitat", "sigma")))
    ),
    signatures = list(data = NULL, covs = NULL)
  )
  class(fit) <- "fitLangevin"
  return(fit)
}

test_that("regionProb catches invalid inputs", {
  env <- get_mock_covs_and_mask()
  fit <- get_mock_fit_for_rp()

  expect_error(regionProb(fit = list(), spatialCovs = env$covs, mask = env$mask), "must be a fitLangevin object")
  expect_error(regionProb(fit = fit, spatialCovs = env$mask, mask = env$mask), "must be a list of SpatRaster")
  expect_error(regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = -10), "single non-negative integer")
})

test_that("regionProb calculates dynamic multi-layer bounds correctly", {
  env <- get_mock_covs_and_mask(dynamic = TRUE)
  fit <- get_mock_fit_for_rp()

  out <- suppressMessages(regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = 10, show_progress = FALSE))

  expect_s3_class(out, "regLangevin")
  expect_length(out$Point_Estimate, 2)
  expect_equal(nrow(out$CI_delta), 2)
  expect_equal(terra::nlyr(out$prob_raster), 2)
})

test_that("regionProb respects level argument", {
  env <- get_mock_covs_and_mask()
  fit <- get_mock_fit_for_rp()

  out_95 <- suppressMessages(regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = 10, level = 0.95, show_progress = FALSE))
  out_50 <- suppressMessages(regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = 10, level = 0.50, show_progress = FALSE))

  width_95 <- out_95$CI_delta[1, 2] - out_95$CI_delta[1, 1]
  width_50 <- out_50$CI_delta[1, 2] - out_50$CI_delta[1, 1]

  expect_true(width_50 < width_95)
})

test_that("plot.regLangevin handles log argument and dynamic plots safely", {
  env <- get_mock_covs_and_mask(dynamic = TRUE)
  fit <- get_mock_fit_for_rp()

  out <- suppressMessages(regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = 0, show_progress = FALSE))

  p_raw <- plot(out, log = FALSE)
  p_log <- plot(out, log = TRUE)

  expect_s3_class(p_raw, "ggplot")
  expect_s3_class(p_log, "ggplot")
  expect_true(inherits(p_raw$facet, "FacetWrap"), info = "Multi-layer plotting should facet")
})

test_that("print.regLangevin formats multi-layer output correctly", {
  env <- get_mock_covs_and_mask(dynamic = TRUE)
  fit <- get_mock_fit_for_rp()

  out_sim <- suppressMessages(regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = 10, level = 0.95, show_progress = FALSE))

  printed_sim <- utils::capture.output(print(out_sim))

  # Check that it splits the output into the 3 distinct blocks
  expect_true(any(grepl("--- Point Estimates ---", printed_sim)))
  expect_true(any(grepl("--- Delta Method Approximation ---", printed_sim)))
  expect_true(any(grepl("--- Monte Carlo Simulation ---", printed_sim)))

  # Check that the dynamic CI headers rendered properly
  expect_true(any(grepl("95%_CI", printed_sim)))
})

test_that("plot.regLangevin handles log argument and dynamic plots safely", {
  env <- get_mock_covs_and_mask(dynamic = TRUE)
  fit <- get_mock_fit_for_rp()

  out <- suppressMessages(regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = 0, show_progress = FALSE))

  p_raw <- plot(out, log = FALSE)
  p_log <- plot(out, log = TRUE)

  expect_s3_class(p_raw, "ggplot")
  expect_s3_class(p_log, "ggplot")
  expect_true(inherits(p_raw$facet, "FacetWrap"), info = "Multi-layer plotting should facet")

  # Safely extract the labeller mapping directly from the facet environment
  label_mapping <- environment(p_raw$facet$params$labeller)$mapping

  # Ensure the mapping successfully captured our 'Prob:' strings
  expect_true(all(grepl("Prob: ", label_mapping)))
})

test_that("plot.regLangevin warns if custom extent crops the active region", {
  env <- get_mock_covs_and_mask(dynamic = FALSE)
  fit <- get_mock_fit_for_rp()

  out <- suppressMessages(regionProb(fit = fit, spatialCovs = env$covs, mask = env$mask, nSims = 0, show_progress = FALSE))

  # aggressively small extent guaranteed to crop out some of the active mock mask
  tiny_extent <- c(0, 2, 0, 2)

  expect_warning(
    plot(out, extent = tiny_extent),
    "crops out parts of the region of interest"
  )

  # massive extent fully encapsulating the active mask should NOT trigger the warning
  massive_extent <- c(-10, 20, -10, 20)
  expect_silent(plot(out, extent = massive_extent))
})
