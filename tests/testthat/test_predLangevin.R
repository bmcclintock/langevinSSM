# tests/testthat/test_predLangevin.R

# Helper to create lightweight mock objects to bypass actual TMB fitting
get_predict_mocks <- function() {
  # Mock dataLangevin object (Pre-padded by formatData)
  mock_data <- data.frame(
    id = factor(c("SealA_1", "SealA_1")),
    date = as.POSIXct(c("2023-01-01 10:00:00", "2023-01-01 11:00:00"), tz = "UTC"),
    dt = c(0, 1),
    x = c(500000, NA), # NA represents the padded prediction row
    y = c(4000000, NA)
  )
  class(mock_data) <- c("dataLangevin", "data.frame")
  attr(mock_data, "time.unit") <- "hours"

  # Mock spatialCovs
  mock_covs <- list(dummy = terra::rast(matrix(0, 10, 10)))
  terra::ext(mock_covs$dummy) <- c(0, 10, 0, 10) # Give it a defined extent

  # Generate exact, valid signatures using the package's internal functions
  ns <- asNamespace("langevinSSM")
  unpadded_mock <- mock_data[!is.na(mock_data$x), ]
  data_sig <- ns$get_data_signature(unpadded_mock, coord = c("x", "y"))
  covs_sig <- ns$get_covs_signature(mock_covs)

  # Mock fitLangevin object
  # (Notice psi is NOT in par, mimicking a default fit)
  mock_fit <- list(
    par = c(beta = 0.5, log_sigma = 1.2, log_gamma = -0.5),
    conditions = list(scaleFactor = 1, model = "underdamped", barrier = NULL, lambda = NULL, coord = c("x", "y")),
    tmb_setup = list(
      parList = list(beta = c(0.5), log_sigma = 1.2, log_gamma = -0.5)
    ),
    signatures = list(data = data_sig, covs = covs_sig)
  )
  class(mock_fit) <- "fitLangevin"

  return(list(data = mock_data, fit = mock_fit, covs = mock_covs))
}

# Wrapper to safely hijack the final fitLangevin call
with_mockery <- function(code) {
  old_fit <- langevinSSM::fitLangevin

  assignInNamespace("fitLangevin", function(...) {
    res <- list(
      success = TRUE,
      args = list(...),
      tmb_setup = list(parList = list(mu = matrix(1, 2, 2)))
    )
    class(res) <- "fitLangevin"
    return(res)
  }, ns = "langevinSSM")

  # Mock plot.fitLangevin so it returns NULL instead of opening a graphics device
  old_plot <- tryCatch(get("plot.fitLangevin", envir = asNamespace("langevinSSM")), error = function(e) NULL)
  if (!is.null(old_plot)) {
    assignInNamespace("plot.fitLangevin", function(...) invisible(NULL), ns = "langevinSSM")
  }

  on.exit({
    assignInNamespace("fitLangevin", old_fit, ns = "langevinSSM")
    if (!is.null(old_plot)) assignInNamespace("plot.fitLangevin", old_plot, ns = "langevinSSM")
  }, add = TRUE)

  force(code)
}

test_that("Base argument validation works correctly", {
  mocks <- get_predict_mocks()

  with_mockery({
    # We use capture.output() to swallow any print(NULL) calls from the mocked plots
    capture.output({
      expect_error(suppressMessages(predLangevin(structure(list(), class="not_a_fit"), data = mocks$data, spatialCovs = mocks$covs)),
                   "must be of class 'fitLangevin'")

      expect_error(suppressMessages(predLangevin(mocks$fit, data = data.frame(x=1, y=2), spatialCovs = mocks$covs)),
                   "must be provided and must be a formatted 'dataLangevin' object")

      expect_error(suppressMessages(predLangevin(mocks$fit, data = mocks$data, spatialCovs = "not_a_list")),
                   "must be provided as a list of SpatRaster")
    })
  })
})

test_that("verify_signatures correctly ignores data mismatch but catches modified covariates", {
  mocks <- get_predict_mocks()

  with_mockery({
    # 1. Baseline: Should pass silently
    capture.output({
      res <- suppressWarnings(suppressMessages(
        predLangevin(mocks$fit, data = mocks$data, spatialCovs = mocks$covs)
      ))
    })
    expect_s3_class(res, "predLangevin")

    # 2. Tampered Covariates: Modify the spatial extent
    bad_covs <- mocks$covs
    terra::ext(bad_covs$dummy) <- c(0, 100, 0, 100)

    capture.output({
      expect_error(suppressMessages(predLangevin(mocks$fit, data = mocks$data, spatialCovs = bad_covs)),
                   "Safeguard triggered: the 'spatialCovs' provided do not match")
    })
  })
})

test_that("predLangevin correctly freezes MLEs and calls inner optimization", {
  mocks <- get_predict_mocks()

  with_mockery({
    capture.output({
      res <- suppressMessages(predLangevin(mocks$fit, data = mocks$data, spatialCovs = mocks$covs))
    })

    expect_s3_class(res, "predLangevin")
    expect_s3_class(res, "fitLangevin")

    passed_args <- res$args
    expect_true("map" %in% names(passed_args))
    frozen_map <- passed_args$map

    # Check that estimated parameters are mapped to NA (frozen)
    expect_true(all(is.na(frozen_map$beta)))
    expect_true(is.na(frozen_map$sigma))
    expect_true(is.na(frozen_map$gamma))

    # Because psi wasn't in the original fit$par, it shouldn't exist in the map
    expect_false("psi" %in% names(frozen_map))

    # Verify the padded data was passed through
    expect_equal(nrow(passed_args$data), 2)
  })
})
