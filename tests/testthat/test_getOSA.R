# tests/testthat/test_getOSA.R

# --- Helper to generate a fast, reusable mock fit ---
get_mock_osa_fit <- function() {
  set.seed(42)
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 100, 0, 100))
  terra::values(r) <- runif(100)
  names(r) <- "habitat"
  spatialCovs <- list(habitat = r)

  # 2 animals, 15 obs each to test multi-track looping
  init_par <- list(beta = 0, sigma = 5)
  dat <- suppressMessages(simLangevin(
    obsPerAnimal = 15,
    nbAnimals = 2,
    model = "overdamped",
    spatialCovs = spatialCovs,
    par = init_par,
    measurementError = list(x.sd = 2, y.sd = 2)
  ))

  # Ensure valid dates
  start_time <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dat$date <- start_time + (1:nrow(dat)) * 3600

  # Fit the model with calcOSA = TRUE so we have a baseline to compare against
  fit <- suppressMessages(suppressWarnings({
    fitLangevin(
      data = dat,
      model = "overdamped",
      spatialCovs = spatialCovs,
      par = init_par,
      calcOSA = TRUE,
      silent = TRUE
    )
  }))

  return(list(fit = fit, dat = dat, covs = spatialCovs))
}

# Generate the mock objects once for the test file
mock_env <- get_mock_osa_fit()
fit_base <- mock_env$fit
dat_base <- mock_env$dat
covs_base <- mock_env$covs

test_that("getOSA catches invalid inputs and legacy fit objects", {

  # 1. Non-fitLangevin object
  bad_fit <- list(par = c(1,2,3))
  expect_error(
    getOSA(bad_fit, dat_base, covs_base),
    "'fit' must be a fitLangevin object."
  )

  # 2. Legacy fit object (missing the tmb_setup blueprint)
  legacy_fit <- fit_base
  legacy_fit$tmb_setup <- NULL
  expect_error(
    getOSA(legacy_fit, dat_base, covs_base),
    "does not contain a 'tmb_setup' blueprint"
  )
})

test_that("getOSA returns correctly formatted output", {

  # Run a fresh post-hoc extraction
  osa_fresh <- suppressMessages(getOSA(fit_base, dat_base, covs_base))

  # Check class
  expect_s3_class(osa_fresh, "osaLangevin")
  expect_s3_class(osa_fresh, "data.frame")

  # Check required columns
  expect_true(all(c("id", "date", "residual.x", "residual.y") %in% names(osa_fresh)))

  # Check rows perfectly match the input data
  expect_equal(nrow(osa_fresh), nrow(dat_base))
})

test_that("getOSA post-hoc reconstruction perfectly matches fitLangevin internal calc", {

  # fit_base$osa was calculated INSIDE fitLangevin using the active TMB obj.
  # osa_fresh is calculated by RECONSTRUCTING the TMB obj from the blueprint.
  osa_fresh <- suppressMessages(getOSA(fit_base, dat_base, covs_base))

  # They should be mathematically identical
  expect_equal(osa_fresh$residual.x, fit_base$osa$residual.x)
  expect_equal(osa_fresh$residual.y, fit_base$osa$residual.y)
})

test_that("getOSA successfully accepts and passes alternative OSA methods", {

  # Reconstruct using "fullGaussian" instead of the default
  osa_full <- suppressMessages(
    getOSA(fit_base, dat_base, covs_base, method = "fullGaussian")
  )

  expect_s3_class(osa_full, "osaLangevin")

  # The method should run without error and produce valid numeric residuals.
  # Note: For this specific Gaussian model, fullGaussian and oneStepGaussianOffMode
  # are mathematically identical, so we just verify it produced valid numerics.
  valid_residuals <- na.omit(osa_full$residual.x)
  expect_true(length(valid_residuals) > 0)
  expect_true(is.numeric(valid_residuals))
})

test_that("getOSA handles individual track failures gracefully without crashing the whole process", {

  bad_dat <- dat_base

  # Sabotage Track 2 by leaving only its FIRST observation intact.
  # This triggers the "not enough valid observations" warning, but anchors
  # the random walk so TMB's precision matrix doesn't become singular and crash Track 1.
  t2_idx <- which(bad_dat$id == "2")
  bad_dat$x[t2_idx[-1]] <- NA
  bad_dat$y[t2_idx[-1]] <- NA

  # need to update the fit object's signature so it allows this
  bad_fit <- fit_base
  bad_fit$signatures$data <- langevinSSM:::get_data_signature(bad_dat)

  # Expect a warning that Track 2 didn't have enough observations
  expect_warning(
    osa_sabotaged <- suppressMessages(getOSA(bad_fit, bad_dat, covs_base)),
    "Track ID 2 does not have enough valid observations"
  )

  # BUT the function should not have crashed. It should still return an object...
  expect_s3_class(osa_sabotaged, "osaLangevin")

  # ...where Animal 1 successfully got its residuals...
  anim_1_resids <- na.omit(osa_sabotaged$residual.x[osa_sabotaged$id == "1"])
  expect_true(length(anim_1_resids) > 0)

  # ...and Animal 2 is entirely NA
  anim_2_resids <- na.omit(osa_sabotaged$residual.x[osa_sabotaged$id == "2"])
  expect_equal(length(anim_2_resids), 0)
})
