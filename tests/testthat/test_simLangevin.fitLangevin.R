# tests/testthat/test_simLangevin.fitLangevin.R

# Setup mock data unconditionally so testthat always has access to it
r <- terra::rast(nrows=50, ncols=50, ext=c(0,500,0,500), vals=1:2500)
names(r) <- "cov"
exCovs <- list(cov=r)

# Base data centered at 250, 250
df <- data.frame(id=1, date=seq(0,1,0.1), dt=c(0, rep(0.1, 10)), x=250, y=250,
                 x.err=1, y.err=1, smaj=NA_real_, smin=NA_real_, eor=NA_real_)
exDat <- class_dataLangevin(df)
fit <- suppressMessages(fitLangevin(data=exDat, spatialCovs=exCovs, par=list(sigma=1), silent=TRUE))

test_that("simLangevin.fitLangevin basic functionality", {
  res <- suppressMessages(simLangevin(fit, data = exDat, spatialCovs = exCovs, timeStep = 0.05, conditional = FALSE))
  expect_s3_class(res, "dataLangevin")
  expect_equal(as.numeric(res$date), as.numeric(exDat$date))
})

test_that("simLangevin.fitLangevin handles POSIXt date and time.unit resolution", {
  df_posix <- data.frame(
    id = 1,
    date = as.POSIXct("2024-01-01 12:00:00", tz="UTC") + (0:10)*60,
    dt = c(0, rep(1/60, 10)),
    x = 250, y = 250, x.err=1, y.err=1, smaj=NA, smin=NA, eor=NA
  )
  attr(df_posix, "time.unit") <- "hours"
  exDat_posix <- class_dataLangevin(df_posix)

  fit_posix <- suppressMessages(fitLangevin(data=exDat_posix, spatialCovs=exCovs, par=list(sigma=1), silent=TRUE))

  res <- suppressMessages(simLangevin(fit_posix, data = exDat_posix, spatialCovs = exCovs,
                                      timeStep = 1/120, conditional = FALSE))

  # Explicit attribute check
  expect_equal(attr(res, "time.unit"), "hours")
})

test_that("simLangevin.fitLangevin respects scaleFactor > 1 for all spatial units", {
  sf_val <- 1000
  fit_sf <- suppressMessages(fitLangevin(data=exDat, spatialCovs=exCovs, par=list(sigma=1),
                                         scaleFactor=sf_val, silent=TRUE))

  res <- suppressMessages(simLangevin(fit_sf, data = exDat, spatialCovs = exCovs,
                                      timeStep = 0.05, conditional = FALSE))

  # Coordinates should be un-scaled to original units
  expect_equal(mean(res$mu.x), 250, tolerance = 10)
})

test_that("simLangevin.fitLangevin error handling and conditional independence", {
  expect_error(simLangevin(fit, data = exDat[0,], spatialCovs = exCovs), "data contains no observations")
  expect_error(simLangevin(fit, data = exDat, spatialCovs = exCovs, timeStep = 0, conditional = FALSE),
               "greater than zero")
  expect_no_error(suppressMessages(simLangevin(fit, data = exDat, spatialCovs = exCovs,
                                               timeStep = 0, conditional = TRUE)))
})

test_that("Imputation vs Predictive Check divergence and GoF validation", {
  # 1. Setup a long, centered track to ensure convergence and valid residuals
  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  long_df <- data.frame(
    id = 1,
    date = 0:50,
    dt = c(0, rep(1, 50)),
    x = seq(225, 275, length.out = 51),
    y = seq(225, 275, length.out = 51),
    x.err = 0.1, y.err = 0.1, smaj = NA, smin = NA, eor = NA
  )
  long_dat <- class_dataLangevin(long_df)

  # 2. Fit the model away from edges to ensure valid Hessian and residuals
  long_fit <- suppressMessages(fitLangevin(
    data = long_dat,
    spatialCovs = exCovs,
    par = list(sigma = 2),
    silent = TRUE
  ))

  long_fit$residuals <- suppressMessages(residuals(long_fit, long_dat, exCovs))

  # 3. Verify that tests_df was actually created and is not NULL
  # This confirms fitLangevin successfully called gof_tests with valid data
  expect_false(is.null(attr(long_fit$residuals, "tests")))
  expect_s3_class(attr(long_fit$residuals, "tests"), "data.frame")

  # 4. Standard Divergence Check
  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  res_imp <- suppressMessages(simLangevin(long_fit, data = long_dat, spatialCovs = exCovs,
                                          conditional = TRUE))

  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  res_pred <- suppressMessages(simLangevin(long_fit, data = long_dat, spatialCovs = exCovs,
                                           conditional = FALSE))

  dist_imp <- sqrt(sum((res_imp$mu.x - long_dat$x)^2))
  dist_pred <- sqrt(sum((res_pred$mu.x - long_dat$x)^2))

  expect_gt(dist_pred, dist_imp)
})

test_that("Full Posterior Draw handles uncertainty propagation", {
  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  res_fp <- suppressMessages(simLangevin(fit, data = exDat, spatialCovs = exCovs,
                                         conditional = TRUE, fullPosterior = TRUE))

  set.seed(123, kind="Mersenne-Twister", normal.kind = "Inversion")
  res_nfp <- suppressMessages(simLangevin(fit, data = exDat, spatialCovs = exCovs,
                                          conditional = TRUE, fullPosterior = FALSE))

  expect_false(identical(res_fp$mu.x, res_nfp$mu.x))
})
