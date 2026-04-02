# tests/testthat/test_fitLangevin.R

# --- Helpers ---
get_valid_raster <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 100, 0, 100))
  terra::values(r) <- 1:100
  names(r) <- "habitat"
  return(r)
}

get_valid_dataLangevin <- function() {
  df <- data.frame(
    id = as.factor(rep("A", 5)), # FIXED: Must be a factor to mimic formatData
    date = as.POSIXct(c("2023-01-01 10:00:00", "2023-01-01 11:00:00", "2023-01-01 12:00:00",
                        "2023-01-01 13:00:00", "2023-01-01 14:00:00"), tz = "UTC"),
    dt = rep(1, 5),
    x = c(40, 45, 50, 55, 60),
    y = c(40, 45, 50, 55, 60),
    x.sd = rep(1, 5),
    y.sd = rep(1, 5),
    smaj = rep(NA, 5),
    smin = rep(NA, 5),
    eor = rep(NA, 5),
    lc = as.factor(rep("G", 5)) # FIXED: Also made factor for consistency
  )
  attr(df, "time.unit") <- "hours"
  class(df) <- append("dataLangevin", class(df))
  return(df)
}

get_valid_par <- function() {
  list(beta = c(0.5), sigma = 5, gamma = 0.5)
}

test_that("fitLangevin rejects unformatted data", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  # Just a regular data frame, not run through formatData()
  bad_data <- data.frame(id = "A", x = 50, y = 50)

  expect_error(fitLangevin(data = bad_data, spatialCovs = r, par = p),
               "is not formatted as a 'dataLangevin' object")
})

test_that("fitLangevin catches missing coordinate columns", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()
  dat <- get_valid_dataLangevin()

  # Tell fitLangevin to look for lon/lat instead of x/y
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p, coord = c("lon", "lat")),
               "coord not found in data")
})

test_that("fitLangevin enforces smoothGradient rules", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()
  dat <- get_valid_dataLangevin()

  # Invalid npoints
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p,
                           smoothGradient = TRUE, npoints = 5),
               "npoints must be 4 or 8")

  # Invalid curweight (must be >= 0 and < 1)
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p,
                           smoothGradient = TRUE, curweight = 1.5),
               "curweight must be >=0 and <1")
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p,
                           smoothGradient = TRUE, curweight = -0.1),
               "curweight must be >=0 and <1")

  # Invalid zetaScale
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p,
                           smoothGradient = TRUE, zetaScale = -2),
               "zetaScale must be >0")
})

test_that("fitLangevin catches spatial overlap errors from prepareRaster", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()
  dat <- get_valid_dataLangevin()

  # Push coordinates completely outside the 0-100 raster extent
  dat$x <- c(200, 205, 210, 215, 220)

  # This tests that fitLangevin correctly passes the data to prepareRaster
  # and fast-fails BEFORE attempting any slow TMB evaluations
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p),
               "overlap with 'spatialCovs'")
})

test_that("fitLangevin validates parameter lengths and map structures", {
  r <- list(habitat = get_valid_raster()) # 1 covariate
  p <- get_valid_par()
  dat <- get_valid_dataLangevin()

  # beta length mismatch (provided 2 coefficients for 1 covariate)
  p_bad_beta <- p
  p_bad_beta$beta <- c(0.5, -0.2)
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p_bad_beta),
               "beta")

  # tau length mismatch (provided a single scalar instead of a 2-vector)
  p_bad_tau <- p
  p_bad_tau$tau <- 1.5
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p_bad_tau),
               "tau")

  # rho_o bounds check (correlation must theoretically be bounded, though TMB maps it)
  p_bad_rho <- p
  p_bad_rho$rho_o <- 1.5
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p_bad_rho),
               "rho_o")

  nonsense_p <- list(beta=0,gamma=0.5,sigma=1, nonsense=123)
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = nonsense_p),
               "par")

  # map mismatch
  # If a user tries to map out tau (which is intrinsically length 2) but only provides a length-1 factor
  bad_map <- list(tau = as.factor(1))
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p, map = bad_map),
               "map")

  # If a user tries to map out nonsense
  nonsense_map <- list(nonsense = as.factor(1))
  expect_error(fitLangevin(data = dat, spatialCovs = r, par = p, map = nonsense_map),
               "map")
})

test_that("fitLangevin rejects observation parameters when data lacks corresponding errors", {
  r <- list(habitat = get_valid_raster())
  dat <- get_valid_dataLangevin() # This mock has x.sd/y.sd, but NO smaj/smin/eor
  p <- get_valid_par()

  # attempt to pass psi (should fail because there is no ellipse data)
  p_bad_psi <- p
  p_bad_psi$psi <- 1.5
  expect_error(suppressMessages(fitLangevin(data = dat, spatialCovs = r, par = p_bad_psi)),
               "error ellipse observations")

  # map psi should also fail
  expect_error(suppressMessages(fitLangevin(data = dat, spatialCovs = r, par = p, map = list(psi = as.factor(NA)))),
               "error ellipse observations")

  # now strip LS/GPS data to test tau/rho rejection
  dat_no_err <- dat
  dat_no_err$x.sd <- NA
  dat_no_err$y.sd <- NA

  p_bad_tau <- p
  p_bad_tau$tau <- c(1, 1)
  expect_error(suppressMessages(fitLangevin(data = dat_no_err, spatialCovs = r, par = p_bad_tau)),
               "standard deviation observations")

  # map rho_o should fail
  expect_error(suppressMessages(fitLangevin(data = dat_no_err, spatialCovs = r, par = p, map = list(rho_o = as.factor(NA)))),
               "standard deviation observations")
})
