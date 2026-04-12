# tests/testthat/test_initialValues.R

# --- Helper functions for testing ---
get_mock_covs <- function() {
  # Mock spatial covariates list of length 2
  list(cov1 = "dummy", cov2 = "dummy")
}

get_mock_dataLangevin <- function() {
  # Carefully constructed to test:
  # 1. Multiple tracks (Track A and Track B)
  # 2. NA gaps in the middle of a track
  # 3. NA trailing at the end of a track (to test rule = 2)
  # 4. dt differences
  df <- data.frame(
    id = as.factor(c("A", "A", "A", "A", "A", "B", "B")),
    date = as.POSIXct(c(
      "2023-01-01 00:00:00", "2023-01-01 01:00:00", "2023-01-01 02:00:00",
      "2023-01-01 03:00:00", "2023-01-01 04:00:00",
      "2023-01-01 00:00:00", "2023-01-01 01:00:00"
    ), tz = "UTC"),
    dt = c(0, 1, 1, 1, 1, 0, 1),
    x = c(0, NA, 4, 6, NA, 10, 10),
    y = c(0, NA, 0, 0, NA, 10, 12)
  )
  attr(df, "time.unit") <- "hours"
  class(df) <- append("dataLangevin", class(df))
  return(df)
}

test_that("initialValues catches invalid input types", {
  dat <- get_mock_dataLangevin()
  covs <- get_mock_covs()

  # data not dataLangevin
  bad_data <- data.frame(id = "A", x = 1, y = 1)
  expect_error(initialValues(data = bad_data, par = list(), spatialCovs = covs),
               "is not formatted as a 'dataLangevin' object")

  # par not a list
  expect_error(initialValues(data = dat, par = c(sigma = 1), spatialCovs = covs),
               "par must be a list")

  # par has invalid names
  expect_error(initialValues(data = dat, par = list(sigma = 1, invalid_param = 5), spatialCovs = covs),
               "is limited to c")
})

test_that("initialValues generates correct defaults for overdamped model", {
  dat <- get_mock_dataLangevin()
  covs <- get_mock_covs()

  par_init <- initialValues(data = dat, model = "overdamped", spatialCovs = covs)

  # Check parameter existence
  expect_true(!is.null(par_init$beta))
  expect_true(!is.null(par_init$sigma))
  expect_true(!is.null(par_init$mu))

  # Check parameters NOT present in overdamped
  expect_true(is.null(par_init$gamma))
  expect_true(is.null(par_init$vel))

  # Check beta length matches covariates
  expect_equal(length(par_init$beta), 2)
  expect_equal(par_init$beta, c(0, 0))

  # Check empirical sigma is positive and numeric
  expect_true(is.numeric(par_init$sigma))
  expect_true(par_init$sigma > 0)
})

test_that("initialValues generates correct defaults for underdamped model", {
  dat <- get_mock_dataLangevin()
  covs <- get_mock_covs()

  par_init <- initialValues(data = dat, model = "underdamped", spatialCovs = covs)

  # Check parameter existence unique to underdamped
  expect_true(!is.null(par_init$gamma))
  expect_true(!is.null(par_init$vel))

  # Check gamma is positive and numeric
  expect_true(is.numeric(par_init$gamma))
  expect_true(par_init$gamma > 0)

  # Check vel matrix dimensions
  expect_equal(nrow(par_init$vel), nrow(dat))
  expect_equal(ncol(par_init$vel), 2)
  expect_true(all(par_init$vel == 0))
})

test_that("User-provided par values override empirical estimates", {
  dat <- get_mock_dataLangevin()

  # inject dummy standard deviations so it is legal to test rho_o
  dat$x.err <- rep(1, nrow(dat))
  dat$y.err <- rep(1, nrow(dat))

  covs <- get_mock_covs()

  user_par <- list(
    sigma = 99.9,
    gamma = 50.5,
    beta = c(1.1, 2.2),
    rho_o = 0.5
  )

  par_init <- initialValues(data = dat, model = "underdamped", par = user_par, spatialCovs = covs)

  # Ensure overrides worked
  expect_equal(par_init$sigma, 99.9)
  expect_equal(par_init$gamma, 50.5)
  expect_equal(par_init$beta, c(1.1, 2.2))
  expect_equal(par_init$rho_o, 0.5)
})

test_that("Missing coordinates in mu are linearly interpolated correctly", {
  dat <- get_mock_dataLangevin()
  covs <- get_mock_covs()

  par_init <- initialValues(data = dat, model = "underdamped", spatialCovs = covs)
  mu <- par_init$mu

  # Check no NAs remain
  expect_false(any(is.na(mu)))

  # Track A: Gap in the middle (t=1). Should interpolate exactly halfway between t=0 (x=0) and t=2 (x=4)
  expect_equal(mu[2, 1], 2)
  expect_equal(mu[2, 2], 0)

  # Track A: NA at the end (t=4). rule=2 should carry the last known observation forward (t=3, x=6)
  expect_equal(mu[5, 1], 6)
  expect_equal(mu[5, 2], 0)

  # Track B should be completely unaffected
  expect_equal(mu[6, 1], 10)
  expect_equal(mu[7, 2], 12)
})

test_that("Empirical sigma and gamma accurately calculate across NA gaps", {
  dat <- get_mock_dataLangevin()
  covs <- get_mock_covs()

  par_init <- initialValues(data = dat, model = "underdamped", spatialCovs = covs)

  # Mathematical manual check of the valid jumps in get_mock_dataLangevin():
  # Jump 1 (Track A): t=0 to t=2 (dt=2, dx=4, dy=0) -> R^2 = 16. Sigma^2_1 = 16 / (2*2) = 4
  # Jump 2 (Track A): t=2 to t=3 (dt=1, dx=2, dy=0) -> R^2 = 4. Sigma^2_2 = 4 / (2*1) = 2
  # Jump 3 (Track B): t=0 to t=1 (dt=1, dx=0, dy=2) -> R^2 = 4. Sigma^2_3 = 4 / (2*1) = 2

  # Mean Sigma^2 = (4 + 2 + 2) / 3 = 8 / 3
  expected_sigma <- sqrt(8 / 3)

  # dt's evaluated: c(2, 1, 1). Median dt = 1.
  # Expected gamma = 1 / median(dt) = 1 / 1 = 1

  expect_equal(par_init$sigma, expected_sigma)
  expect_equal(par_init$gamma, 1)
})

test_that("initialValues validates user-provided par$mu", {
  dat <- get_mock_dataLangevin()
  covs <- get_mock_covs()

  # wrong Dimensions (Too few rows)
  bad_mu_dim <- matrix(0, nrow = nrow(dat) - 1, ncol = 2)
  expect_error(initialValues(data = dat, par = list(mu = bad_mu_dim), spatialCovs = covs),
               "must be a matrix with the same number of rows as 'data' and 2 columns")

  # wrong Format (Vector instead of Matrix)
  bad_mu_vec <- rep(0, nrow(dat) * 2)
  expect_error(initialValues(data = dat, par = list(mu = bad_mu_vec), spatialCovs = covs),
               "must be a matrix")

  # contains NAs
  bad_mu_na <- matrix(0, nrow = nrow(dat), ncol = 2)
  bad_mu_na[1, 1] <- NA
  expect_error(initialValues(data = dat, par = list(mu = bad_mu_na), spatialCovs = covs),
               "cannot contain missing values")
})

test_that("initialValues rejects observation parameters when data lacks corresponding errors", {
  dat <- get_mock_dataLangevin()
  covs <- get_mock_covs()

  # attempt to pass psi
  expect_error(initialValues(data = dat, par = list(psi = 1.5), spatialCovs = covs),
               "error ellipse observations")

  # attempt to pass tau
  expect_error(initialValues(data = dat, par = list(tau = c(1,1)), spatialCovs = covs),
               "standard error observations")

  # attempt to pass rho_o
  expect_error(initialValues(data = dat, par = list(rho_o = 0.5), spatialCovs = covs),
               "standard error observations")
})

test_that("User-provided par$mu and par$vel override empirical estimates", {
  dat <- get_mock_dataLangevin()
  covs <- get_mock_covs()

  # Create custom matrices
  custom_mu <- matrix(99, nrow = nrow(dat), ncol = 2)
  custom_vel <- matrix(88, nrow = nrow(dat), ncol = 2)

  user_par <- list(mu = custom_mu, vel = custom_vel)

  par_init <- initialValues(data = dat, model = "underdamped", par = user_par, spatialCovs = covs)

  # Ensure the custom matrices bypassed the NA-interpolation and zero-filling
  expect_equal(par_init$mu, custom_mu)
  expect_equal(par_init$vel, custom_vel)
})

test_that("initialValues rejects gamma and vel for the overdamped model", {
  dat <- get_mock_dataLangevin()
  covs <- get_mock_covs()

  # Attempt to pass gamma
  expect_error(initialValues(data = dat, model = "overdamped", par = list(gamma = 1), spatialCovs = covs),
               "model = 'overdamped'")

  # Attempt to pass vel
  custom_vel <- matrix(0, nrow = nrow(dat), ncol = 2)
  expect_error(initialValues(data = dat, model = "overdamped", par = list(vel = custom_vel), spatialCovs = covs),
               "model = 'overdamped'")
})
