# tests/testthat/test_simLangevin.R

# --- Helper to create a valid SpatRaster ---
get_valid_raster <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 100, 0, 100))
  terra::values(r) <- 1:100
  names(r) <- "habitat"
  return(r)
}

# --- Helper to create valid parameters ---
get_valid_par <- function() {
  list(beta = c(0.5), sigma = 5, gamma = 0.5)
}

test_that("simLangevin catches invalid numerical arguments", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  # Invalid nbAnimals
  expect_error(simLangevin(par = p, spatialCovs = r, nbAnimals = 0), "should be at least 1")
  expect_error(simLangevin(par = p, spatialCovs = r, nbAnimals = -5), "should be at least 1")
  expect_error(simLangevin(par = p, spatialCovs = r, nbAnimals = Inf), "should be at least 1")

  # Invalid obsPerAnimal
  expect_error(simLangevin(par = p, spatialCovs = r, obsPerAnimal = 0), "should be at least 1")
  expect_error(simLangevin(par = p, spatialCovs = r, obsPerAnimal = -10), "should be at least 1")

  # Invalid timeStep
  expect_error(simLangevin(par = p, spatialCovs = r, timeStep = 0), "greater than zero")
  expect_error(simLangevin(par = p, spatialCovs = r, timeStep = -0.5), "greater than zero")
})

test_that("simLangevin catches invalid model types", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  # match.arg will catch this natively
  expect_error(simLangevin(model = "invalid_model", par = p, spatialCovs = r), "should be one of")
})

test_that("simLangevin returns a correctly structured dataLangevin object", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  # Run a very tiny simulation (5 steps) to test structure without slowing down tests
  res <- simLangevin(model = "underdamped", par = p, spatialCovs = r,
                     nbAnimals = 1, obsPerAnimal = 5, timeStep = 0.1)

  expect_s3_class(res, "dataLangevin")
  expect_true(all(c("id", "date", "dt", "x", "y", "mu.x", "mu.y", "vel.x", "vel.y") %in% names(res)))
  expect_equal(nrow(res), 5)

  # Test overdamped output structure (should lack vel.x and vel.y)
  res_over <- simLangevin(model = "overdamped", par = p, spatialCovs = r,
                          nbAnimals = 1, obsPerAnimal = 5, timeStep = 0.1)
  expect_true(all(c("mu.x", "mu.y") %in% names(res_over)))
  expect_false("vel.x" %in% names(res_over))
})

test_that("simLangevin validates natural scale parameters in par", {
  r <- list(habitat = get_valid_raster()) # 1 covariate

  # Missing beta entirely
  expect_error(simLangevin(par = list(sigma = 5, gamma = 0.5), spatialCovs = r),
               "beta")

  # beta length mismatch (provided 2 coefficients, but only 1 spatial covariate exists)
  expect_error(simLangevin(par = list(beta = c(0.5, -0.2), sigma = 5, gamma = 0.5), spatialCovs = r),
               "beta")

  # Missing sigma
  expect_error(simLangevin(par = list(beta = c(0.5), gamma = 0.5), spatialCovs = r),
               "sigma")

  # Missing gamma (required for the default underdamped model)
  expect_error(simLangevin(model = "underdamped", par = list(beta = c(0.5), sigma = 5), spatialCovs = r),
               "gamma")

  # Invalid tau length (must be a 2-vector for x and y standard deviations)
  expect_error(simLangevin(par = list(beta = c(0.5), sigma = 5, gamma = 0.5, tau = 1.5), spatialCovs = r),
               "tau")
})
