# tests/testthat/test_rasterOverlap.R

get_valid_raster <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 100, 0, 100))
  terra::values(r) <- 1:100
  names(r) <- "habitat"
  return(r)
}

get_valid_par <- function() {
  list(beta = c(0.5), sigma = 5, gamma = 0.5)
}

test_that("simLangevin catches invalid numerical arguments", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  expect_error(simLangevin(par = p, spatialCovs = r, nbAnimals = 0), "should be at least 1")
  expect_error(simLangevin(par = p, spatialCovs = r, nbAnimals = -5), "should be at least 1")
  expect_error(simLangevin(par = p, spatialCovs = r, nbAnimals = Inf), "should be at least 1")

  expect_error(simLangevin(par = p, spatialCovs = r, obsPerAnimal = 0), "should be at least 1")
  expect_error(simLangevin(par = p, spatialCovs = r, obsPerAnimal = -10), "should be at least 1")

  expect_error(simLangevin(par = p, spatialCovs = r, timeStep = 0), "greater than zero")
  expect_error(simLangevin(par = p, spatialCovs = r, timeStep = -0.5), "greater than zero")
})

test_that("simLangevin catches invalid model types", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  expect_error(simLangevin(model = "invalid_model", par = p, spatialCovs = r), "should be one of")
})

test_that("simLangevin returns a correctly structured dataLangevin object", {
  r <- list(habitat = get_valid_raster())
  p <- get_valid_par()

  res <- simLangevin(model = "underdamped", par = p, spatialCovs = r,
                     nbAnimals = 1, obsPerAnimal = 5, timeStep = 0.1,
                     initialPosition = c(50, 50))

  expect_s3_class(res, "dataLangevin")
  expect_true(all(c("id", "date", "dt", "x", "y", "mu.x", "mu.y", "vel.x", "vel.y") %in% names(res)))
  expect_equal(nrow(res), 5)

  res_over <- simLangevin(model = "overdamped", par = p, spatialCovs = r,
                          nbAnimals = 1, obsPerAnimal = 5, timeStep = 0.1,
                          initialPosition = c(50, 50))
  expect_true(all(c("mu.x", "mu.y") %in% names(res_over)))
  expect_false("vel.x" %in% names(res_over))
})

test_that("simLangevin validates natural scale parameters in par", {
  r <- list(habitat = get_valid_raster())

  expect_error(simLangevin(par = list(sigma = 5, gamma = 0.5), spatialCovs = r),
               "beta")

  expect_error(simLangevin(par = list(beta = c(0.5, -0.2), sigma = 5, gamma = 0.5), spatialCovs = r),
               "beta")

  expect_error(simLangevin(par = list(beta = c(0.5), gamma = 0.5), spatialCovs = r),
               "sigma")

  expect_error(simLangevin(model = "underdamped", par = list(beta = c(0.5), sigma = 5), spatialCovs = r),
               "gamma")

  expect_error(simLangevin(par = list(beta = c(0.5), sigma = 5, gamma = 0.5, tau = 1.5), spatialCovs = r),
               "tau")
})
