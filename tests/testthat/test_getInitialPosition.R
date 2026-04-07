# tests/testthat/test_getInitialPosition.R

# --- Helper to create a mock spatial covariate list ---
get_mock_covs <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, ext = c(0, 10, 0, 10))
  terra::values(r) <- runif(100)
  names(r) <- "cov1"
  return(list(cov1 = r))
}

test_that("getInitialPosition random sampling works when missing", {
  covs <- get_mock_covs()
  beta <- c(0.5)
  nbAnimals <- 5

  # Should output message and return 5x2 matrix
  expect_message(
    pos <- getInitialPosition(nbAnimals = nbAnimals, spatialCovs = covs, beta = beta),
    "Randomly drawing initial positions"
  )

  expect_true(is.matrix(pos))
  expect_equal(nrow(pos), nbAnimals)
  expect_equal(ncol(pos), 2)

  # Check that coordinates fall within the raster extent
  ext_r <- terra::ext(covs[[1]])
  expect_true(all(pos[, 1] >= ext_r$xmin & pos[, 1] <= ext_r$xmax))
  expect_true(all(pos[, 2] >= ext_r$ymin & pos[, 2] <= ext_r$ymax))
})

test_that("getInitialPosition handles a single vector input", {
  covs <- get_mock_covs()
  beta <- c(0.5)
  nbAnimals <- 3

  # Single coordinate pair within the 0-10 extent
  init_vec <- c(5, 5)

  pos <- getInitialPosition(nbAnimals = nbAnimals, initialPosition = init_vec, spatialCovs = covs, beta = beta)

  expect_true(is.matrix(pos))
  expect_equal(nrow(pos), 3)

  # All rows should identical to the input vector
  expect_equal(pos[1, ], init_vec)
  expect_equal(pos[2, ], init_vec)
  expect_equal(pos[3, ], init_vec)
})

test_that("getInitialPosition handles a properly formatted list input", {
  covs <- get_mock_covs()
  beta <- c(0.5)
  nbAnimals <- 2

  init_list <- list(c(2, 2), c(8, 8))

  pos <- getInitialPosition(nbAnimals = nbAnimals, initialPosition = init_list, spatialCovs = covs, beta = beta)

  expect_true(is.matrix(pos))
  expect_equal(nrow(pos), 2)
  expect_equal(pos[1, ], c(2, 2))
  expect_equal(pos[2, ], c(8, 8))
})

test_that("getInitialPosition catches invalid list and vector formats", {
  covs <- get_mock_covs()
  beta <- c(0.5)

  # List of wrong length
  bad_list_len <- list(c(2, 2))
  expect_error(
    getInitialPosition(nbAnimals = 2, initialPosition = bad_list_len, spatialCovs = covs, beta = beta),
    "must be a list of length 2"
  )

  # List containing invalid vectors (length 3 instead of 2)
  bad_list_vec <- list(c(2, 2, 2), c(8, 8))
  expect_error(
    getInitialPosition(nbAnimals = 2, initialPosition = bad_list_vec, spatialCovs = covs, beta = beta),
    "must be a finite numeric vector of length 2"
  )

  # Bad single vector (NAs or characters)
  expect_error(
    getInitialPosition(nbAnimals = 2, initialPosition = c(NA, 5), spatialCovs = covs, beta = beta),
    "must be a finite numeric vector of length 2"
  )
})

test_that("getInitialPosition catches coordinates outside the spatial extent", {
  covs <- get_mock_covs()
  beta <- c(0.5)

  # The mock extent is 0 to 10. This is outside.
  out_of_bounds_vec <- c(15, 15)

  expect_error(
    getInitialPosition(nbAnimals = 2, initialPosition = out_of_bounds_vec, spatialCovs = covs, beta = beta),
    "not within the spatial extent"
  )
})
