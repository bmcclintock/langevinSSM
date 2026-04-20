#' Simulate spatial covariate
#'
#' Simulate a spatial covariate using a Gaussian random field with a Matérn covariance function. The resulting covariate is on a grid of size (2*sca+1) x (2*sca+1) and has spatial range equal to \code{irange} * \code{sca}.
#' The simulation is performed using the \code{\link[fields]{matern.image.cov}} function from the \pkg{fields} package, which uses FFT padding to efficiently simulate large spatial fields.
#' @param sca Numeric value for the spatial scale of the covariate. The resulting covariate will be on a grid of size (2*sca+1) x (2*sca+1). Default: 100.
#' @param irange Numeric value for the spatial range of the covariate, expressed as a proportion of \code{sca}. The resulting covariate will have spatial range equal to \code{irange} * \code{sca}. Default: 0.3.
#' @param sigma2 Numeric value for the variance of the covariate. Default: 0.1.
#' @param kappa Numeric value for the smoothness of the Matérn covariance function. Default: 0.5 (exponential covariance).
#' @param M Numeric value for the number of rows in the grid used for FFT padding. Default: NULL, in which case the optimal number of columns is dynamically calculated.
#' @param N Numeric value for the number of columns in the grid used for FFT padding. Default: NULL, in which case the optimal number of columns is dynamically calculated.
#' @return A \code{\link[terra]{SpatRaster-class}} object containing the simulated spatial covariate.
#' @export
simCov <- function(sca = 100,
                   irange = 0.3,
                   sigma2 = 0.1,
                   kappa = 0.5,
                   M = NULL,
                   N = NULL) {

  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Package \"fields\" needed for generating spatial covariates. Please install it.", call. = FALSE)
  }

  # --- Parameters ---
  phi <- irange * sca
  n_grid <- 2 * sca + 1

  grid_list <- list(x = seq(-sca - 0.5, sca + 0.5, length.out = n_grid),
                    y = seq(-sca - 0.5, sca + 0.5, length.out = n_grid))

  # Dynamically calculate initial baseline FFT padding
  decay_distance <- 4 * phi
  min_pad <- max(2 * n_grid, n_grid + (2 * decay_distance))
  base_pad <- 2^ceiling(log2(min_pad))

  current_M <- if (is.null(M)) base_pad else M
  current_N <- if (is.null(N)) base_pad else N

  success <- FALSE
  attempt <- 1
  max_attempts <- 4

  # Retry loop to exponentially expand padding if negative FFT eigenvalues occur
  while (!success && attempt <= max_attempts) {
    obj <- fields::matern.image.cov(setup = TRUE,
                                    grid = grid_list,
                                    theta = phi,
                                    smoothness = kappa,
                                    M = current_M,
                                    N = current_N)

    grf_raw <- try(fields::sim.rf(obj), silent = TRUE)

    if (!inherits(grf_raw, "try-error")) {
      success <- TRUE
    } else {
      # If the FFT has negative values, double the padding grid and try again
      attempt <- attempt + 1
      if (is.null(M)) current_M <- current_M * 2
      if (is.null(N)) current_N <- current_N * 2
    }
  }

  if (!success) {
    stop("fields::sim.rf failed: FFT of covariance has negative values even after maximum padding. Try reducing 'irange' or 'sca'.")
  }

  # Scale by the standard deviation
  grf_fields <- sqrt(sigma2) * grf_raw

  # Transpose the matrix to match spatial orientation
  grf_matrix <- t(grf_fields)
  grf_flipped <- grf_matrix[nrow(grf_matrix):1, ]

  # Convert directly to a standalone SpatRaster
  spatialCov <- terra::rast(
    grf_flipped,
    extent = terra::ext(min(grid_list$x), max(grid_list$x),
                        min(grid_list$y), max(grid_list$y))
  )

  return(spatialCov)
}
