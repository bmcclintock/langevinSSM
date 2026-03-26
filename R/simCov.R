#' Simulate spatial covariate
#'
#' Simulate a spatial covariate using a Gaussian random field with a Matérn covariance function. The resulting covariate is on a grid of size (2*sca+1) x (2*sca+1) and has spatial range equal to \code{irange} * \code{sca}.
#' The simulation is performed using the \code{\link[fields]{matern.image.cov}} function from the \pkg{fields} package, which uses FFT padding to efficiently simulate large spatial fields.
#' @param sca Numeric value for the spatial scale of the covariate. The resulting covariate will be on a grid of size (2*sca+1) x (2*sca+1). Default: 100.
#' @param irange Numeric value for the spatial range of the covariate, expressed as a proportion of \code{sca}. The resulting covariate will have spatial range equal to \code{irange} * \code{sca}. Default: 0.3.
#' @param sigma2 Numeric value for the variance of the covariate. Default: 0.1.
#' @param kappa Numeric value for the smoothness of the Matérn covariance function. Default: 0.5 (exponential covariance).
#' @param M Numeric value for the number of rows in the grid used for FFT padding. Default: NULL, in which case the optimal number of columns is dynamically calculated as the next power of 2 greater than or equal to \code{max(2 * n_grid, n_grid + (2 * decay_distance)}, where \code{n_grid} is the number of rows/columns in the original grid (i.e., \code{2*sca+1}) and \code{decay_distance} is the distance at which the covariance decays to near zero (i.e., \code{4 * phi}, where \code{phi = irange * sca} is the spatial range).
#' @param N Numeric value for the number of columns in the grid used for FFT padding. Default: NULL, in which case the optimal number of columns is dynamically calculated as the next power of 2 greater than or equal to \code{max(2 * n_grid, n_grid + (2 * decay_distance)}, where \code{n_grid} is the number of rows/columns in the original grid (i.e., \code{2*sca+1}) and \code{decay_distance} is the distance at which the covariance decays to near zero (i.e., \code{4 * phi}, where \code{phi = irange * sca} is the spatial range).
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

  # Dynamically calculate optimal FFT padding
  if (is.null(M) || is.null(N)) {

    decay_distance <- 4 * phi

    min_pad <- max(2 * n_grid, n_grid + (2 * decay_distance))

    optimal_pad <- 2^ceiling(log2(min_pad))

    if (is.null(M)) M <- optimal_pad
    if (is.null(N)) N <- optimal_pad
  }

  # Define the grid
  grid_list <- list(x = seq(-sca - 0.5, sca + 0.5, length.out = n_grid),
                    y = seq(-sca - 0.5, sca + 0.5, length.out = n_grid))

  # Setup the Matérn covariance object with FFT padding
  obj <- fields::matern.image.cov(setup = TRUE,
                                  grid = grid_list,
                                  theta = phi,
                                  smoothness = kappa,
                                  M = M,
                                  N = N)

  # Simulate and scale by the standard deviation
  grf_fields <- sqrt(sigma2) * fields::sim.rf(obj)

  # Convert to SpatRaster and orient correctly
  spatialCov <- terra::flip(
    terra::rast(
      t(grf_fields),
      extent = terra::ext(min(grid_list$x), max(grid_list$x),
                          min(grid_list$y), max(grid_list$y))
    ),
    direction = "vertical"
  )

  return(spatialCov)
}
