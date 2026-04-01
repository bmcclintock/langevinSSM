#' Compute Utilization Distribution
#' @param spatialCovs List of named \code{\link[terra]{SpatRaster-class}} objects containing the spatial covariates. The covariates must be on the same spatial grid and have the same spatial extent.
#' @param beta Numeric vector of habitat selection coefficients for the spatial covariates. The order of the coefficients must match the order of the covariates in \code{spatialCovs}.
#' @param log Logical indicating whether or not to return the log of the utilization distribution. Default: \code{TRUE}.
#' @return A \code{\link[terra]{SpatRaster-class}} object containing the (log) utilization distribution. If the covariates in \code{spatialCovs} have time information (see \code{\link[terra]{time}}), the resulting time-dependent ``utilization distribution'' will also have time information, and each layer will correspond to a different time point.
#' @seealso \code{\link{plotRaster}} for plotting the utilization distribution and covariates.
#'
#' @importFrom terra global nlyr
#' @export
getUD <- function(spatialCovs, beta, log = TRUE) {
  if(length(spatialCovs) != length(beta)) stop("length(spatialCovs) must equal length(beta)")

  ud_rast <- spatialCovs[[1]] * beta[1]
  if(length(spatialCovs) > 1) {
    for (j in 2:length(spatialCovs)) {
      ud_rast <- ud_rast + (spatialCovs[[j]] * beta[j])
    }
  }

  if (!log) {
    # Log-Sum-Exp trick to prevent exp() overflow
    max_log <- terra::global(ud_rast, "max", na.rm = TRUE)$max

    for(k in 1:terra::nlyr(ud_rast)) {
      ud_rast[[k]] <- exp(ud_rast[[k]] - max_log[k])
    }

    # Normalize so the probabilities sum to exactly 1
    layer_sums <- terra::global(ud_rast, "sum", na.rm = TRUE)$sum

    for(k in 1:terra::nlyr(ud_rast)) {
      ud_rast[[k]] <- ud_rast[[k]] / layer_sums[k]
    }
    # ---------------------------------------------------
  }

  n_layers <- sapply(spatialCovs, terra::nlyr)
  max_layers <- max(n_layers)

  if (max_layers > 1) {
    dyn_idx <- which(n_layers == max_layers)[1]
    z_times <- terra::time(spatialCovs[[dyn_idx]])

    if (!is.null(z_times)) {
      names(ud_rast) <- paste0(ifelse(log, "log_", ""),"UD_", z_times)
      terra::time(ud_rast) <- z_times
    } else {
      names(ud_rast) <- paste0(ifelse(log, "log_", ""),"UD_layer_", 1:max_layers)
    }
  } else {
    names(ud_rast) <- paste0(ifelse(log, "log_", ""), "UD")
  }

  return(ud_rast)
}
