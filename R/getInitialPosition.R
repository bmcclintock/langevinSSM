#' @importFrom terra ncell values xyFromCell cellFromXY
getInitialPosition <- function(nbAnimals, initialPosition, spatialCovs, beta, replace = FALSE) {

  if (missing(initialPosition)) {
    message("   Randomly drawing initial positions from UD...")
    UD <- getUD(spatialCovs, beta = beta, log = TRUE)

    log_ud_vals <- as.numeric(terra::values(UD[[1]]))
    max_log_ud <- max(log_ud_vals, na.rm = TRUE)

    prob_vals <- exp(log_ud_vals - max_log_ud)
    prob_vals[is.na(prob_vals)] <- 0
    prob_vals <- prob_vals / sum(prob_vals)

    if(terra::ncell(UD[[1]]) < nbAnimals) replace = TRUE
    sampled_cells <- sample(terra::ncell(UD[[1]]), nbAnimals, replace = replace, prob = prob_vals)
    initialPosition <- terra::xyFromCell(UD[[1]], sampled_cells)

  } else {

    if (is.list(initialPosition)) {
      if (length(initialPosition) != nbAnimals) stop("initialPosition must be a list of length ", nbAnimals)
      for (i in 1:nbAnimals) {
        if (length(initialPosition[[i]]) != 2 | !is.numeric(initialPosition[[i]]) | any(!is.finite(initialPosition[[i]]))) {
          stop("each element of initialPosition must be a finite numeric vector of length 2")
        }
      }
      initialPosition <- do.call(rbind, initialPosition)

    } else {
      if (length(initialPosition) != 2 | !is.numeric(initialPosition) | any(!is.finite(initialPosition))) {
        stop("initialPosition must be a finite numeric vector of length 2")
      }
      initialPosition <- matrix(initialPosition, nrow = nbAnimals, ncol = 2, byrow = TRUE)
    }

    cells <- terra::cellFromXY(spatialCovs[[1]], initialPosition)
    if (any(is.na(cells))) {
      bad_idx <- which(is.na(cells))[1]
      stop("initialPosition for individual ", bad_idx, " is not within the spatial extent of the rasters")
    }
  }

  # ensure matrix columns aren't named to prevent downstream TMB/C++ mismatches
  colnames(initialPosition) <- NULL
  return(initialPosition)
}
