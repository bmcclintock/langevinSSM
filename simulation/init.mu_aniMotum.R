init.mu_aniMotum <- function(subDat,model="rw",timeSteps){

  aniDat <- data.frame(
    id = as.character(subDat$id),
    date = as.POSIXlt(subDat$date * 1000), # fit_ssm doesn't like very small \Delta_t
    x = subDat$x,
    y = subDat$y,
    lc = "G",
    smaj = subDat$smaj,
    smin = subDat$smin,
    eor = subDat$eor,
    x.sd = NA,
    y.sd = NA
  )
  aniDat <- sf::st_as_sf(
    x = aniDat,
    coords = c("x", "y"),
    crs = 3416,
    remove = FALSE
  )

  # prevent automatic km conversion by fit_ssm
  suppressWarnings(sf::st_crs(aniDat) <- "+proj=lcc +lat_0=47.5 +lon_0=13.3333333333333 +lat_1=49 +lat_2=46 +x_0=400 +y_0=400 +ellps=GRS80 +units=km +no_defs")

  # Get unique IDs
  unique_ids <- unique(aniDat$id)

  # Create a list to store all the data needed for each ID
  id_data_list <- list()

  for (id in unique_ids) {
    # Filter data for this ID
    id_indices <- which(aniDat$id == id)

    id_data_list[[id]] <- aniDat[id_indices,]
  }

  # Function to fit SSM for a single ID
  fit_single_ssm <- function(id_data) {

    # Convert to sf
    id_data_sf <- sf::st_as_sf(
      x = id_data,
      coords = c("x", "y"),
      crs = 3416,
      remove = FALSE
    )

    # Fit the model
    tryCatch({
      result <- aniMotum::fit_ssm(
        id_data_sf,
        spdf = TRUE,
        model = model,
        time.step = timeSteps,
        map = list(psi = factor(NA))
      )
      return(result)
    }, error = function(e) {
      message("Error processing ID: ", unique(id_data$id), " - ", e$message)
      return(NULL)
    })
  }

  # Fit models in parallel
  ssm_results <- tryCatch(furrr::future_map(unique_ids,function(x) fit_single_ssm(id_data_list[[x]]), .options = furrr::furrr_options(seed = TRUE)),error=function(e) e)

  # Combine results
  aniFit <- list()
  for (i in seq_along(unique_ids)) {
    if(inherits(ssm_results[[i]],"ssm_df")) {
      aniFit[[i]] <- matrix(unlist(ssm_results[[i]]$ssm[[1]]$predicted$geometry),nrow=2)
    } else {
      message("      aniMotum failed for individual ",unique_ids[i],"; using true locations instead instead")
      aniFit[[i]] <- matrix(unlist(ssm_results[[i]][,c("mu.x","mu.y")]),nrow=2)
    }
  }

  # Extract initial mu as before
  init.mu <- t(do.call(cbind, aniFit))

  return(init.mu)
}
