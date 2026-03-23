Rcpp::sourceCpp("src/simulateLangevin.cpp")

rasterList <- function (rast) 
{
  lim <- as.vector(terra::ext(rast))
  res <- terra::res(rast)
  xgrid <- seq(lim[1] + res[1]/2, lim[2] - res[1]/2, by = res[1])
  ygrid <- seq(lim[3] + res[2]/2, lim[4] - res[2]/2, by = res[2])
  
  # Add wide = TRUE so terra returns [nrow, ncol] instead of [ncells, nlayers]
  z_mat <- terra::as.matrix(rast, wide = TRUE)
  z <- t(apply(z_mat, 2, rev))
  
  return(list(x = xgrid, y = ygrid, z = z))
}

getUD <- function (spatialCovs, beta, log = F) 
{
  if(length(spatialCovs)!=length(beta)) stop("length(spatialCovs) must equal length(beta)")
  covariates <- lapply(spatialCovs,rasterList)
  ud_rast <- covariates[[1]]
  dx <- diff(ud_rast$x)[1]
  dy <- diff(ud_rast$y)[1]
  J <- length(covariates)
  ud_rast$z <- Reduce("+", lapply(1:J, function(j) dx * dy * 
                                    beta[j] * covariates[[j]]$z))
  if (!log) {
    ud_rast$z <- exp(ud_rast$z)
    ud_rast$z <- ud_rast$z/sum(ud_rast$z)
  }
  terra::rast(ud_rast)
}

plotRaster <- function (rast, norm = FALSE, log = FALSE, scale.name = "", light = FALSE) 
{
  # FIX 1: Use as.numeric() to ensure 'val' is a simple vector, not a matrix column
  covmap <- data.frame(terra::crds(rast), val = as.numeric(terra::values(rast)))
  
  if (norm) {
    # FIX 2: Use native terra::res() instead of raster::xres/yres
    r_res <- terra::res(rast)
    s <- sum(covmap$val) * r_res[1] * r_res[2]
    covmap$val <- covmap$val / s
  }
  
  if (log) {
    covmap$val <- log(covmap$val)
  }
  
  # FIX 3: Replace aes_string() with standard aes()
  p <- ggplot(covmap, aes(x = x, y = y)) + 
    geom_raster(aes(fill = val)) + 
    coord_equal()
  
  if (light) {
    p <- p + scale_fill_viridis(guide = "none") + 
      theme(axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
  } else {
    p <- p + scale_fill_viridis(name = scale.name)
  }
  
  return(p)
}

simCov <- function(sca = 200, irange = 0.3, sigma2 = 0.1, kappa = 0.5, M = 2048, N = 2048) {
  
  # --- Parameters ---
  phi <- irange * sca  
  n_grid <- 2 * sca + 1
  
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

radian <- function(degree) {
  radian <- degree * (pi/180)
  return(radian)
}

formatData <- function(data, id = "id", date = "date", coord = c("x", "y"), lc = "lc", epar = c("smaj", "smin", "eor"), sderr = c("x.sd", "y.sd"), time.unit = "hours", tz = "UTC"){
  
  if(is.null(data[[lc]])) stop("data$lc is missing")
  if(!all(data[[lc]] %in% c(3,2,1,0,"A","B","Z","G"))) stop("data$lc can only be '3', '2', '1', '0', 'A', 'B', and 'Z' for location quality classes or 'G' for GPS data")
  if(!inherits(data[[date]],"POSIXt")) stop("data$date must be of class 'POSIXt'")
  out <- aniMotum::format_data(x = data, id = id, date = date, coord = coord, lc = lc, epar = epar, sderr = sderr, tz = tz)
  out <- out %>% dplyr::rename(x=lon,y=lat) %>% 
                 dplyr::arrange(id,date) %>%
                 dplyr::mutate(dt=do.call(c,mapply(function(x) c(0,diff(out$date[which(out$id==x)],units=time.unit)), unique(out$id), SIMPLIFY = FALSE))) %>% 
                 dplyr::select(id,date,dt,x,y,smaj,smin,eor,x.sd,y.sd)
  attr(out,"row.names") <- NULL
  attr(out,"time.unit") <- time.unit
  return(out)
}

prepareRaster <- function(spatialCovs,scale_factor=1,time.unit="hours"){
  
  if(!is.list(spatialCovs)) stop('spatialCovs must be a list')
  spatialcovnames <- names(spatialCovs)
  if(is.null(spatialcovnames)) stop('spatialCovs must be a named list')
  nbSpatialCovs <- length(spatialcovnames)
  
  # 1. Update package check to terra
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package \"terra\" needed for spatial covariates. Please install it.", call. = FALSE)
  }
  
  for(j in 1:nbSpatialCovs){
    # 2. Update class check (everything is a SpatRaster)
    if(!inherits(spatialCovs[[j]], "SpatRaster")) {
      stop("spatialCovs$", spatialcovnames[j], " must be of class 'SpatRaster'")
    }
    
    # 3. Update value extraction
    if(any(is.na(terra::values(spatialCovs[[j]])))) {
      stop("missing values are not permitted in spatialCovs$", spatialcovnames[j])
    }
    
    # 4. Update multi-layer (stack/brick) check and Z-value logic
    if(terra::nlyr(spatialCovs[[j]]) > 1) {
      
      t_vals <- terra::time(spatialCovs[[j]])
      
      # Check if time/Z values are set
      if(is.null(t_vals) || all(is.na(t_vals))) {
        stop("spatialCovs$", spatialcovnames[j], " is a multi-layer raster that must have time values set (see ?terra::time)")
      } 
      
      # Since terra doesn't store a "Z name" mapping, check directly against the expected 
      # time column in your data. (Change "time" to whatever your standard column name is).
      else if(!("time" %in% names(data))) {
        stop("spatialCovs$", spatialcovnames[j], " requires a 'time' column in 'data' to match the raster's dynamic layers")
      }
    }
  }
  
  if(any(spatialcovnames %in% names(data))) stop("spatialCovs cannot have same names as data")
  if(anyDuplicated(spatialcovnames)) stop("spatialCovs must have unique names")
  
  # Extract the 3D array [nrow, ncol, nlayer]
  vals_array <- as.array(raster_stack) 
  
  # CORRECTED: Permute to [ncol, nrow, nlayer] so it perfectly matches the C++ idx formula
  vals_array <- aperm(vals_array, c(2, 1, 3))
  
  # Extract and convert raster times
  times_list <- lapply(spatialCovs, function(r) {
    t_vals <- terra::time(r)
    
    if (is.null(t_vals) || all(is.na(t_vals))) {
      if (terra::nlyr(r) == 1) {
        return(0) # Static covariate dummy value
      } else {
        # Fallback for dynamic rasters without time metadata
        return(as.numeric(1:terra::nlyr(r)))
      }
    } else {
      # Check if the time values are datetime or date objects
      if (inherits(t_vals, "POSIXt") || inherits(t_vals, "Date")) {
        return(as.numeric(difftime(t_vals, as.POSIXct("1970-01-01 00:00:00", tz = "UTC"), units = time.unit)))
      } else {
        # If they are already numeric (like in your simulations), return as-is
        return(as.numeric(t_vals))
      }
    }
  })
  
  all_z_values_R <- unlist(times_list)
  
  n_zvals_cov <- sapply(spatialCovs, nlyr) 
  cov_offset_R <- c(0, cumsum(n_zvals_cov)[-length(n_zvals_cov)])
  
  all_z_values_R <- unlist(times_list)
  
  raster_data <- list(
    raster_vals = vals_array,
    raster_coords = terra::crds(raster_stack)/scale_factor, 
    raster_resolution = terra::res(raster_stack)/scale_factor,
    raster_extent = as.vector(terra::ext(raster_stack)/scale_factor),
    n_covs = length(n_zvals_cov), 
    all_z_values = as.numeric(all_z_values_R), # The flattened raster slice times
    n_zvals_cov = as.integer(n_zvals_cov),
    cov_offset = as.integer(cov_offset_R)
  )
  return(raster_data)
}

measurementError <- function(data,M,m,c,psi){
  M <- tmpM <- abs(rnorm(nbAnimals*obsPerAnimal,0,sd=M))
  m <- tmpm <- abs(rnorm(nbAnimals*obsPerAnimal,0,sd=m))
  M[which(tmpM < tmpm)] <- tmpm[which(tmpM < tmpm)]
  m[which(tmpM < tmpm)] <- tmpM[which(tmpM < tmpm)]
  c <- radian(runif(nbAnimals*obsPerAnimal,c[1],c[2]))
  z = sqrt(2);
  s2c = sin(c) * sin(c);
  c2c = cos(c) * cos(c);
  M2  = (M / z) * (M / z);
  m2 = (m * psi / z) * (m * psi / z);
  
  data$mux <- data$mu.x
  data$muy <- data$mu.y
  data$mu.x <- NA
  data$mu.y <- NA
  data$error_semimajor_axis <- M
  data$error_semiminor_axis <- m
  data$error_ellipse_orientation <- c
  
  for(i in 1:nrow(data)){
    cov_obs <- matrix(0,2,2)
    cov_obs[1,1] = (M2[i] * s2c[i] + m2[i] * c2c[i]);
    cov_obs[2,2] = (M2[i] * c2c[i] + m2[i] * s2c[i]);
    cov_obs[1,2] = (0.5 * (M[i] * M[i] - (m[i] * psi * m[i] * psi))) * cos(c[i]) * sin(c[i]);
    cov_obs[2,1] = cov_obs[1,2];
    mu <- mvtnorm::rmvnorm(1,cbind(data$mux[i],data$muy[i]),cov_obs)
    data$mu.x[i] <- mu[1]
    data$mu.y[i] <- mu[2]
  }
  return(data)
}

init.mu_aniMotum <- function(subDat,model="rw",timeSteps){
  
  aniDat <- data.frame(
    id = as.character(subDat$ID),
    date = as.POSIXlt(subDat$time * 1/mean(subDat$dt)), # fit_ssm doesn't like very small \Delta_t
    x = subDat$mu.x,
    y = subDat$mu.y,
    lc = 3,
    smaj = subDat$error_semimajor_axis,
    smin = subDat$error_semiminor_axis,
    eor = subDat$error_ellipse_orientation,
    x.sd = NA,
    y.sd = NA
  )
  aniDat <- sf::st_as_sf(
    x = aniDat,
    coords = c("x", "y"),
    crs = 3416,  
    remove = FALSE
  )
  
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
  for(j in which(unlist(lapply(ssm_results,is.null)) | unlist(lapply(ssm_results,function(x) !isTRUE(x$converged))))){
    message("      aniMotum failed for individual ",unique_ids[j],"; trying crawl instead")
    locErr <- crawl::argosDiag2Cov(id_data_list[[j]]$smaj,id_data_list[[j]]$smin,id_data_list[[j]]$eor/(pi/180))
    id_data_list[[j]]$ln.sd.x <- locErr$ln.sd.x
    id_data_list[[j]]$ln.sd.y <- locErr$ln.sd.y
    id_data_list[[j]]$error.corr <- locErr$error.corr
    predTime <- list()
    predTime[[unique_ids[j]]] <- timeSteps$date[which(timeSteps$id==j)]
    crfit <- tryCatch(momentuHMM::crawlWrap(id_data_list[[j]] %>% rename(ID=id,time=date),predTime=predTime,err.model = list(x =  ~ ln.sd.x - 1,y =  ~ ln.sd.y - 1, rho =  ~ error.corr),fixPar=c(1,1,NA,NA),initialSANN=NULL,retryFits=5),error=function(e) e)
    if(!inherits(crfit,"error")) ssm_results[[j]] <- crfit$crwPredict[which(crfit$crwPredict$locType=="p"),]
    else stop("      crawl also failed for individual ",unique_ids[j],": ",crfit$message)
  }
  
  # Combine results
  aniFit <- list()
  for (i in seq_along(unique_ids)) {
    if(inherits(ssm_results[[i]],"ssm_df")) aniFit[[i]] <- matrix(unlist(ssm_results[[i]]$ssm[[1]]$predicted$geometry),nrow=2)
    else aniFit[[i]] <- matrix(unlist(ssm_results[[i]][,c("mu.x","mu.y")]),nrow=2)
  }

  # Extract initial mu as before
  init.mu <- do.call(cbind, aniFit)
  
  return(init.mu)
}
