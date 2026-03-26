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

#' @export
getUD <- function(spatialCovs, beta, log = FALSE) {
  if(length(spatialCovs) != length(beta)) stop("length(spatialCovs) must equal length(beta)")
  
  # Get cell area (dx * dy)
  r_res <- terra::res(spatialCovs[[1]])
  cell_area <- r_res[1] * r_res[2]
  
  ud_rast <- spatialCovs[[1]] * beta[1] * cell_area
  for (j in 2:length(spatialCovs)) {
    ud_rast <- ud_rast + (spatialCovs[[j]] * beta[j] * cell_area)
  }
  
  if (!log) {
    ud_rast <- exp(ud_rast)
    
    layer_sums <- terra::global(ud_rast, "sum", na.rm = TRUE)$sum
    
    for(k in 1:terra::nlyr(ud_rast)) {
      ud_rast[[k]] <- ud_rast[[k]] / layer_sums[k]
    }
  }
  
  n_layers <- sapply(spatialCovs, terra::nlyr)
  max_layers <- max(n_layers)
  
  if (max_layers > 1) {
    
    dyn_idx <- which(n_layers == max_layers)[1]
    z_times <- terra::time(spatialCovs[[dyn_idx]])
    
    if (!is.null(z_times)) {
      names(ud_rast) <- paste0("UD_time_", z_times)
      terra::time(ud_rast) <- z_times
    } else {
      names(ud_rast) <- paste0("UD_layer_", 1:max_layers)
    }
  } else {

    names(ud_rast) <- "UD_static"
  }
  
  return(ud_rast)
}

#' @export
plotRaster <- function (rast, norm = FALSE, log = FALSE, scale.name = expression(pi), light = FALSE) 
{
  covmap <- data.frame(terra::crds(rast), val = as.numeric(terra::values(rast)))
  
  if (norm) {
    r_res <- terra::res(rast)
    s <- sum(covmap$val) * r_res[1] * r_res[2]
    covmap$val <- covmap$val / s
  }
  
  if (log) {
    covmap$val <- log(covmap$val)
  }
  
  p <- ggplot2::ggplot(covmap, aes(x = x, y = y)) + 
    ggplot2::geom_raster(aes(fill = val)) + 
    ggplot2::coord_equal()
  
  if (light) {
    p <- p + viridis::scale_fill_viridis(guide = "none") + 
      ggplot2::theme(axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
  } else {
    p <- p + viridis::scale_fill_viridis(name = scale.name)
  }
  
  return(p)
}

#' @export
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

#' @export
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
  class(out) <- append(class(out),"dataLangevin")
  return(out)
}

prepareRaster <- function(spatialCovs, scaleFactor=1, time.unit="hours", data = NULL) {
  
  if(!is.list(spatialCovs)) stop('spatialCovs must be a list')
  spatialcovnames <- names(spatialCovs)
  if(is.null(spatialcovnames)) stop('spatialCovs must be a named list')
  nbSpatialCovs <- length(spatialcovnames)
  
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package \"terra\" needed for spatial covariates. Please install it.", call. = FALSE)
  }
  
  raster_stack <- terra::rast(spatialCovs)
  
  for(j in 1:nbSpatialCovs) {
    if(!inherits(spatialCovs[[j]], "SpatRaster")) {
      stop("spatialCovs$", spatialcovnames[j], " must be of class 'SpatRaster'")
    }
    
    if(any(is.na(terra::values(spatialCovs[[j]])))) {
      stop("missing values are not permitted in spatialCovs$", spatialcovnames[j])
    }
    
    if(terra::nlyr(spatialCovs[[j]]) > 1) {
      
      t_vals <- terra::time(spatialCovs[[j]])
      
      # Check if time/Z values are set
      if(is.null(t_vals) || all(is.na(t_vals))) {
        stop("spatialCovs$", spatialcovnames[j], " is a multi-layer raster that must have time values set (see ?terra::time)")
      } 
      
      else if(!is.null(data) && !("date" %in% names(data))) {
        stop("spatialCovs$", spatialcovnames[j], " requires a 'date' column in 'data' to match the raster's dynamic layers")
      }
    }
  }
  
  if(!is.null(data) && any(spatialcovnames %in% names(data))) stop("spatialCovs cannot have same names as data")
  if(anyDuplicated(spatialcovnames)) stop("spatialCovs must have unique names")
  
  vals_array <- as.array(raster_stack) 
  
  # Permute to [ncol, nrow, nlayer] so it matches the C++ idx formula
  vals_array <- aperm(vals_array, c(2, 1, 3))
  
  # Extract and convert raster times
  times_list <- lapply(spatialCovs, function(r) {
    t_vals <- terra::time(r)
    
    if (terra::nlyr(r) == 1 && (is.null(t_vals) || all(is.na(t_vals)))) {
      return(0) 
    } else {
      if (inherits(t_vals, "POSIXt") || inherits(t_vals, "Date")) {
        return(as.numeric(difftime(t_vals, as.POSIXct("1970-01-01 00:00:00", tz = "UTC"), units = time.unit)))
      } else {
        return(as.numeric(t_vals))
      }
    }
  })
  
  n_zvals_cov <- sapply(spatialCovs, terra::nlyr) 
  cov_offset_R <- c(0, cumsum(n_zvals_cov)[-length(n_zvals_cov)])
  
  all_z_values_R <- unlist(times_list)
  
  raster_data <- list(
    raster_vals = vals_array,
    raster_coords = terra::crds(raster_stack)/scaleFactor, 
    raster_resolution = terra::res(raster_stack)/scaleFactor,
    raster_extent = as.vector(terra::ext(raster_stack)/scaleFactor),
    n_covs = length(n_zvals_cov), 
    all_z_values = as.numeric(all_z_values_R), # The flattened raster slice times
    n_zvals_cov = as.integer(n_zvals_cov),
    cov_offset = as.integer(cov_offset_R)
  )
  return(raster_data)
}

checkPar <- function(par, model, map=NULL, dat=NULL){
  
  if(is.null(map)) map <- list()
  if(!is.list(par)) stop("par must be a list.")
  if(!all(names(par) %in% c("beta","sigma","gamma","mu","v_mu","psi","tau","rho_o"))) stop("names(par) is limited to c('beta','sigma','gamma','mu','v_mu','psi','tau','rho_o')")
  if(model=="underdamped"){
    if(!is.null(par$gamma)) gamma <- par$gamma
    else stop("par$gamma is missing, with no default.")
    if(!is.finite(gamma) || gamma<=0){
      stop("gamma should be greater than zero.")
    } else {
      par$log_gamma <- log(par$gamma)
      par$gamma <- NULL
    }
    if(!is.null(map$gamma)){
      if(length(map$gamma)!=length(gamma)) stop("map$gamma should be of length 1.")
      map$log_gamma <- map$gamma
      map$gamma <- NULL
    }
  } else {
    map$log_gamma <- factor(NA)
    map$gamma <- NULL
    par$log_gamma <- 0
  }
  if(!is.null(par$sigma)){
    sigma <- par$sigma
    if(!is.finite(sigma) || sigma<=0){
      stop("sigma should be greater than zero.")
    } else {
      par$log_sigma <- log(par$sigma)
      par$sigma <- NULL
    }
    if(!is.null(map$sigma)){
      if(length(map$sigma)!=length(sigma)) stop("map$sigma should be of length 1.")
      map$log_sigma <- map$sigma
      map$sigma <- NULL
    }
  } else stop("par$sigma is missing, with no default.")
  
  if(!is.null(par$beta)){
    beta <- par$beta
    if(any(!is.finite(beta))){
      stop("beta must be finite.")
    }
    if(!is.null(map$beta)){
      if(length(map$beta)!=length(beta)) stop("map$beta should be of length ",length(beta),".")
    }
  } else stop("par$beta is missing, with no default.")
  
  if(!is.null(par$psi)){
    if(!is.finite(par$psi) || psi<=0)
      stop("par$psi should be greater than zero.")
    par$l_psi <- log(par$psi)
    par$psi <- NULL
    if(!is.null(map$psi)){
      if(length(map$psi)!=length(par$l_psi)) stop("map$psi should be of length 1.")
      map$l_psi <- map$psi
      map$psi <- NULL
    }
  } else {
    par$l_psi <- 0
    map$l_psi <- factor(NA)
  }
  if(!is.null(par$tau)){
    if(length(par$tau)!=2) stop("par$tau should be of length 2.")
    if(!is.finite(par$tau) || tau<=0)
      stop("par$tau should be greater than zero.")
    par$l_tau <- log(par$tau)
    par$tau <- NULL
    if(!is.null(map$tau)){
      if(length(map$tau)!=length(par$l_tau)) stop("map$tau should be of length 2.")
      map$l_tau <- map$tau
      map$tau <- NULL
    }
  } else {
    par$l_tau <- c(0,0)
    map$l_tau <- factor(rep(NA,2))
  }
  if(!is.null(par$rho_o)){
    if(!is.finite(par$rho_o)) stop("par$rho_o must be finite.")
    par$l_rho_o <- log(par$rho_o)
    par$rho_o <- NULL
    if(!is.null(map$rho_o)){
      if(length(map$rho_o)!=length(par$l_rho_o)) stop("map$rho_o should be of length 1.")
      map$l_rho_o <- map$rho_o
      map$rho_o <- NULL
    }
  } else {
    par$l_rho_o <- 0
    map$l_rho_o <- factor(NA)
  }
  
  if(!is.null(dat)){
    
    if(!is.null(par$mu)) par$mu <- t(par$mu)
    if(!is.null(par$v_mu)) par$v_mu <- t(par$v_mu)
    
    if(model=="overdamped"){
      re <- "mu"
      par$v_mu <- matrix(0,2,ncol(dat$Y)) 
      map$v_mu <- factor(rep(NA,length(dat$Y)))
      map$log_gamma <- factor(NA)
    } else {
      re <- c("mu","v_mu")
    }
    
    if(all(is.na(dat$obs_mod))){
      if(model=="overdamped") re <- NULL
      else re <- "v_mu"
      par$mu <- dat$Y
      map$mu <- factor(rep(NA,length(dat$Y)))
    } 
    
    if(is.null(par$mu)) stop("par$mu is missing, with no default.")
    if(is.null(par$v_mu)) stop("par$v_mu is missing, with no default.")
    if(any(dim(par$mu)!=dim(dat$Y))) stop("par$mu must have ",ncol(dat$Y)," rows and 2 columns")
    if(any(dim(par$v_mu)!=dim(dat$Y))) stop("par$v_mu must have ",ncol(dat$Y)," rows and 2 columns")
    if(any(!is.finite(par$mu)) | any(!is.finite(par$v_mu))) stop("par$mu and/or par$v_mu must be finite")
    
    par <- par[c("beta","log_sigma","log_gamma","mu","v_mu","l_psi","l_tau","l_rho_o")]
    out <- list(par=par,map=map,re=re)
  } else out <- list(par=par,map=map)
  
  return(out)
}

getInitialPosition <- function(nbAnimals,initialPosition,spatialCovs){
  
  spatialcovnames <- names(spatialCovs)
  
  if(missing(initialPosition)){
    message("   Randomly drawing initial positions from UD...")
    UD <- getUD(spatialCovs, beta=beta,log=TRUE)
    initPos <- matrix(sample(terra::ncell(UD[[1]]),nbAnimals,replace=FALSE,prob=exp(terra::values(UD[[1]]))/sum(exp(terra::values(UD[[1]])))),
                      1,nbAnimals,byrow=TRUE)
    initialPosition <- t(mapply(function(x) xyFromCell(UD,initPos[,x]),1:nbAnimals,SIMPLIFY = FALSE))
  } else {
    if(is.list(initialPosition)){
      if(length(initialPosition)!=nbAnimals) stop("initialPosition must be a list of length ",nbAnimals)
      for(i in 1:nbAnimals){
        if(length(initialPosition[[i]])!=2 | !is.numeric(initialPosition[[i]]) | any(!is.finite(initialPosition[[i]]))) stop("each element of initialPosition must be a finite numeric vector of length 2")
      }
    } else {
      if(length(initialPosition)!=2 | !is.numeric(initialPosition) | any(!is.finite(initialPosition))) stop("initialPosition must be a finite numeric vector of length 2")
      tmpPos<-initialPosition
      initialPosition<-vector('list',nbAnimals)
      for(i in 1:nbAnimals){
        initialPosition[[i]]<-tmpPos
      }
    }
    for(i in 1:length(initialPosition)){
      for(j in 1:length(spatialCovs)){
        if(is.na(terra::cellFromXY(spatialCovs[[j]],matrix(initialPosition[[i]],ncol=2)))) stop("initialPosition for individual ",i," is not within the spatial extent of the ",spatialcovnames[j]," raster")
      }
    }
  }
  initialPosition <- do.call(rbind,initialPosition)
  return(initialPosition)
}

addMeasurementError <- function(model, out, par, measurementError){
  if(!is.null(measurementError)){
    if(!is.list(measurementError)) stop("'measurementError' must be a list.")
    if(!is.null(measurementError$M)){
      M <- measurementError$M
      if(!is.finite(M) || M<=0)
        stop("measurementerror$M should be greater than zero.")
    } 
    if(!is.null(measurementError$m)){
      m <- measurementError$m
      if(!is.finite(m) || m<=0)
        stop("measurementerror$m should be greater than zero.")
    } 
    if(!is.null(measurementError$c)){
      r <- measurementError$c
      if(length(r)!=2) stop("measurementError$c should be of length 2.")
      if(any(!is.finite(r))) stop("measurementError$c must be finite.")
      if(any(r < 0 | r > 180)) stop("measurementError$c must be in degrees between 0 and 180.")
      if(r[2] < r[1]) stop("measurementError$c[2] must be greater than measurementError$c[1].")
    } 
    if(!is.null(par$l_psi)){
      psi <- exp(par$l_psi)
    } else psi <- 1
    out <- measurementError_rcpp(out,
                                 M=M,
                                 m=m,
                                 c=r,
                                 psi=psi,
                                 model=ifelse(model=="underdamped",1,0))
    out$sd.x <- NA
    out$sd.y <- NA
  } else {
    out <- out %>% dplyr::mutate(x=mu.x,y=mu.y,smaj=NA,smin=NA,eor=NA,sd.x=NA,sd.y=NA)
  }
  return(out)
}

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
