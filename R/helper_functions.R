Rcpp::sourceCpp("src/simulateLangevin.cpp")

rasterList <- function (rast) 
{
  lim <- as.vector(extent(rast))
  res <- raster::res(rast)
  xgrid <- seq(lim[1] + res[1]/2, lim[2] - res[1]/2, by = res[1])
  ygrid <- seq(lim[3] + res[2]/2, lim[4] - res[2]/2, by = res[2])
  z <- t(apply(as.matrix(rast), 2, rev))
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
  raster::raster(ud_rast)
}

plotRaster <- function (rast, norm = FALSE, log = FALSE, scale.name = "", light = FALSE) 
{
  covmap <- data.frame(coordinates(rast), val = values(rast))
  if (norm) {
    s <- sum(covmap$val) * xres(rast) * yres(rast)
    covmap$val <- covmap$val/s
  }
  if (log) {
    covmap$val <- log(covmap$val)
  }
  p <- ggplot(covmap, aes_string(x = "x", y = "y")) + geom_raster(aes_string(fill = "val")) + 
    coord_equal()
  if (light) {
    p <- p + scale_fill_viridis(guide = "none") + theme(axis.title = element_blank(), 
                                                        axis.text = element_blank(), axis.ticks = element_blank())
  }
  else {
    p <- p + scale_fill_viridis(name = scale.name)
  }
  return(p)
}

simCov <- function(sca = 200, irange=0.3, sigma2 = 0.1, kappa = 0.5, M = 2048, N = 2048) {
  
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
  
  # Convert to raster and orient correctly
  spatialCov <- raster::flip(
    raster::raster(
      t(grf_fields), 
      xmn = min(grid_list$x), xmx = max(grid_list$x),
      ymn = min(grid_list$y), ymx = max(grid_list$y)
    ),
    direction = 'y' 
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
