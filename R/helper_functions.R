Rcpp::sourceCpp("src/simulateLangevin.cpp")

measurementError <- function(data,M,m,c,psi){
  M <- tmpM <- abs(rnorm(nbAnimals*obsPerAnimal,0,sd=M))
  m <- tmpm <- abs(rnorm(nbAnimals*obsPerAnimal,0,sd=m))
  M[which(tmpM < tmpm)] <- tmpm[which(tmpM < tmpm)]
  m[which(tmpM < tmpm)] <- tmpM[which(tmpM < tmpm)]
  c <- momentuHMM:::radian(runif(nbAnimals*obsPerAnimal,c[1],c[2]))
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
  for(j in which(unlist(lapply(ssm_results,is.null)))){
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

simulate_udLangevin <- function(nbAnimals,obsPerAnimal,lambda,gamma,sigma,beta,spatialCovs,initialPosition=matrix(0,nbAnimals,2),min_dt=4.e-5,ncores=1,progress=FALSE,UD=NULL){
  
  # simulate individuals in parallel
  if(ncores>1) future::plan(future::multisession, workers = ncores)
  
  if(!is.null(UD) & ncores==1) {
    par(mfrow=c(1,1))
    plot(rhabitToRaster(UD))
  }
  bb <- bbox(spatialCovs[[1]])

  s2 <- sigma^2
  
  simDat <- foreach(i = 1:nbAnimals, .combine = "rbind") %dorng% {
    if(progress) message("Individual ",i,'\n')
    iDat <- data.frame(ID=rep(i,obsPerAnimal),time = NA, dt = NA, mu.x = NA, mu.y = NA, v_mux = NA, v_muy = NA)
    for(cov in names(spatialCovs)){
      iDat[[paste0(cov,".x")]] <- NA
      iDat[[paste0(cov,".y")]] <- NA
    }
    waitTimes <- rexp(obsPerAnimal-1,lambda)
    while(any(waitTimes < min_dt)){ # get rid of tiny wait times because they can cause numerical issues
      tInd <- which(waitTimes < min_dt)
      waitTimes[tInd] <- rexp(length(tInd),lambda)
    }
    iDat$time <- cumsum(c(0,waitTimes))
    iDat$dt <- c(0,diff(iDat$time))
    iDat$mu.x[1] <- initialPosition[i,1]
    iDat$mu.y[1] <- initialPosition[i,2]
    iDat$v_mux[1] <- rnorm(1,0,sigma)
    iDat$v_muy[1] <- rnorm(1,0,sigma) 
    iDat[1,] <- momentuHMM:::getGradients(iDat[1,],spatialCovs,coordNames = c("mu.x","mu.y"))
    for(t in 1:(obsPerAnimal-1)){
      if(progress) cat("    Iteration ",t+1,"\r")
      dt_step <- iDat$dt[t+1]
      exp_gdt = exp(-gamma * dt_step);
      exp_2gdt = exp(-2 * gamma * dt_step);
      h = numeric(2)
      h[1] <- sum(iDat[t,paste0(names(spatialCovs),".x")] * beta)
      h[2] <- sum(iDat[t,paste0(names(spatialCovs),".y")] * beta)
      pred_mux <- iDat$mu.x[t] + 
        iDat$v_mux[t]/gamma * (1 - exp_gdt) +
        s2*h[1]/gamma * (dt_step - (1 - exp_gdt)/gamma); 
      pred_muy <- iDat$mu.y[t] + 
        iDat$v_muy[t]/gamma * (1 - exp_gdt) +
        s2*h[2]/gamma * (dt_step - (1 - exp_gdt)/gamma); 
      pred_v_mux = iDat$v_mux[t] * exp_gdt +
        s2*h[1]/gamma * (1 - exp_gdt);
      pred_v_muy = iDat$v_muy[t] * exp_gdt +
        s2*h[2]/gamma * (1 - exp_gdt);
      
      var_x = s2/(gamma*gamma) * 
        (2*gamma*dt_step - 3 + 
         4*exp_gdt - exp_2gdt);
      var_v = s2 * (1 - exp_2gdt);
      cov_xv = s2/gamma * (1 - 2*exp_gdt + exp_2gdt);
      Sigma <- matrix(0,4,4)
      Sigma[1,1] <- Sigma[3,3] <- var_x
      Sigma[1,2] <- Sigma[2,1] <- Sigma[3,4] <- Sigma[4,3] <- cov_xv
      Sigma[2,2] <- Sigma[4,4] <- var_v
  
      mu_v <- mvtnorm::rmvnorm(1,c(pred_mux,pred_v_mux,pred_muy,pred_v_muy),Sigma)
      iDat$mu.x[t+1] <- mu_v[1]
      iDat$mu.y[t+1] <- mu_v[3]
      iDat$v_mux[t+1] <- mu_v[2]
      iDat$v_muy[t+1] <- mu_v[4]
      if(iDat$mu.x[t+1]<bb[1,1] | iDat$mu.x[t+1]>bb[1,2] | iDat$mu.y[t+1]<bb[2,1] | iDat$mu.y[t+1]>bb[2,2]) stop("movement is beyond the extent of the raster(s)")
      iDat[t+1,] <- momentuHMM:::getGradients(iDat[t+1,],spatialCovs,coordNames = c("mu.x","mu.y"))
      if(!is.null(UD) & ncores==1) points(iDat$mu.x[(i-1)*obsPerAnimal+t],iDat$mu.y[(i-1)*obsPerAnimal+t],col=i,type="o",pch=20)
    }
    if(progress) cat('\n')
    return(iDat)
  }
  if(ncores>1) future::plan(future::sequential)
  return(simDat)
}
