checkPar <- function(par, model, map=NULL, dat=NULL, spatialCovs = NULL){

  if(!is.list(par)) stop("par must be a list.")
  if(!all(names(par) %in% c("beta","sigma","gamma","mu","vel","psi","tau","rho_o"))) {
    stop("names(par) is limited to c('beta','sigma','gamma','mu','vel','psi','tau','rho_o')")
  }

  if(!is.null(map)){
    if(!is.list(map)) stop("map must be a list.")
    if(!all(names(map) %in% c("beta","sigma","gamma","mu","vel","psi","tau","rho_o"))) {
      stop("names(map) is limited to c('beta','sigma','gamma','mu','vel','psi','tau','rho_o')")
    }
  } else map <- list()

  if(model=="underdamped"){
    if(!is.null(par$gamma)){
      if(any(!is.finite(par$gamma)) || any(par$gamma <= 0)){
        stop("gamma should be greater than zero.")
      } else {
        par$log_gamma <- log(par$gamma)
        par$gamma <- NULL
      }
    } else {
      stop("par$gamma is missing, with no default.")
    }

    if(!is.null(map$gamma)){
      if(length(map$gamma) != length(par$log_gamma)) stop("map$gamma should be of length 1.")
      map$log_gamma <- map$gamma
      map$gamma <- NULL
    }
  } else {
    map$log_gamma <- factor(NA)
    map$gamma <- NULL
    par$log_gamma <- 0
  }

  if(!is.null(par$sigma)){
    if(any(!is.finite(par$sigma)) || any(par$sigma <= 0)){
      stop("sigma should be greater than zero.")
    } else {
      par$log_sigma <- log(par$sigma)
      par$sigma <- NULL
    }

    if(!is.null(map$sigma)){
      if(length(map$sigma) != length(par$log_sigma)) stop("map$sigma should be of length 1.")
      map$log_sigma <- map$sigma
      map$sigma <- NULL
    }
  } else {
    stop("par$sigma is missing, with no default.")
  }

  if(!is.null(par$beta)){
    if(any(!is.finite(par$beta))){
      stop("beta must be finite.")
    }
    if(!is.null(map$beta)){
      if(length(map$beta) != length(par$beta)) stop("map$beta should be of length ", length(par$beta), ".")
    }
  } else {
    stop("par$beta is missing, with no default.")
  }

  # Safely check spatialCovs length
  if(!is.null(spatialCovs)){
    if(length(par$beta) != length(spatialCovs)) {
      stop("par$beta is of length ", length(par$beta), ", but spatialCovs is of length ", length(spatialCovs), ". They must be the same length.")
    }
  }

  if(!is.null(par$psi)){
    if(length(par$psi) != 1) stop("par$psi should be of length 1.")
    if(!is.finite(par$psi) || par$psi <= 0){
      stop("par$psi should be greater than zero.")
    }
    par$l_psi <- log(par$psi)
    par$psi <- NULL

    if(!is.null(map$psi)){
      if(length(map$psi) != length(par$l_psi)) stop("map$psi should be of length 1.")
      map$l_psi <- map$psi
      map$psi <- NULL
    }
  } else {
    par$l_psi <- 0
    if(!is.null(map$psi)){
      if(length(map$psi) != length(par$l_psi)) stop("map$psi should be of length 1.")
      map$l_psi <- map$psi
      map$psi <- NULL
    } else map$l_psi <- factor(NA)
  }

  if(!is.null(par$tau)){
    if(length(par$tau) != 2) stop("par$tau should be of length 2.")
    if(any(!is.finite(par$tau)) || any(par$tau <= 0)){
      stop("par$tau should be greater than zero.")
    }
    par$l_tau <- log(par$tau)
    par$tau <- NULL

    if(!is.null(map$tau)){
      if(length(map$tau) != length(par$l_tau)) stop("map$tau should be of length 2.")
      map$l_tau <- map$tau
      map$tau <- NULL
    }
  } else {
    par$l_tau <- c(0,0)
    if(!is.null(map$tau)){
      if(length(map$tau) != length(par$l_tau)) stop("map$tau should be of length 2.")
      map$l_tau <- map$tau
      map$tau <- NULL
    }
    map$l_tau <- factor(rep(NA,2))
  }

  if(!is.null(par$rho_o)){
    if(length(par$rho_o) != 1) stop("par$rho_o should be of length 1.")
    if(!is.finite(par$rho_o) || ((par$rho_o <= -1) | par$rho_o >= 1)){
      stop("par$rho_o must be > -1 and < 1.")
    }
    par$l_rho_o <- log((1+par$rho_o)/(1-par$rho_o))
    par$rho_o <- NULL

    if(!is.null(map$rho_o)){
      if(length(map$rho_o) != length(par$l_rho_o)) stop("map$rho_o should be of length 1.")
      map$l_rho_o <- map$rho_o
      map$rho_o <- NULL
    }
  } else {
    par$l_rho_o <- 0
    if(!is.null(map$rho_o)){
      if(length(map$rho_o) != length(par$l_rho_o)) stop("map$rho_o should be of length 1.")
      map$l_rho_o <- map$rho_o
      map$rho_o <- NULL
    } else map$l_rho_o <- factor(NA)
  }

  if(!is.null(dat)){
    if(!is.null(par$mu)) par$mu <- t(par$mu)
    if(!is.null(par$vel)) par$vel <- t(par$vel)

    if(model=="overdamped"){
      re <- "mu"
      par$vel <- matrix(0, 2, ncol(dat$Y))
      map$vel <- factor(rep(NA, length(dat$Y)))
      map$log_gamma <- factor(NA)
    } else {
      re <- c("mu","vel")
    }

    if(all(is.na(dat$obs_mod))){
      if(model=="overdamped") re <- NULL
      else re <- "vel"
      par$mu <- dat$Y
      map$mu <- factor(rep(NA, length(dat$Y)))
    }

    if(is.null(par$mu)) stop("par$mu is missing, with no default.")
    if(is.null(par$vel)) stop("par$vel is missing, with no default.")
    if(any(dim(par$mu) != dim(dat$Y))) stop("par$mu must have ", ncol(dat$Y), " rows and 2 columns")
    if(any(dim(par$vel) != dim(dat$Y))) stop("par$vel must have ", ncol(dat$Y), " rows and 2 columns")
    if(any(!is.finite(par$mu)) || any(!is.finite(par$vel))) stop("par$mu and/or par$vel must be finite")

    par <- par[c("beta","log_sigma","log_gamma","mu","vel","l_psi","l_tau","l_rho_o")]
    out <- list(par=par, map=map, re=re)
  } else {
    out <- list(par=par, map=map)
  }

  return(out)
}
