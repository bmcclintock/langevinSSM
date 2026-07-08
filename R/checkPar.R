# --- Helper: Validate Observation Error Parameters ---
validate_obs_error <- function(par, map, dat) {
  has_ee <- any(!is.na(dat$smaj))
  has_ls <- any(!is.na(dat$K[,1]))

  if (!has_ee && ("psi" %in% names(par) || "psi" %in% names(map))) {
    stop("Cannot specify 'psi' in 'par' or 'map' because the data does not contain error ellipse observations.")
  }
  if (!has_ls && (any(c("tau", "rho_o") %in% names(par)) || any(c("tau", "rho_o") %in% names(map)))) {
    stop("Cannot specify 'tau' or 'rho_o' in 'par' or 'map' because the data does not contain standard error observations ('x.err', 'y.err').")
  }
}

# --- Helper: Parse Bayesian Priors ---
parse_priors <- function(par, prior, spatialCovs) {
  pr_mean <- unlist(par) * 0 + NA_real_
  pr_sd   <- unlist(par) * 0 + NA_real_

  nice_names <- unlist(lapply(names(par), function(nm) {
    n_elem <- length(par[[nm]])
    if (n_elem == 1) return(nm)
    if (nm == "beta" && !is.null(names(spatialCovs)) && length(spatialCovs) == n_elem) {
      return(paste0("beta_", names(spatialCovs)))
    }
    if (nm %in% c("mu", "vel")) {
      num_steps <- n_elem / 2
      return(paste0(nm, c(".x_", ".y_"), rep(seq_len(num_steps), each = 2)))
    }
    return(paste0(nm, "_", seq_len(n_elem)))
  }))

  names(pr_mean) <- nice_names
  names(pr_sd) <- nice_names

  if (!is.null(prior)) {
    if (!is.data.frame(prior) || ncol(prior) != 2) {
      stop("'prior' must be a 2-column data frame with mean in col 1 and sd in col 2.")
    }
    for (rn in rownames(prior)) {
      if (rn %in% nice_names) {
        pr_mean[rn] <- prior[rn, 1]
        pr_sd[rn]   <- prior[rn, 2]
      } else if (rn %in% names(par)) {
        idx <- which(startsWith(nice_names, paste0(rn, "_")) | nice_names == rn)
        pr_mean[idx] <- prior[rn, 1]
        pr_sd[idx]   <- prior[rn, 2]
      } else {
        beta_names <- nice_names[startsWith(nice_names, "beta_")]
        if (length(beta_names) == 0) {
          example_beta_str <- "'beta_1'"
        } else if (length(beta_names) < 3) {
          example_beta_str <- paste0("'", beta_names, "'", collapse = ", ")
        } else {
          example_beta_str <- paste0("'", beta_names[1], "', ..., '", beta_names[length(beta_names)], "'")
        }
        stop("Prior parameter name '", rn, "' is invalid. Acceptable base parameters for this model are: ",
             paste(names(par), collapse = ", "),
             "\n  To target specific elements, use their exact working names (e.g., ",
             example_beta_str, ", 'mu.x_1', 'vel.y_10').")
      }
    }
  }

  pr_mean_list <- utils::relist(unname(pr_mean), skeleton = par)
  pr_sd_list   <- utils::relist(unname(pr_sd), skeleton = par)

  priors <- list()

  priors$has_prior_beta <- as.integer(any(!is.na(pr_mean_list$beta)))
  priors$prior_mean_beta <- as.numeric(pr_mean_list$beta)
  priors$prior_sd_beta   <- as.numeric(pr_sd_list$beta)

  priors$has_prior_log_sigma <- as.integer(any(!is.na(pr_mean_list$log_sigma)))
  priors$prior_mean_log_sigma <- as.numeric(pr_mean_list$log_sigma)
  priors$prior_sd_log_sigma   <- as.numeric(pr_sd_list$log_sigma)

  priors$has_prior_log_gamma <- as.integer(!is.null(pr_mean_list$log_gamma) && any(!is.na(pr_mean_list$log_gamma)))
  priors$prior_mean_log_gamma <- if (!is.null(pr_mean_list$log_gamma)) as.numeric(pr_mean_list$log_gamma) else NA_real_
  priors$prior_sd_log_gamma   <- if (!is.null(pr_sd_list$log_gamma)) as.numeric(pr_sd_list$log_gamma) else NA_real_

  priors$has_prior_l_psi <- as.integer(!is.null(pr_mean_list$l_psi) && any(!is.na(pr_mean_list$l_psi)))
  priors$prior_mean_l_psi <- if (!is.null(pr_mean_list$l_psi)) as.numeric(pr_mean_list$l_psi) else NA_real_
  priors$prior_sd_l_psi   <- if (!is.null(pr_sd_list$l_psi)) as.numeric(pr_sd_list$l_psi) else NA_real_

  priors$has_prior_l_tau <- as.integer(!is.null(pr_mean_list$l_tau) && any(!is.na(pr_mean_list$l_tau)))
  priors$prior_mean_l_tau <- if (!is.null(pr_mean_list$l_tau)) as.numeric(pr_mean_list$l_tau) else rep(NA_real_, 2)
  priors$prior_sd_l_tau   <- if (!is.null(pr_sd_list$l_tau)) as.numeric(pr_sd_list$l_tau) else rep(NA_real_, 2)

  priors$has_prior_l_rho_o <- as.integer(!is.null(pr_mean_list$l_rho_o) && any(!is.na(pr_mean_list$l_rho_o)))
  priors$prior_mean_l_rho_o <- if (!is.null(pr_mean_list$l_rho_o)) as.numeric(pr_mean_list$l_rho_o) else NA_real_
  priors$prior_sd_l_rho_o   <- if (!is.null(pr_sd_list$l_rho_o)) as.numeric(pr_sd_list$l_rho_o) else NA_real_

  priors$has_prior_mu <- as.integer(!is.null(pr_mean_list$mu) && any(!is.na(pr_mean_list$mu)))
  if (priors$has_prior_mu == 1L) {
    valid_mu <- which(!is.na(pr_mean_list$mu))
    priors$prior_idx_mu <- as.integer(valid_mu - 1L) # 0-based for C++
    priors$prior_mean_mu_val <- as.numeric(pr_mean_list$mu[valid_mu])
    priors$prior_sd_mu_val   <- as.numeric(pr_sd_list$mu[valid_mu])
  } else {
    priors$prior_idx_mu <- integer(0)
    priors$prior_mean_mu_val <- numeric(0)
    priors$prior_sd_mu_val   <- numeric(0)
  }

  priors$has_prior_vel <- as.integer(!is.null(pr_mean_list$vel) && any(!is.na(pr_mean_list$vel)))
  if (priors$has_prior_vel == 1L) {
    valid_vel <- which(!is.na(pr_mean_list$vel))
    priors$prior_idx_vel <- as.integer(valid_vel - 1L) # 0-based for C++
    priors$prior_mean_vel_val <- as.numeric(pr_mean_list$vel[valid_vel])
    priors$prior_sd_vel_val   <- as.numeric(pr_sd_list$vel[valid_vel])
  } else {
    priors$prior_idx_vel <- integer(0)
    priors$prior_mean_vel_val <- numeric(0)
    priors$prior_sd_vel_val   <- numeric(0)
  }

  return(priors)
}

# --- Main checkPar Function ---
checkPar <- function(par, model, map=NULL, dat=NULL, spatialCovs = NULL, prior = NULL){

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

  if (model == "overdamped") {
    if (any(c("gamma", "vel") %in% names(par)) || any(c("gamma", "vel") %in% names(map))) {
      stop("Cannot specify 'gamma' or 'vel' in 'par' or 'map' when model = 'overdamped'.")
    }
  }

  if(!is.null(dat)) {
    validate_obs_error(par, map, dat)
  }

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

    if (all(is.na(dat$obs_mod))) {
      has_nas <- any(is.na(dat$Y))

      if (model == "overdamped") {
        re <- if (has_nas) "mu" else NULL
      } else {
        re <- if (has_nas) c("mu", "vel") else "vel"
      }

      valid_Y <- !is.na(dat$Y)
      par$mu[valid_Y] <- dat$Y[valid_Y] * dat$scale_factor

      if (has_nas) {
        mu_map <- rep(NA, length(dat$Y))
        na_idx <- which(is.na(dat$Y))
        mu_map[na_idx] <- 1:length(na_idx)
        map$mu <- factor(mu_map)
      } else {
        map$mu <- factor(rep(NA, length(dat$Y)))
      }
    }

    if(is.null(par$mu)) stop("par$mu is missing, with no default.")
    if(is.null(par$vel)) stop("par$vel is missing, with no default.")
    if(any(dim(par$mu) != dim(dat$Y))) stop("par$mu must have ", ncol(dat$Y), " rows and 2 columns")
    if(any(dim(par$vel) != dim(dat$Y))) stop("par$vel must have ", ncol(dat$Y), " rows and 2 columns")
    if(any(!is.finite(par$mu)) || any(!is.finite(par$vel))) stop("par$mu and/or par$vel must be finite")

    par <- par[c("beta","log_sigma","log_gamma","mu","vel","l_psi","l_tau","l_rho_o")]
  }

  priors <- parse_priors(par, prior, spatialCovs)

  if(!is.null(dat)){
    out <- list(par=par, map=map, re=re, priors=priors)
  } else {
    out <- list(par=par, map=map, priors=priors)
  }

  return(out)
}
