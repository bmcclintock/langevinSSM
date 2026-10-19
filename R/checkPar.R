# --- Helper: Validate Observation Error Parameters ---
validate_obs_error <- function(par, map, dat) {
  has_ee <- any(!is.na(dat$smaj))
  has_ls <- any(!is.na(dat$K[,1]))

  if (!has_ee && ("psi" %in% names(par) || "psi" %in% names(map))) {
    stop("Cannot specify 'psi' in 'par' or 'map' because the data does not contain error ellipse observations.")
  }
  if (!has_ee && any(.ell_par_names %in% c(names(par), names(map)))) {
    stop("Cannot specify 'psi_maj', 'b_maj' or 'b_min' in 'par' or 'map' because the data does not contain error ellipse observations.")
  }
  if (!has_ls && (any(c("tau", "rho_o") %in% names(par)) || any(c("tau", "rho_o") %in% names(map)))) {
    stop("Cannot specify 'tau' or 'rho_o' in 'par' or 'map' because the data does not contain standard error observations ('x.err', 'y.err').")
  }
  # dat$obs_dist = any_t + 2 * any_Huber
  if (!dat$obs_dist %in% c(1, 3) && (("df" %in% names(par)) || ("df" %in% names(map)))) {
    stop("Cannot specify 'df' in 'par' or 'map' unless obs_dist = 't' (for at least some observations).")
  }
  if (dat$obs_dist < 2 && (("k" %in% names(par)) || ("k" %in% names(map)))) {
    stop("Cannot specify 'k' in 'par' or 'map' unless obs_dist = 'Huber' (for at least some observations).")
  }
}

# --- Helper: Parse Bayesian Priors ---
# Each row of 'prior' targets elements of the working-scale parameters by
#   * the exact element name ('beta_cov1', 'l_tau_3', 'l_df_2', 'mu.x_1', ...),
#   * a group name of a grouped observation parameter ('l_tau.G', 'l_tau.x.G', 'l_tau.x', 'l_df.B', 'l_k.epar'; see .prior_aliases()),
#   * or the base name ('beta', 'l_tau', 'l_df', ...; all of its elements).
# Priors are only applied to estimated parameters: a row whose elements are all fixed (by 'map', or by default because the
# parameter is not in 'par' or not used by the model) is an error, and a base or group name only applies to the estimated
# elements. Elements tied by 'map' are one estimated parameter and receive the prior once (with a warning if the user gave
# them different priors).
parse_priors <- function(par, prior, spatialCovs, map = NULL, obs_spec = NULL) {
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
    base_of <- rep(names(par), lengths(par))   # working parameter of each element
    aliases <- .prior_aliases(par, obs_spec)
    lev <- .prior_levels(par, map)             # estimated parameter of each element (NA = fixed; "re" = random effects)
    src <- rep(NA_integer_, length(pr_mean))   # row of 'prior' that set each element
    for (i in seq_len(nrow(prior))) {
      rn <- rownames(prior)[i]
      if (rn %in% nice_names) {
        idx <- match(rn, nice_names); kind <- "element"
      } else if (rn %in% names(aliases)) {
        idx <- aliases[[rn]]; kind <- "group"
      } else if (rn %in% names(par)) {
        idx <- which(base_of == rn); kind <- "base"
      } else {
        beta_names <- nice_names[startsWith(nice_names, "beta_")]
        if (length(beta_names) == 0) {
          example_beta_str <- "'beta_1'"
        } else if (length(beta_names) < 3) {
          example_beta_str <- paste0("'", beta_names, "'", collapse = ", ")
        } else {
          example_beta_str <- paste0("'", beta_names[1], "', ..., '", beta_names[length(beta_names)], "'")
        }
        grp_str <- if (length(aliases)) paste0("\n  Groups of observation parameters for this model: ", paste0("'", names(aliases), "'", collapse = ", "), ".") else ""
        stop("Prior parameter name '", rn, "' is invalid. Acceptable base parameters for this model are: ",
             paste(names(par), collapse = ", "),
             "\n  To target specific elements, use their exact working names (e.g., ",
             example_beta_str, ", 'mu.x_1', 'vel.y_10').", grp_str)
      }
      fixed <- is.na(lev[idx])
      if (all(fixed)) stop(.prior_fixed_msg(rn, base_of[idx[1]], kind), call. = FALSE)
      idx <- idx[!fixed]
      pr_mean[idx] <- prior[i, 1]
      pr_sd[idx]   <- prior[i, 2]
      src[idx] <- i
    }

    # elements tied by 'map' are a single estimated parameter: apply its prior once
    has <- !is.na(src) & !is.na(lev) & lev != "re"
    for (L in unique(lev[has])) {
      idx <- which(has & lev == L)
      if (length(idx) < 2) next
      rows <- unique(src[idx])
      if (length(rows) > 1) {
        warning("Elements ", paste0("'", nice_names[idx], "'", collapse = ", "), " are tied by 'map' (one estimated parameter) but ",
                "have priors from different rows of 'prior' (", paste0("'", rownames(prior)[rows], "'", collapse = ", "),
                "). Only the prior of '", nice_names[idx[1]], "' (from row '", rownames(prior)[src[idx[1]]], "') is used.", call. = FALSE)
      }
      pr_mean[idx[-1]] <- NA_real_
      pr_sd[idx[-1]]   <- NA_real_
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

  priors$has_prior_log_lambda <- as.integer(!is.null(pr_mean_list$log_lambda) && any(!is.na(pr_mean_list$log_lambda)))
  priors$prior_mean_log_lambda <- if (!is.null(pr_mean_list$log_lambda)) as.numeric(pr_mean_list$log_lambda) else NA_real_
  priors$prior_sd_log_lambda   <- if (!is.null(pr_sd_list$log_lambda)) as.numeric(pr_sd_list$log_lambda) else NA_real_

  priors$has_prior_l_psi <- as.integer(!is.null(pr_mean_list$l_psi) && any(!is.na(pr_mean_list$l_psi)))
  priors$prior_mean_l_psi <- if (!is.null(pr_mean_list$l_psi)) as.numeric(pr_mean_list$l_psi) else NA_real_
  priors$prior_sd_l_psi   <- if (!is.null(pr_sd_list$l_psi)) as.numeric(pr_sd_list$l_psi) else NA_real_

  for (nm in c("l_psi_maj", "l_b_maj", "l_b_min")) {
    priors[[paste0("has_prior_", nm)]] <- as.integer(!is.null(pr_mean_list[[nm]]) && any(!is.na(pr_mean_list[[nm]])))
    priors[[paste0("prior_mean_", nm)]] <- if (!is.null(pr_mean_list[[nm]])) as.numeric(pr_mean_list[[nm]]) else NA_real_
    priors[[paste0("prior_sd_", nm)]]   <- if (!is.null(pr_sd_list[[nm]])) as.numeric(pr_sd_list[[nm]]) else NA_real_
  }

  priors$has_prior_l_tau <- as.integer(!is.null(pr_mean_list$l_tau) && any(!is.na(pr_mean_list$l_tau)))
  priors$prior_mean_l_tau <- if (!is.null(pr_mean_list$l_tau)) as.numeric(pr_mean_list$l_tau) else rep(NA_real_, 2)
  priors$prior_sd_l_tau   <- if (!is.null(pr_sd_list$l_tau)) as.numeric(pr_sd_list$l_tau) else rep(NA_real_, 2)

  priors$has_prior_l_rho_o <- as.integer(!is.null(pr_mean_list$l_rho_o) && any(!is.na(pr_mean_list$l_rho_o)))
  priors$prior_mean_l_rho_o <- if (!is.null(pr_mean_list$l_rho_o)) as.numeric(pr_mean_list$l_rho_o) else NA_real_
  priors$prior_sd_l_rho_o   <- if (!is.null(pr_sd_list$l_rho_o)) as.numeric(pr_sd_list$l_rho_o) else NA_real_

  priors$has_prior_l_df <- as.integer(!is.null(pr_mean_list$l_df) && any(!is.na(pr_mean_list$l_df)))
  priors$prior_mean_l_df <- if (!is.null(pr_mean_list$l_df)) as.numeric(pr_mean_list$l_df) else NA_real_
  priors$prior_sd_l_df   <- if (!is.null(pr_sd_list$l_df)) as.numeric(pr_sd_list$l_df) else NA_real_

  priors$has_prior_l_k <- as.integer(!is.null(pr_mean_list$l_k) && any(!is.na(pr_mean_list$l_k)))
  priors$prior_mean_l_k <- if (!is.null(pr_mean_list$l_k)) as.numeric(pr_mean_list$l_k) else NA_real_
  priors$prior_sd_l_k   <- if (!is.null(pr_sd_list$l_k)) as.numeric(pr_sd_list$l_k) else NA_real_

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
# lambda_mode: 0 = fixed or no barrier, 1 = tied, 2 = estimate (par$lambda/map$lambda only allowed for 2)
checkPar <- function(par, model, map=NULL, dat=NULL, spatialCovs = NULL, prior = NULL, lambda_mode = 0L, obs_spec = NULL){

  if(!is.list(par)) stop("par must be a list.")
  if(!all(names(par) %in% .par_names_all)) {
    stop("names(par) is limited to c(", paste0("'", .par_names_all, "'", collapse = ","), ")")
  }
  if(!is.null(map)){
    if(!is.list(map)) stop("map must be a list.")
    if(!all(names(map) %in% .par_names_all)) {
      stop("names(map) is limited to c(", paste0("'", .par_names_all, "'", collapse = ","), ")")
    }
  } else map <- list()

  # barrier strength (natural scale); a free parameter only when lambda = "estimate"
  if (lambda_mode != 2L && ("lambda" %in% names(par) || "lambda" %in% names(map))) {
    stop("'lambda' can only be specified in 'par' or 'map' when a barrier is specified and lambda = \"estimate\".")
  }
  if (!is.null(par$lambda)) {
    if (length(par$lambda) != 1 || !is.numeric(par$lambda) || !is.finite(par$lambda) || par$lambda <= 0) {
      stop("par$lambda should be a single value greater than zero.")
    }
    par$log_lambda <- log(par$lambda)
    par$lambda <- NULL
    if (!is.null(map$lambda)) {
      if (length(map$lambda) != 1) stop("map$lambda should be of length 1.")
      map$log_lambda <- map$lambda
      map$lambda <- NULL
    }
  } else {
    if (lambda_mode == 2L) stop("par$lambda is missing, with no default.")
    par$log_lambda <- 0
    map$log_lambda <- factor(NA)
  }

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

  # error ellipse power law (scale of the semi-major axis and exponents of both axes); default 1, fixed
  for (nm in .ell_par_names) {
    wn <- paste0("l_", nm)
    if (!is.null(par[[nm]])) {
      if (length(par[[nm]]) != 1) stop("par$", nm, " should be of length 1.")
      if (!is.finite(par[[nm]]) || par[[nm]] <= 0) stop("par$", nm, " should be greater than zero.")
      par[[wn]] <- log(par[[nm]])
      par[[nm]] <- NULL
    } else {
      par[[wn]] <- 0
      if (is.null(map[[nm]])) map[[wn]] <- factor(NA)
    }
    if (!is.null(map[[nm]])) {
      if (length(map[[nm]]) != 1) stop("map$", nm, " should be of length 1.")
      map[[wn]] <- map[[nm]]
      map[[nm]] <- NULL
    }
  }

  # tau: a 2-vector (x, y) shared by all LS/GPS observations, or a 2-row matrix (or named list) with one
  # column per location-class group; working parameter l_tau is the column-major vector
  if (!is.null(map$tau)) map$tau <- .norm_tau_map(map$tau)
  if(!is.null(par$tau)){
    par$tau <- .norm_tau(par$tau)
    if(any(!is.finite(par$tau)) || any(par$tau <= 0)){
      stop("par$tau should be greater than zero.")
    }
    par$l_tau <- as.vector(log(par$tau))
    par$tau <- NULL

    if(!is.null(map$tau)){
      if(length(map$tau) != length(par$l_tau)) stop("map$tau should be of length ", length(par$l_tau), ".")
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

  # df: one value shared by all observations, or a named vector with one value per group
  # (error types and/or location classes; see .obs_spec())
  if(!is.null(par$df)){
    if (is.list(par$df)) par$df <- unlist(par$df)
    if(any(!is.finite(par$df)) || any(par$df <= 0)){
      stop("par$df should be greater than zero.")
    }
    par$l_df <- unname(log(par$df))
    par$df <- NULL

    if(!is.null(map$df)){
      if(length(map$df) != length(par$l_df)) stop("map$df should be of length ", length(par$l_df), ".")
      map$l_df <- map$df
      map$df <- NULL
    }
  } else {
    par$l_df <- log(5)
    if(!is.null(map$df)){
      if(length(map$df) != length(par$l_df)) stop("map$df should be of length 1.")
      map$l_df <- map$df
      map$df <- NULL
    } else map$l_df <- factor(NA)
  }

  # k: one value shared by all observations, or a named vector with one value per group
  # (error types and/or location classes; see .obs_spec())
  if(!is.null(par$k)){
    if (is.list(par$k)) par$k <- unlist(par$k)
    if(any(!is.finite(par$k)) || any(par$k <= 0)){
      stop("par$k should be greater than zero.")
    }
    par$l_k <- unname(log(par$k))
    par$k <- NULL

    if(!is.null(map$k)){
      if(length(map$k) != length(par$l_k)) stop("map$k should be of length ", length(par$l_k), ".")
      map$l_k <- map$k
      map$k <- NULL
    }
  } else {
    par$l_k <- log(2)
    if(!is.null(map$k)){
      if(length(map$k) != length(par$l_k)) stop("map$k should be of length 1.")
      map$l_k <- map$k
      map$k <- NULL
    } else map$l_k <- factor(NA)
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

    par <- par[.tmb_par_order]
  }

  priors <- parse_priors(par, prior, spatialCovs, map = map, obs_spec = obs_spec)

  if(!is.null(dat)){
    out <- list(par=par, map=map, re=re, priors=priors)
  } else {
    out <- list(par=par, map=map, priors=priors)
  }

  return(out)
}

# Group names of grouped observation parameters that can be used in 'prior' (global element indices of unlist(par)):
#   'l_tau.x', 'l_tau.y'                 x (y) element of tau of every group
#   'l_tau.<g>', 'l_tau.x.<g>', 'l_tau.y.<g>'   tau of group <g> (x and y, x, y), e.g. 'l_tau.G', 'l_tau.x.3'
#   'l_df.<g>', 'l_k.<g>'                df (k) of group <g>, e.g. 'l_df.B', 'l_k.epar'
# '.' separates group names because positional element names ('l_tau_1', 'l_df_2') use '_' and classes include "1", "2", ...
.prior_aliases <- function(par, obs_spec = NULL) {
  off <- cumsum(c(0, lengths(par)))[seq_along(par)]
  names(off) <- names(par)
  al <- list()
  if (!is.null(par$l_tau) && length(par$l_tau) %% 2 == 0) {
    n <- length(par$l_tau)
    xi <- unname(off["l_tau"]) + seq(1, n, by = 2)
    yi <- xi + 1
    al[["l_tau.x"]] <- xi
    al[["l_tau.y"]] <- yi
    labs <- obs_spec$tau_labels
    if (!is.null(labs) && length(labs) == n / 2) {
      for (j in seq_along(labs)) {
        al[[paste0("l_tau.", labs[j])]] <- c(xi[j], yi[j])
        al[[paste0("l_tau.x.", labs[j])]] <- xi[j]
        al[[paste0("l_tau.y.", labs[j])]] <- yi[j]
      }
    }
  }
  for (nm in c("df", "k")) {
    wn <- paste0("l_", nm)
    labs <- obs_spec[[paste0(nm, "_labels")]]
    if (!is.null(par[[wn]]) && !is.null(labs) && length(labs) == length(par[[wn]])) {
      for (j in seq_along(labs)) al[[paste0(wn, ".", labs[j])]] <- unname(off[wn]) + j
    }
  }
  al
}

# Estimated parameter of each element of unlist(par): "<name>:<map level>" (elements with the same value are tied),
# NA if fixed by 'map', and "re" for the random effects (mu, vel), which are not checked.
.prior_levels <- function(par, map = NULL) {
  unlist(lapply(names(par), function(nm) {
    n <- length(par[[nm]])
    if (nm %in% c("mu", "vel")) return(rep("re", n))
    m <- map[[nm]]
    if (is.null(m)) return(paste0(nm, ":", seq_len(n)))
    m <- as.character(m)
    ifelse(is.na(m), NA_character_, paste0(nm, ":", m))
  }), use.names = FALSE)
}

.prior_fixed_msg <- function(rn, wn, kind) {
  what <- if (kind == "base" || identical(rn, wn)) paste0("'", wn, "'") else paste0("'", rn, "' (of '", wn, "')")
  nat <- c(log_sigma = "sigma", log_gamma = "gamma", log_lambda = "lambda", l_psi = "psi", l_psi_maj = "psi_maj", l_b_maj = "b_maj",
           l_b_min = "b_min", l_tau = "tau", l_rho_o = "rho_o", l_df = "df", l_k = "k", beta = "beta")[wn]
  hint <- switch(wn,
    l_df = " 'df' is only estimated if it is included in 'par' and some observations have obs_dist = \"t\".",
    l_k = " 'k' is only estimated if it is included in 'par' and some observations have obs_dist = \"Huber\".",
    l_psi_maj = , l_b_maj = , l_b_min = " The error ellipse power law parameters are only estimated if they are included in 'par'.",
    l_psi = , l_tau = , l_rho_o = paste0(" '", nat, "' is only estimated if it is included in 'par'."),
    log_gamma = " 'gamma' is not a parameter of the overdamped model.",
    "")
  paste0("A prior is specified for ", what, ", which is not estimated (fixed by 'map', or by default).", hint,
         " Priors can only be specified for estimated parameters.")
}

# Names of the natural-scale parameters accepted in 'par' and 'map'
.par_names_all <- c("beta","sigma","gamma","lambda","mu","vel","psi","psi_maj","b_maj","b_min","tau","rho_o","df","k")

# Order of the working parameters of the TMB template
.tmb_par_order <- c("beta","log_sigma","log_gamma","log_lambda","mu","vel","l_psi","l_psi_maj","l_b_maj","l_b_min",
                    "l_tau","l_rho_o","l_df","l_k")

# Fits saved by earlier versions lack parameters that were added to the template later (l_k for
# obs_dist = "Huber"; l_psi_maj, l_b_maj, l_b_min for the error ellipse power law) and their prior
# entries. Add them at their neutral values, fixed via map, so that the fit can be rebuilt with
# the current template.
.upgrade_tmb_setup <- function(ts) {
  if (is.null(ts) || is.null(ts$parList) || is.null(ts$parList$mu)) return(ts)
  defaults <- list(l_psi_maj = 0, l_b_maj = 0, l_b_min = 0, l_k = log(2))
  missing_par <- setdiff(names(defaults), names(ts$parList))
  if (!length(missing_par)) return(ts)
  if (is.null(ts$map)) ts$map <- list()
  for (nm in missing_par) {
    ts$parList[[nm]] <- defaults[[nm]]
    ts$map[[nm]] <- factor(NA)
    if (!is.null(ts$priors) && is.null(ts$priors[[paste0("has_prior_", nm)]])) {
      ts$priors[[paste0("has_prior_", nm)]] <- 0L
      ts$priors[[paste0("prior_mean_", nm)]] <- NA_real_
      ts$priors[[paste0("prior_sd_", nm)]] <- NA_real_
    }
  }
  ts$parList <- ts$parList[c(intersect(.tmb_par_order, names(ts$parList)), setdiff(names(ts$parList), .tmb_par_order))]
  ts
}
