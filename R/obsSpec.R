# Observation-model structure: the observation distribution, tau, df and k by error type or location class,
# and the error ellipse power law. Internal helpers shared by fitLangevin(), build_tmb_data(),
# addMeasurementError() and the downstream functions.
#
# Groups are named with the same terms as formatData():
#   "epar"   observations with an error ellipse (smaj, smin, eor)
#   "sderr"  observations with x- and y-axis standard errors (x.err, y.err)
#   "3", "2", "1", "0", "A", "B", "Z", "G", "GL"   observations of that location quality class ('lc')
# An observation takes the entry of its location class if there is one, and otherwise the entry of its error type.

.lc_all <- c("3", "2", "1", "0", "A", "B", "Z", "G", "GL")
.err_types <- c("epar", "sderr")
.dist_choices <- c("Gaussian", "t", "Huber")
.ell_par_names <- c("psi_maj", "b_maj", "b_min")

# Error type of each row: "epar" (error ellipse), "sderr" (x.err/y.err) or NA (no measurement error / missing coordinates)
.obs_type <- function(data, coord = c("x", "y")) {
  ok <- !is.na(data[[coord[1]]]) & !is.na(data[[coord[2]]])
  has <- function(v) if (v %in% names(data)) !is.na(data[[v]]) else rep(FALSE, nrow(data))
  ty <- rep(NA_character_, nrow(data))
  ty[ok & has("x.err") & has("y.err")] <- "sderr"
  ty[ok & has("smaj") & has("smin") & has("eor")] <- "epar"
  ty
}

.data_lc <- function(data) if ("lc" %in% names(data)) as.character(data$lc) else rep(NA_character_, nrow(data))

# Group index (1-based; NA = unmatched) of each observation: the entry of its location class, otherwise of its error type
.group_index <- function(labels, lc, type) {
  idx <- match(as.character(lc), labels)
  ti <- match(type, labels)
  # fits from development version r29 may have an "Argos" group (all Argos location classes)
  if ("Argos" %in% labels) {
    ar <- is.na(idx) & !is.na(lc) & !lc %in% c("G", "GL")
    idx[ar] <- match("Argos", labels)
  }
  idx[is.na(idx)] <- ti[is.na(idx)]
  idx
}

.check_labels <- function(labels, what, types = .err_types) {
  ok <- c(types, .lc_all)
  allowed <- paste0(if (length(types)) paste0("error type (", paste0("'", types, "'", collapse = ", "), ") and/or ") else "",
                    "location class (", paste(.lc_all, collapse = ", "), ")")
  if (is.null(labels) || any(is.na(labels)) || any(labels == "")) {
    stop(what, " must be named by ", allowed, ".", call. = FALSE)
  }
  bad <- setdiff(labels, ok)
  if (length(bad)) stop("Invalid name(s) for ", what, ": ", paste(bad, collapse = ", "), ". Names must be ", allowed, ".", call. = FALSE)
  if (anyDuplicated(labels)) stop("Duplicated name(s) for ", what, ".", call. = FALSE)
  invisible(TRUE)
}

# Assign each observation in 'rows' to a group; error if any cannot be assigned. A group without
# observations is harmless if its parameter is fixed, but is not identifiable if estimated ('free'), so warn.
.assign_groups <- function(labels, data, rows, type, what, kind, free = rep(TRUE, length(labels)), types = .err_types) {
  lc <- .data_lc(data)
  cls <- intersect(labels, .lc_all)
  if (length(cls) && any(rows) && all(is.na(lc[rows]))) {
    stop(what, " has entries for location classes (", paste(cls, collapse = ", "), "), but the data have no location classes ('lc') for the ",
         kind, ". Use ", paste0("'", types, "'", collapse = " and/or "), " instead.", call. = FALSE)
  }
  g <- .group_index(labels, lc, type)
  miss <- rows & is.na(g)
  if (any(miss)) {
    desc <- unique(ifelse(is.na(lc[miss]), paste0("'", type[miss], "' observations without a location class"),
                          paste0("'", type[miss], "' observations of location class ", lc[miss])))
    stop(what, " has no entry for ", paste(desc, collapse = "; "), ". Add an entry for the error type (",
         paste0("'", types, "'", collapse = " or "), ") or for the location class.", call. = FALSE)
  }
  empty <- !seq_along(labels) %in% unique(g[rows]) & free
  if (any(empty)) warning(what, " has estimated entries for ", paste(labels[empty], collapse = ", "),
                          ", but none of the ", kind, " belong to ", if (sum(empty) > 1) "them" else "it",
                          "; ", if (sum(empty) > 1) "they are" else "it is", " not identifiable.", call. = FALSE)
  g
}

# Which groups of a parameter are estimated, given its natural-scale 'map' entry ('per' map entries per group)
.free_groups <- function(map_el, n_groups, per = 1) {
  if (is.null(map_el)) return(rep(TRUE, n_groups))
  m <- matrix(!is.na(map_el), nrow = per)
  as.vector(colSums(m) > 0)
}

# obs_dist: a single distribution (default "Gaussian"), or a named character vector (or list) by error type and/or location class
.norm_obs_dist <- function(obs_dist) {
  if (is.null(obs_dist)) return("Gaussian")
  if (is.list(obs_dist)) {
    if (!all(lengths(obs_dist) == 1)) stop("Each element of 'obs_dist' must be a single distribution.", call. = FALSE)
    obs_dist <- unlist(obs_dist)
  }
  if (!is.character(obs_dist) || !length(obs_dist)) {
    stop("'obs_dist' must be \"Gaussian\", \"t\" or \"Huber\", or a named vector of these.", call. = FALSE)
  }
  if (is.null(names(obs_dist))) {
    if (identical(obs_dist, .dist_choices)) return("Gaussian")   # default argument
    if (length(obs_dist) != 1) {
      stop("'obs_dist' must be a single distribution, or a named vector with names by error type ('epar', 'sderr') and/or location class, ",
           "e.g. obs_dist = c(epar = \"Huber\", sderr = \"Gaussian\").", call. = FALSE)
    }
    return(match.arg(obs_dist, .dist_choices))
  }
  .check_labels(names(obs_dist), "obs_dist")
  bad <- !obs_dist %in% .dist_choices
  if (any(bad)) stop("Invalid distribution(s) in 'obs_dist': ", paste(unique(obs_dist[bad]), collapse = ", "),
                     ". Must be \"Gaussian\", \"t\" or \"Huber\".", call. = FALSE)
  obs_dist
}

.fmt_obs_dist <- function(obs_dist) {
  if (is.null(obs_dist)) return("Gaussian")
  if (is.null(names(obs_dist))) return(paste(obs_dist, collapse = ", "))
  paste0(names(obs_dist), " = ", obs_dist, collapse = ", ")
}

# Does any observation (possibly) use a non-Gaussian distribution?
.non_gaussian <- function(obs_dist) !is.null(obs_dist) && any(obs_dist %in% c("t", "Huber"))

# Distribution code of each row (0 = Gaussian, 1 = t, 2 = Huber). With a single distribution every row gets it
# (as before); with groups, rows without measurement error get 0.
.dist_per_obs <- function(obs_dist, data, coord = c("x", "y")) {
  obs_dist <- .norm_obs_dist(obs_dist)
  n <- nrow(data)
  if (is.null(names(obs_dist))) return(rep(match(obs_dist, .dist_choices) - 1L, n))
  type <- .obs_type(data, coord)
  g <- .group_index(names(obs_dist), .data_lc(data), type)
  rows <- !is.na(type)
  if (any(rows & is.na(g))) .assign_groups(names(obs_dist), data, rows, type, "obs_dist", "observations with measurement error")
  code <- match(unname(obs_dist), .dist_choices)[g] - 1L
  code[!rows | is.na(code)] <- 0L
  as.integer(code)
}

# Bring 'tau' to its canonical form: a numeric vector of length 2 (shared) or a 2-row matrix with
# one column per group (rows x, y). A named list of scalars or 2-vectors is converted.
.norm_tau <- function(tau) {
  if (is.null(tau)) return(NULL)
  if (is.list(tau)) {
    if (is.null(names(tau))) stop("par$tau given as a list must be named by location class (and/or 'sderr').", call. = FALSE)
    if (!all(lengths(tau) %in% c(1, 2))) stop("Each element of par$tau must be of length 1 or 2 (x, y).", call. = FALSE)
    nm <- names(tau)
    tau <- matrix(unlist(lapply(tau, function(v) rep(as.numeric(v), length.out = 2))), nrow = 2)
    colnames(tau) <- nm
  }
  if (is.matrix(tau)) {
    if (nrow(tau) != 2) stop("par$tau given as a matrix must have 2 rows (x, y) and one column per location class.", call. = FALSE)
    rownames(tau) <- c("x", "y")
  } else if (length(tau) != 2) {
    stop("par$tau should be of length 2, or a 2-row matrix (or named list) by location class.", call. = FALSE)
  }
  tau
}

# Natural-scale map entry of a grouped tau: a named list (one entry of length 1 or 2 per class) is
# flattened to the column-major order of the tau matrix.
.norm_tau_map <- function(m) {
  if (is.null(m) || !is.list(m)) return(m)
  factor(unlist(lapply(m, function(v) rep(as.character(v), length.out = 2)), use.names = FALSE))
}

# Reference lengths (natural units) of the error ellipse power law: geometric means of the reported axes.
.ell_ref <- function(data, rows) {
  smaj <- data$smaj[rows]; smin <- data$smin[rows]
  if (!length(smaj)) return(c(maj = 1, min = 1))
  if (any(smaj <= 0) || any(smin <= 0)) {
    stop("The error ellipse power law ('psi_maj', 'b_maj', 'b_min') requires all 'smaj' and 'smin' to be > 0.", call. = FALSE)
  }
  c(maj = exp(mean(log(smaj))), min = exp(mean(log(smin))))
}

.default_obs_spec <- function() list(dist_labels = NULL, tau_labels = NULL, df_labels = NULL, k_labels = NULL,
                                     ell_power = FALSE, ell_ref = c(maj = 1, min = 1))

# Observation-model specification of a fit, derived from obs_dist and the natural-scale 'par' and 'map'.
#   dist_labels  groups of obs_dist (NULL = one distribution for all observations)
#   tau_labels   groups of tau (NULL = one tau shared by all 'sderr' observations)
#   df_labels    groups of df over the observations with obs_dist "t" (NULL = shared)
#   k_labels     groups of k over the observations with obs_dist "Huber" (NULL = shared)
#   ell_power    use the error ellipse power law (psi_maj, b_maj, b_min)
#   ell_ref      reference axis lengths (natural units) of the power law
.obs_spec <- function(par, map, data, obs_dist, coord = c("x", "y")) {
  if (is.null(par)) par <- list()
  if (is.null(map)) map <- list()
  spec <- .default_obs_spec()
  type <- .obs_type(data, coord)
  me <- !is.na(type)
  obs_dist <- .norm_obs_dist(obs_dist)

  if (!is.null(names(obs_dist))) {
    if (!any(me)) stop("'obs_dist' is specified by error type or location class, but the data contain no observations with measurement error.", call. = FALSE)
    .assign_groups(names(obs_dist), data, me, type, "obs_dist", "observations with measurement error", free = rep(FALSE, length(obs_dist)))
    spec$dist_labels <- names(obs_dist)
  }
  code <- .dist_per_obs(obs_dist, data, coord)

  tau <- .norm_tau(par$tau)
  if (is.matrix(tau)) {
    tl <- colnames(tau)
    .check_labels(tl, "par$tau", types = "sderr")
    mt <- .norm_tau_map(map$tau)
    if (!is.null(mt) && length(mt) != length(tau)) stop("map$tau should be of length ", length(tau), ".", call. = FALSE)
    .assign_groups(tl, data, type %in% "sderr", type, "par$tau", "'sderr' (x.err, y.err) observations",
                   free = .free_groups(mt, ncol(tau), 2), types = "sderr")
    spec$tau_labels <- tl
  }

  for (nm in c("df", "k")) {
    d <- if (nm == "df") "t" else "Huber"
    x <- par[[nm]]
    if (is.list(x)) x <- unlist(x)
    if (is.null(x) || (is.null(names(x)) && length(x) == 1)) next
    labels <- names(x)
    .check_labels(labels, paste0("par$", nm))
    mx <- map[[nm]]
    if (!is.null(mx) && length(mx) != length(x)) stop("map$", nm, " should be of length ", length(x), ".", call. = FALSE)
    rows <- me & code == match(d, .dist_choices) - 1L
    if (!any(rows)) stop("par$", nm, " is specified by group, but no observations use obs_dist = \"", d, "\".", call. = FALSE)
    .assign_groups(labels, data, rows, type, paste0("par$", nm), paste0("observations with obs_dist = \"", d, "\""),
                   free = .free_groups(mx, length(x), 1))
    spec[[paste0(nm, "_labels")]] <- labels
  }

  if (any(.ell_par_names %in% c(names(par), names(map)))) {
    if (!any(type %in% "epar")) {
      stop("Cannot specify 'psi_maj', 'b_maj' or 'b_min' in 'par' or 'map' because the data do not contain error ellipse observations.", call. = FALSE)
    }
    spec$ell_power <- TRUE
    ref <- attr(data, "ell_ref")   # set internally by predLangevin()
    spec$ell_ref <- if (!is.null(ref)) ref else .ell_ref(data, type %in% "epar")
  }
  spec
}

# Observation-model specification stored in a fit (fits made before it existed use the defaults)
.cond_obs_spec <- function(cond) {
  spec <- cond$obs_spec
  if (is.null(spec)) return(.default_obs_spec())
  # development versions stored one set of tail labels and used the names "KF"/"LS"
  if (!is.null(spec$tail_labels)) {
    od <- .norm_obs_dist(cond$obs_dist)
    if ("t" %in% od) spec$df_labels <- spec$tail_labels
    if ("Huber" %in% od) spec$k_labels <- spec$tail_labels
    spec$tail_labels <- NULL
  }
  for (nm in c("df_labels", "k_labels")) if (!is.null(spec[[nm]])) {
    l <- spec[[nm]]; l[l == "KF"] <- "epar"; l[l == "LS"] <- "sderr"; spec[[nm]] <- l
  }
  utils::modifyList(.default_obs_spec(), spec, keep.null = TRUE)
}

# TMB data entries of the observation-model specification (per-observation codes and group indices, 0-based)
.obs_spec_data <- function(spec, data, coord, scaleFactor, obs_dist = "Gaussian") {
  if (is.null(spec)) spec <- .default_obs_spec()
  n <- nrow(data)
  type <- .obs_type(data, coord)
  me <- !is.na(type)
  lc <- .data_lc(data)
  code <- .dist_per_obs(obs_dist, data, coord)
  grp <- function(labels, rows, what) {
    if (is.null(labels)) return(integer(n))
    g <- .group_index(labels, lc, type)
    if (any(rows & is.na(g))) stop("The data contain observations that do not belong to any ", what,
                                   " group of the model (", paste(labels, collapse = ", "), ").", call. = FALSE)
    g[is.na(g) | !rows] <- 1L
    as.integer(g - 1L)
  }
  if (isTRUE(spec$ell_power) && any(type %in% "epar" & (data$smaj <= 0 | data$smin <= 0))) {
    stop("The error ellipse power law ('psi_maj', 'b_maj', 'b_min') requires all 'smaj' and 'smin' to be > 0.", call. = FALSE)
  }
  list(
    dist_obs = code,
    tau_group = grp(spec$tau_labels, type %in% "sderr", "tau"),
    df_group = grp(spec$df_labels, me & code == 1L, "df"),
    k_group = grp(spec$k_labels, me & code == 2L, "k"),
    ell_power = as.integer(isTRUE(spec$ell_power)),
    ell_ref_maj = unname(spec$ell_ref[["maj"]]) / scaleFactor,
    ell_ref_min = unname(spec$ell_ref[["min"]]) / scaleFactor
  )
}

# Attach group labels to a natural-scale observation parameter of a fit (tau as a 2-row matrix by
# group; df/k named by group), as accepted by 'par'.
.label_obs_par <- function(name, value, spec) {
  if (is.null(spec)) return(value)
  if (name == "tau" && !is.null(spec$tau_labels) && length(value) == 2 * length(spec$tau_labels)) {
    value <- matrix(value, nrow = 2, dimnames = list(c("x", "y"), spec$tau_labels))
  } else if (name %in% c("df", "k")) {
    l <- spec[[paste0(name, "_labels")]]
    if (!is.null(l) && length(value) == length(l)) names(value) <- l
  }
  value
}

# Which tail parameters are used by the model: any_t + 2 * any_Huber (0 = all Gaussian)
.obs_dist_flag <- function(obs_dist) {
  obs_dist <- .norm_obs_dist(obs_dist)
  as.integer(("t" %in% obs_dist) + 2L * ("Huber" %in% obs_dist))
}
