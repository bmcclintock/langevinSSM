#remotes::install_github('bmcclintock/langevinSSM')
library(langevinSSM)
library(terra)
library(sf)
library(readr)
library(dplyr)
library(ggplot2)

set.seed(1,kind="Mersenne-Twister",normal.kind = "Inversion")

source("applications/narwhal/getCovs.R")
origextent <- ext(c(469606.176308853, 628137.927473616, 7968123.93718224, 8212082.12420781))

scaleFactor <- 1000

model <- "underdamped"

# expand raster boundaries by 200 kilometers in every direction
covs_buffered <- getCovs(buffer_meters = 200000, target_res = 821.408)

bathy <- covs_buffered$bathy

coast_file <- "applications/narwhal/gshhg_f_L1_utm17.gpkg"
if (!file.exists(coast_file)) {
  gz <- file.path(tempdir(), "gshhg-shp-2.3.7.zip")
  options(timeout = max(1800, getOption("timeout")))
  download.file("https://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip", gz, mode = "wb")
  shp <- unzip(gz, list = TRUE)$Name
  unzip(gz, files = shp[grepl("GSHHS_shp/f/GSHHS_f_L1\\.", shp)], exdir = tempdir())
  box <- st_as_sfc(st_bbox(as.vector(ext(bathy) + 20000), crs = st_crs(crs(bathy))))
  land <- st_read(file.path(tempdir(), "GSHHS_shp/f/GSHHS_f_L1.shp"), wkt_filter = st_as_text(st_transform(box, 4326)), quiet = TRUE)
  land <- st_crop(st_make_valid(st_transform(land, st_crs(box))), box)
  st_write(land, coast_file, quiet = TRUE)
}
land <- st_read(coast_file, quiet = TRUE)

# signed distance to the coastline: positive in water (allowed), negative on land (restricted)
barrier <- prepBarrier(land, template = bathy, restricted = TRUE)

# binary mask on the same grid (water = 1, land = 0), for routeTracks and plotting
mask <- rasterize(vect(land), rast(bathy), field = 0, background = 1)

# bathymetry, with land elevations set to 0 m
covs <- list(bathy=scale(clamp(bathy, upper = 0)), barrier = barrier)

# hard wall (narwhals cannot go on land): with lambda = "tied" the barrier strength is set at every time
# step to the largest numerically stable value, and getUD() gives zero use of land. The barrier force
# ramps in over one cell (821 m) inland, about the accuracy of the coastline and the Fastloc fixes.
lambda <- "tied"
barrierOnset <- "res"

narRaw <- read_csv("https://raw.githubusercontent.com/Fanny-Dupont/Langevin_SSM/284de22e4746cc393a9019f0f13c8ee6fef45dcc/Code_Case_Study/Narwhal_Case_Study.csv")
narRaw <- narRaw %>% mutate(X=X*1000,Y=Y*1000,smaj=smaj*1000,smin=smin*1000)

# Set Fastloc to LS with 50m error and Argos to KF
narRaw$x.err <- narRaw$y.err <- NA
floc_idx <- which(narRaw$loc_class=="G")
narRaw$smaj[floc_idx] <- narRaw$smin[floc_idx] <- narRaw$eor[floc_idx] <- NA
narRaw$x.err[floc_idx] <- narRaw$y.err[floc_idx] <- 50

# individuals (track segments "17001-1", "17001-2", ... belong to animal "17001")
narRaw$animal <- sub("-.*", "", narRaw$ID)
animals <- sort(unique(narRaw$animal))

########################################################################
# Stage I: candidate observation models (Fastloc GPS always Gaussian)
#  A: Argos ellipse Gaussian, psi;                      GPS tau (x = y)
#  B: Argos ellipse t,        psi, df;                  GPS tau (x = y)
#  C: B + error ellipse power law (psi_maj, b_maj, b_min)
#  D: C + separate GPS tau_x, tau_y + correlation rho_o
#  E: C + separate t df for the good (3/2/1) and poor (0/A/B) Argos classes
#  F: D + separate t df for the good (3/2/1) and poor (0/A/B) Argos classes
#
# In E and F, 'df' is given for 'epar' (used by classes 3, 2 and 1) and for
# classes 0, A and B; map = factor(c(1,2,2,2)) ties 0/A/B to one df. The poor
# classes have heavier tails than the good ones (df ~2-3 vs ~5-10).
########################################################################

obsModels <- list(
  A = list(obs_dist = c(sderr = "Gaussian", epar = "Gaussian"),
           par = list(psi = 1.5, tau = c(1.27, 1.27)),
           map = list(tau = factor(c(1, 1)))),
  B = list(obs_dist = c(sderr = "Gaussian", epar = "t"),
           par = list(psi = 1.5, tau = c(1.27, 1.27), df = 2.25),
           map = list(tau = factor(c(1, 1)))),
  C = list(obs_dist = c(sderr = "Gaussian", epar = "t"),
           par = list(psi = 1.5, psi_maj = 0.311, b_maj = 0.287, b_min = 0.331, tau = c(1.27, 1.27), df = 2.25),
           map = list(tau = factor(c(1, 1)))),
  D = list(obs_dist = c(sderr = "Gaussian", epar = "t"),
           par = list(psi = 1.5, psi_maj = 0.311, b_maj = 0.287, b_min = 0.331, tau = c(1.19, 1.35), rho_o = 0.065, df = 2.25),
           map = NULL),
  E = list(obs_dist = c(sderr = "Gaussian", epar = "t"),
           par = list(psi = 1.5, psi_maj = 0.311, b_maj = 0.287, b_min = 0.331, tau = c(1.27, 1.27),
                      df = c(epar = 5, "0" = 2.25, A = 2.25, B = 2.25)),
           map = list(tau = factor(c(1, 1)), df = factor(c(1, 2, 2, 2)))),
  F = list(obs_dist = c(sderr = "Gaussian", epar = "t"),
           par = list(psi = 1.5, psi_maj = 0.311, b_maj = 0.287, b_min = 0.331, tau = c(1.19, 1.35), rho_o = 0.065,
                      df = c(epar = 5, "0" = 2.25, A = 2.25, B = 2.25)),
           map = list(df = factor(c(1, 2, 2, 2))))
)

# starting values for the movement parameters and psi
moveStarts <- list(list(beta = -0.37, sigma = 3000, gamma = 1.07),
                   list(beta = -0.40, sigma = 2000, gamma = 0.80),
                   list(beta = -0.20, sigma = 3500, gamma = 1.20))
psiStarts <- c(1.5, 4, 8)

# a fit is usable if all standard errors are finite and the covariance matrix
# of the movement parameters (the ones pooled in stage II) is positive definite
movePars <- c("beta_bathy", "log_sigma", "log_gamma")
usable <- function(fit){
  if(inherits(fit,"error")) return(FALSE)
  est <- fit$estimates$working
  i <- match(movePars, rownames(est))
  V <- as.matrix(fit$covariance$working)[i, i]
  all(is.finite(est[,2])) && all(is.finite(V)) && all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0)
}

# fit observation model 'om' to data 'dat' from each set of starting values; keep the usable fit with the highest log-likelihood
fitObsModel <- function(dat, om){
  best <- NULL
  for(ms in moveStarts){
    for(psi in psiStarts){
      par <- c(ms, obsModels[[om]]$par)
      par$psi <- psi
      par$mu <- cbind(dat$mu.x_pr, dat$mu.y_pr)

      fit <- tryCatch(suppressWarnings(suppressMessages(fitLangevin(dat, model = model, obs_dist = obsModels[[om]]$obs_dist,
                                                                    spatialCovs = covs, barrier = "barrier",
                                                                    lambda = lambda, barrierOnset = barrierOnset,
                                                                    par = par, map = obsModels[[om]]$map,
                                                                    scaleFactor = scaleFactor,
                                                                    control = list(trace = 0, iter.max = 2000, eval.max = 2000, rel.tol = 1e-8),
                                                                    silent = TRUE))), error = function(e) e)

      if(usable(fit) && (is.null(best) || fit$objective < best$objective)) best <- fit
    }
  }
  best
}

AICtab <- NULL
bestFits <- res <- dats <- list()

for(a in animals){
  message("Individual ",a)

  dat <- formatData(narRaw[narRaw$animal==a,], id="ID", date="datetime_UTC", coord=c("X","Y"), lc = "loc_class")
  dat <- routeTracks(dat, mask)

  fits <- list()
  for(om in names(obsModels)){
    message("   observation model ",om)
    fits[[om]] <- fitObsModel(dat, om)
  }
  fits <- fits[!sapply(fits, is.null)]
  if(!length(fits)) { message("   no usable fit for ",a); next }

  tab <- data.frame(id = a,
                    model = names(fits),
                    npar = sapply(fits, function(f) nrow(f$estimates$working)),
                    logLik = sapply(fits, function(f) -f$objective))
  tab$AIC <- 2 * tab$npar - 2 * tab$logLik
  tab$deltaAIC <- tab$AIC - min(tab$AIC)

  # best-supported model: the simplest model within 2 AIC units of the lowest AIC
  close <- tab[tab$deltaAIC <= 2,]
  bestModel <- close$model[which.min(close$npar)]
  tab$best <- tab$model == bestModel
  print(tab, row.names = FALSE)

  AICtab <- rbind(AICtab, tab)
  bestFits[[a]] <- fits[[bestModel]]
  dats[[a]] <- dat
  res[[a]] <- residuals(bestFits[[a]],data=dat,spatialCovs=covs)
  print(plot(res[[a]]))
}

AICtab
sapply(bestFits, function(f) paste(f$conditions$obs_dist, collapse = "/"))
lapply(bestFits, function(f) round(f$estimates$natural, 3))

resAll <- do.call(rbind, res)
datAll <- do.call(rbind, dats)
attr(resAll, "tests") <- NULL
attr(resAll, "obs_dist") <- unique(unlist(lapply(res, attr, "obs_dist")))
resAll$id <- datAll$id <- sub("-.*", "", datAll$id)   # 'id' = individual instead of track segment
lcGroups <- list("Argos 3/2/1" = c("3", "2", "1"), "Argos 0/A/B" = c("0", "A", "B"), GPS = "G")
narResFig <- plot(resAll, by = list("id", lc = lcGroups), data = datAll)
for(g in names(narResFig)) { print(narResFig[[g]]$qq_x); print(narResFig[[g]]$qq_y) }

saveRDS(list(AICtab = AICtab, bestFits = bestFits), "applications/narwhal/narwhalHier_stage1.rds")

########################################################################
# Stage II: pool the movement parameters (beta, sigma, gamma) across
# individuals; each individual keeps its own observation model
########################################################################

narHier <- hierLangevin(bestFits, params = "movement")
narHier
narHier$estimates$natural

popUD <- getUD(spatialCovs=covs,fit=narHier,nSims=1000,maskRast=mask,extent=origextent)

# individual-level (shrunk) estimates
narHier$estimates$random
