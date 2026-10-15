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

# ---------------------------------------------------------------------------------------------------
# Covariates and barrier
#   * Grid: 821.408 m cells (the native ~30" ETOPO spacing in UTM) extending 60 km beyond the original
#     extent (every location and every fitted mu lies > 60 km inside the raster edge; see the check
#     after fit0). Finer grids (250 m) only bilinearly interpolate the same 60" data, which flattens
#     the cell-scale gradient and drives beta_bathy towards 0.
#   * Bathymetry: NOAA ETOPO 2022 from getCovs(), raw (NOT smoothed: smoothing flips beta_bathy positive).
#   * Coast: GSHHG v2.3.7 full-resolution shoreline (L1 = land), used as-is (no buffering or editing).
#     prepBarrier() computes the exact signed distance from each cell centre to the shoreline polygons,
#     so the zero contour of the barrier is the coastline itself (Delporte et al. 2025).
#     The ETOPO land mask (bathy > 0) is only ~1.8 km resolution and is not used for the barrier.
# ---------------------------------------------------------------------------------------------------
target_res <- 821.408
covs_buffered <- getCovs(buffer_meters = 200000, target_res = target_res)
bathy <- covs_buffered$bathy
template <- rast(bathy)

# GSHHG land polygons, cropped and projected to the template (downloaded once [142 MB], then cached)
coast_file <- "applications/narwhal/gshhg_f_L1_utm17.gpkg"
if (!file.exists(coast_file)) {
  gz <- file.path(tempdir(), "gshhg-shp-2.3.7.zip")
  options(timeout = max(1800, getOption("timeout")))
  download.file("https://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip", gz, mode = "wb")
  shp <- unzip(gz, list = TRUE)$Name
  shp <- shp[grepl("GSHHS_shp/f/GSHHS_f_L1\\.", shp)]
  unzip(gz, files = shp, exdir = tempdir())
  unlink(gz)
  e20 <- as.vector(ext(template) + 20000)   # template + 20 km
  box <- st_as_sfc(st_bbox(c(xmin = e20[["xmin"]], ymin = e20[["ymin"]], xmax = e20[["xmax"]], ymax = e20[["ymax"]]),
                           crs = st_crs(crs(template))))
  land <- st_read(file.path(tempdir(), "GSHHS_shp/f/GSHHS_f_L1.shp"),
                  wkt_filter = st_as_text(st_transform(box, 4326)), quiet = TRUE)   # read only the study area
  land <- st_make_valid(st_transform(land, st_crs(box)))
  land <- st_crop(land, box)
  st_write(land, coast_file, quiet = TRUE)
}
land <- st_read(coast_file, quiet = TRUE)

# signed distance to the shoreline (m): positive in water (allowed), negative on land (restricted)
barrier <- prepBarrier(land, template = template)   # ~15 s at 821 m

# binary mask on the same grid (water = 1 allowed, land = 0 restricted), for plotting and maskLeakage
mask <- rasterize(vect(land), template, field = 0, background = 1)
names(mask) <- "mask"
plotRaster(mask,extent=origextent)

d2c <- barrier

# Optional Gaussian smoothing of bathymetry (sd in m). Keep at 0: a 5 km kernel removes the
# fine-scale depth gradient that drives the narwhal signal and makes beta_bathy positive (and the fit
# clearly worse: ~17-24 log-likelihood units).
bathySmooth <- 0
if (bathySmooth > 0) {
  w <- focalMat(bathy, bathySmooth, "Gauss")
  bathy <- focal(bathy, w, fun = "sum", na.rm = TRUE, expand = TRUE) /
    focal(bathy * 0 + 1, w, fun = "sum", na.rm = TRUE, expand = TRUE)
  names(bathy) <- "bathy"
}

plot(bathy, main = "bathymetry")
plotRaster(d2c,extent=origextent)
plot(barrier, main = "signed distance field")

narRaw <- read_csv("https://raw.githubusercontent.com/Fanny-Dupont/Langevin_SSM/284de22e4746cc393a9019f0f13c8ee6fef45dcc/Code_Case_Study/Narwhal_Case_Study.csv")
narRaw <- narRaw %>% mutate(X=X*1000,Y=Y*1000,smaj=smaj*1000,smin=smin*1000)

# Set Fastloc to LS with 50m error and Argos to KF
narRaw$x.err <- narRaw$y.err <- NA
floc_idx <- which(narRaw$loc_class=="G")
narRaw$smaj[floc_idx] <- narRaw$smin[floc_idx] <- narRaw$eor[floc_idx] <- NA
narRaw$x.err[floc_idx] <- narRaw$y.err[floc_idx] <- 50

narDat <- formatData(narRaw, id="ID", date="datetime_UTC", coord=c("X","Y"), lc = "loc_class")

covs <- list(bathy=scale(bathy), barrier = barrier)

leaksDat <- maskLeakage(narDat,mask, level=1, tolerance=res(mask)[1],coord=c("x","y"))

plotRaster(mask,extent = origextent)+geom_point(aes(x=x,y=y),data=narDat,col="#E69F00")

# ---------------------------------------------------------------------------------------------------
# Observation model: Gaussian. With t errors the full-data fit fails even without a barrier
# (false convergence, inner Hessian not positive definite at the optimum), because of conflicting
# Argos outliers (e.g. track 17015-20). The Gaussian model is log-concave in the latent states;
# psi, tau and rho_o still rescale the Argos ellipses and the 50 m Fastloc error.
# ---------------------------------------------------------------------------------------------------
obs_dist <- "Gaussian"

# Starting values near the full-data optimum. From par = list(psi=1, tau=c(1,1), rho_o=0) the optimiser
# stops in a poor local optimum (false convergence, tau_2 -> 0.01, log-likelihood ~110 units lower).
fit0 <- fitLangevin(narDat, model = model, spatialCovs = covs["bathy"],
                    par=list(beta=0, sigma=3000, gamma=1.15, psi=21, tau=c(1.1,1.3), rho_o=0.06),
                    obs_dist = obs_dist,
                    scaleFactor=scaleFactor,
                    control=list(trace=1),
                    silent = TRUE)
fit0
# Covariates are 0 (with zero gradient) outside the hull of raster cell centres, so check that the
# fitted locations stay well inside the raster (the prepareRaster "dangerously close" warning is
# based on obs +/- 3 error SDs, i.e. the huge Argos ellipses, not on the fitted locations).
muEdge <- function(fit, r) {
  m <- fit$estimates$random$mu$est; e <- as.vector(ext(r))
  min(m$mu.x - e[1], e[2] - m$mu.x, m$mu.y - e[3], e[4] - m$mu.y) / 1000
}
cat(sprintf("min distance from fitted mu to raster edge: %.0f km\n", muEdge(fit0, bathy)))
leaks0 <- maskLeakage(fit0,mask, level=1, tolerance=res(m))

# ---------------------------------------------------------------------------------------------------
# Barrier fit
#   * barrierOnset: the force ramps in smoothly (C4) over the first barrierOnset metres inland, then
#     continues as the usual quadratic wall. "res" = one cell (821 m here), about the coastline accuracy.
#     Choose it for the species: for animals that never go on land (narwhals) it is a small numerical
#     device; for animals that haul out it should be larger, or use a soft (small fixed / estimated)
#     lambda. barrierOnset = 0 is Delporte et al.'s penalty.
#   * lambda (natural scale, 1/m^2): fixed at kappa * gamma^2 / sigma^2 from the barrier-free fit, with
#     kappa = 1. kappa = sigma^2 * lambda / gamma^2 is the dimensionless stiffness of the wall in units
#     of the movement process (kappa = 1/4 is critical damping). The package default lambda = "tied"
#     (the per-step stability limit) grows like 1/dt and breaks the Laplace approximation at the
#     2-3 min Argos/Fastloc steps in these data.
#     Set lambdaMode <- "estimate" to estimate lambda instead (starting from the kappa = 1 value).
#   * Starting values: all of fit0's parameters, including the smoothed locations (mu) and velocities
#     (vel). Dropping mu/vel restarts the latent path from the raw observations, and the fit then ends in
#     false convergence at a much worse optimum (~400 log-likelihood units lower).
# ---------------------------------------------------------------------------------------------------
narPar <- getPar(fit0)
kappa <- 1
lambda_kappa <- kappa * narPar$gamma^2 / narPar$sigma^2
lambdaMode <- "fixed"   # "fixed" (kappa = 1) or "estimate"
if (lambdaMode == "estimate") {
  narPar$lambda <- lambda_kappa
  lambdaArg <- "estimate"
} else {
  lambdaArg <- lambda_kappa
}
narFit <- fitLangevin(narDat, model = model, spatialCovs = covs, barrier = "barrier",
                      par = narPar,
                      obs_dist = obs_dist,
                      scaleFactor=scaleFactor,
                      lambda = lambdaArg,
                      barrierOnset = "res",
                      control=list(trace=1),
                      silent = TRUE)
narFit
cat(sprintf("min distance from fitted mu to raster edge: %.0f km\n", muEdge(narFit, bathy)))
plot(narFit,spatialCovs=covs,data=narDat,extent=origextent,normalize=TRUE,maskRast=mask)
narLeaks <- maskLeakage(narFit,mask,level=1,tolerance = res(mask)[1])

narPar_t <- narPar
narPar_t$df <- 30
narFit_t <- fitLangevin(narDat, model = model, spatialCovs = covs, barrier = "barrier",
                      par = narPar_t,
                      obs_dist = "t",
                      scaleFactor=scaleFactor,
                      lambda = lambdaArg,
                      barrierOnset = "res",
                      control=list(trace=1),
                      silent = TRUE)
narFit_t
cat(sprintf("min distance from fitted mu to raster edge: %.0f km\n", muEdge(narFit_t, bathy)))
plot(narFit_t,spatialCovs=covs,data=narDat,extent=origextent,normalize=TRUE,maskRast=mask)
narLeaks_t <- maskLeakage(narFit_t,mask,level=1,tolerance = res(mask)[1])

AIC(fit0,narFit,narFit_t)
