
<!-- README.md is generated from README.Rmd. Please edit that file -->

# {langevinSSM}

#### Habitat-driven Langevin Diffusion with Spatial Uncertainty

`{langevinSSM}` is an R package for simulating and fitting the
habitat-driven Langevin diffusion to animal tracking data subject to
location measurement error and temporal irregularity. The habitat-driven
Langevin diffusion can provide inferences about habitat selection and
utilization distributions. The package provides tools for simulating
animal movement paths (`simLangevin`) and fitting the Langevin diffusion
model to observed tracking data (`fitLangevin`). Location measurement
error can take the form of either (older) Argos Least Squares-based
locations or (newer) Argos Kalman Filter-based locations with error
ellipse information. The Langevin diffusion is a continuous-time model
in state-space form that estimates the underlying movement process while
accounting for location measurement error and associated uncertainty in
the spatial (habitat) covariates. Template Model Builder {TMB} is used
for fast estimation.

## Installation

One can install the `langevinSSM` package from CRAN using the following
command:

``` r
install.packages("langevinSSM") 
```

Alternatively, one can install the package from GitHub using the
`remotes` package:

``` r
remotes::install_github("bmcclintock/langevinSSM")
```

## Usage

To simulate animal movement paths using the Langevin diffusion model,
one can use the `simLangevin` function. For example:

``` r
library(langevinSSM)
library(ggplot2)
library(terra)
library(patchwork)

# Simulate an underdamped Langevin diffusion path

par <- list(beta = c(-4, 6, 5, -0.1), # habitat selection coefficients
            sigma = 5, # diffusion (or speed) parameter
            gamma = 0.5) # autocorrelation parameter

## exampleCovs is a list of four spatial covariates (e.g., habitat features) that loads with the package
simDat <- simLangevin(model = "underdamped",
                      par = par,
                      spatialCovs = exampleCovs,
                      nbAnimals = 3)

# Simulate an underdamped Langevin diffusion path with measurement error
measurementError <- list(smaj.sd = 1.5,      # sd of semi-major axis of error ellipse
                         smin.sd = 0.75,     # sd of semi-minor axis of error ellipse
                         eor = c(0,180)) # range of ellipse orientation (in degrees from north)
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")
exampleDat <- simLangevin(model = "underdamped",
                         par = par,
                         spatialCovs = exampleCovs,
                         nbAnimals = 3,
                         obsPerAnimal = 500,
                         measurementError = measurementError)
```

To fit the Langevin diffusion model to observed tracking data, one can
use the `fitLangevin` function. For example:

``` r
# Fit the underdamped Langevin diffusion model to simulated data with measurement error
## setting calcOSA = TRUE will calculate one-step-ahead residuals for model diagnostics
fit <- fitLangevin(model = "underdamped",
                   data = exampleDat,
                   spatialCovs = exampleCovs,
                   silent = TRUE,
                   calcOSA = TRUE)  

fit
#> 
#> Habitat-Driven Langevin Diffusion Model
#> =======================================
#> Model type:        Underdamped 
#> Convergence:       Successful 
#> Max Log-Likelihood: -2061.97 
#> Optimization time:  1.76 seconds
#> 
#> Parameter Estimates (Natural Scale):
#> ---------------------------------------
#>        Estimate Std. Error
#> beta_1  -3.5438      1.250
#> beta_2   8.6136      2.295
#> beta_3   7.5114      1.947
#> beta_4  -0.3084      0.324
#> sigma    4.3092      0.500
#> gamma    0.5850      0.145
#> rho_o    0.0000      0.000
#> tau_1    1.0000      0.000
#> tau_2    1.0000      0.000
#> psi      1.0000      0.000

# calculate the estimated UD
UD <- getUD(spatialCovs = exampleCovs, beta = fit$estimates$natural[1:length(exampleCovs)])

# plot the estimated (log) UD with the observed and estimated locations
plot(fit, spatialCovs = exampleCovs, data = exampleDat)
```

![](man/figures/README-unnamed-chunk-5-1.png)<!-- -->

``` r

# plot residuals to check model fit
p <- plotResiduals(fit)
p$qq_x + p$qq_y + p$acf_x + p$acf_y + plot_layout(ncol=2)
```

![](man/figures/README-unnamed-chunk-5-2.png)<!-- -->

``` r

# calculate the true utiliziation distribution
trueUD <- getUD(spatialCovs = exampleCovs, beta = par$beta)

# calculate similarity of true and estimated UDs using Bhattacharyya's affinity
rasterOverlap(exp(UD), exp(trueUD))
#> [1] 0.9088076
```

## Citation

If you use `{langevinSSM}` in your research, please cite it as follows:

    To cite package 'langevinSSM' in publications use:

      Dupont, F., McClintock, B.T., Fischer, J.-O., Marcoux, M., Hussey,
      N., and Auger-Méthé, M. (2025). Inferring resource selection and
      utilization distributions from irregular and error-prone animal
      tracking data using the habitat-driven Langevin diffusion.

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {Inferring resource selection and utilization distributions from irregular and error-prone animal tracking data using the habitat-driven Langevin diffusion},
        author = {Fanny Dupont and Brett T. McClintock and Jan-Ole Fischer and Marianne Marcoux and Nigel Hussey and Marie Auger-Méthé},
        journal = {TBD},
        year = {2025},
      }

    Additions and modifications to langevinSSM are frequent, to help with
    reproducibility of output please cite its version number. This is
    'langevinSSM' version 0.0.1
