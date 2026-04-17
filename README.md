
<!-- README.md is generated from README.Rmd. Please edit that file -->

# {langevinSSM}

#### Habitat-driven Langevin Diffusion with Spatial Uncertainty

`{langevinSSM}` is an R package for simulating and fitting the
habitat-driven Langevin diffusion to animal tracking data subject to
location measurement error and temporal irregularity. The habitat-driven
Langevin diffusion provides inferences about both habitat selection and
utilization distributions. The package provides tools for simulating
animal movement paths (`simLangevin`) and fitting the Langevin diffusion
model to observed tracking data (`fitLangevin`). Location measurement
error can take the form of (older) Argos Least Squares-based locations,
(newer) Argos Kalman Filter-based locations with error ellipse
information, or general x- and y-axis errors (e.g. for GPS data). The
Langevin diffusion is a continuous-time model in state-space form that
estimates the underlying movement process while accounting for location
measurement error and associated uncertainty in the spatial (habitat)
covariates. Template Model Builder {TMB} is used for fast estimation.

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

# calculate the true utiliziation distribution
## exampleCovs is a list of four spatial covariates (e.g., habitat features) that loads with the package
trueUD <- getUD(spatialCovs = exampleCovs, beta = par$beta)
```

![](man/figures/README-sim-1.png)<!-- -->

``` r

simDat <- simLangevin(model = "underdamped",
                      par = par,
                      spatialCovs = exampleCovs,
                      nbAnimals = 3)

head(simDat)
#>   id date   dt        x        y smaj smin eor x.err y.err     mu.x     mu.y
#> 1  1 0.00 0.00 1024.000 990.0000   NA   NA  NA    NA    NA 1024.000 990.0000
#> 2  1 0.01 0.01 1024.067 990.0602   NA   NA  NA    NA    NA 1024.067 990.0602
#> 3  1 0.02 0.01 1024.130 990.1191   NA   NA  NA    NA    NA 1024.130 990.1191
#> 4  1 0.03 0.01 1024.195 990.1747   NA   NA  NA    NA    NA 1024.195 990.1747
#> 5  1 0.04 0.01 1024.256 990.2287   NA   NA  NA    NA    NA 1024.256 990.2287
#> 6  1 0.05 0.01 1024.315 990.2791   NA   NA  NA    NA    NA 1024.315 990.2791
#>      vel.x    vel.y
#> 1 6.648996 6.362147
#> 2 6.341045 5.739271
#> 3 6.837852 5.723017
#> 4 6.165345 5.344009
#> 5 5.949173 5.076844
#> 6 5.844776 5.189753

# Simulate an underdamped Langevin diffusion path with measurement error
measurementError <- list(smaj.sd = 1.5,      # sd of semi-major axis of error ellipse
                         smin.sd = 0.75,     # sd of semi-minor axis of error ellipse
                         eor.lim = c(0,180)) # range of ellipse orientation (in degrees from north)

exampleDat <- simLangevin(model = "underdamped",
                         par = par,
                         spatialCovs = exampleCovs,
                         nbAnimals = 3,
                         obsPerAnimal = 500,
                         measurementError = measurementError)

head(exampleDat)
#>   id date   dt        x        y      smaj      smin       eor x.err y.err
#> 1  1 0.00 0.00 1023.610 988.7225 3.1904213 0.8238468 0.4407540    NA    NA
#> 2  1 0.01 0.01 1024.056 990.0353 0.2602149 0.1929297 0.8419944    NA    NA
#> 3  1 0.02 0.01 1022.943 991.0823 1.9072127 0.6294832 2.2160981    NA    NA
#> 4  1 0.03 0.01 1024.243 990.1201 0.1952989 0.0469546 2.2100417    NA    NA
#> 5  1 0.04 0.01 1024.996 990.5625 1.4343749 1.1175133 2.8685380    NA    NA
#> 6  1 0.05 0.01 1024.489 990.9043 0.8349596 0.1392169 0.4690632    NA    NA
#>       mu.x     mu.y    vel.x    vel.y
#> 1 1024.000 990.0000 6.648996 6.362147
#> 2 1024.067 990.0602 6.341045 5.739271
#> 3 1024.130 990.1191 6.837852 5.723017
#> 4 1024.195 990.1747 6.165345 5.344009
#> 5 1024.256 990.2287 5.949173 5.076844
#> 6 1024.315 990.2791 5.844776 5.189753
```

To fit the Langevin diffusion model to observed tracking data, one can
use the `formatData` and `fitLangevin` functions. For example:

``` r
# unformatDat is example data appropriate for formatData that loads with the package
head(unformatDat)
#>   id                date        x        y      smaj      smin       eor x.err
#> 1  1 2026-04-14 00:00:00 1023.610 988.7225 3.1904213 0.8238468  25.25334    NA
#> 2  1 2026-04-14 00:00:36 1024.056 990.0353 0.2602149 0.1929297  48.24272    NA
#> 3  1 2026-04-14 00:01:12 1022.943 991.0823 1.9072127 0.6294832 126.97307    NA
#> 4  1 2026-04-14 00:01:48 1024.243 990.1201 0.1952989 0.0469546 126.62606    NA
#> 5  1 2026-04-14 00:02:24 1024.996 990.5625 1.4343749 1.1175133 164.35512    NA
#> 6  1 2026-04-14 00:03:00 1024.489 990.9043 0.8349596 0.1392169  26.87534    NA
#>   y.err
#> 1    NA
#> 2    NA
#> 3    NA
#> 4    NA
#> 5    NA
#> 6    NA

# format the data for fitLangevin
exampleDat <- formatData(unformatDat, time.unit = "hours")

head(exampleDat)
#>   id                date   dt        x        y   lc      smaj      smin
#> 1  1 2026-04-14 00:00:00 0.00 1023.610 988.7225 <NA> 3.1904213 0.8238468
#> 2  1 2026-04-14 00:00:36 0.01 1024.056 990.0353 <NA> 0.2602149 0.1929297
#> 3  1 2026-04-14 00:01:12 0.01 1022.943 991.0823 <NA> 1.9072127 0.6294832
#> 4  1 2026-04-14 00:01:48 0.01 1024.243 990.1201 <NA> 0.1952989 0.0469546
#> 5  1 2026-04-14 00:02:24 0.01 1024.996 990.5625 <NA> 1.4343749 1.1175133
#> 6  1 2026-04-14 00:03:00 0.01 1024.489 990.9043 <NA> 0.8349596 0.1392169
#>         eor x.err y.err
#> 1 0.4407540    NA    NA
#> 2 0.8419944    NA    NA
#> 3 2.2160981    NA    NA
#> 4 2.2100417    NA    NA
#> 5 2.8685380    NA    NA
#> 6 0.4690632    NA    NA

# Fit the overdamped Langevin diffusion model to simulated data with measurement error
fit_over <- fitLangevin(model = "overdamped",
                   data = exampleDat,
                   spatialCovs = exampleCovs,
                   silent = TRUE)  

fit_over
#> 
#> Habitat-Driven Langevin Diffusion Model
#> =======================================
#> Model type:        Overdamped 
#> Convergence:       Successful 
#> Max Log-Likelihood: -2428.077 
#> Optimization time:  0.35 seconds
#> 
#> Parameter Estimates (Natural Scale):
#> ---------------------------------------
#>           Estimate Std. Error
#> beta_cov1   23.597      7.324
#> beta_cov2   15.901      7.255
#> beta_cov3   -4.850      6.411
#> beta_d2c   -13.782      2.137
#> sigma        1.143      0.035
#> rho_o        0.000      0.000
#> tau_1        1.000      0.000
#> tau_2        1.000      0.000
#> psi          1.000      0.000

# Fit the underdamped Langevin diffusion model
fit_under <- fitLangevin(model = "underdamped",
                   data = exampleDat,
                   spatialCovs = exampleCovs,
                   silent = TRUE)  

fit_under
#> 
#> Habitat-Driven Langevin Diffusion Model
#> =======================================
#> Model type:        Underdamped 
#> Convergence:       Successful 
#> Max Log-Likelihood: -2061.97 
#> Optimization time:  0.61 seconds
#> 
#> Parameter Estimates (Natural Scale):
#> ---------------------------------------
#>           Estimate Std. Error
#> beta_cov1  -3.5437      1.250
#> beta_cov2   8.6136      2.295
#> beta_cov3   7.5113      1.948
#> beta_d2c   -0.3084      0.324
#> sigma       4.3092      0.500
#> gamma       0.5850      0.145
#> rho_o       0.0000      0.000
#> tau_1       1.0000      0.000
#> tau_2       1.0000      0.000
#> psi         1.0000      0.000
```

### Post-processing functions

#### Utilization distribution

``` r
# calculate the estimated UD
UD <- getUD(spatialCovs = exampleCovs, 
            fit = fit_under, 
            nSims = 1000, # Monte Carlo simulation
            show_progress = FALSE)
```

![](man/figures/README-ud-1.png)<!-- -->

``` r

p_UD <- plotUD(UD)

# UD relative uncertainty (Delta method approximation)
p_UD$CV_delta
```

![](man/figures/README-ud-2.png)<!-- -->

``` r

# UD relative uncertainty (Monte Carlo simulation)
p_UD$CV_sim
```

![](man/figures/README-ud-3.png)<!-- -->

``` r

# plot the estimated (log) UD with the observed and estimated locations
plot(fit_under, spatialCovs = exampleCovs, data = exampleDat)
```

![](man/figures/README-ud-4.png)<!-- -->

#### Other S3 methods for `fitLangevin` objects

``` r
# fixed effect estimates
coef(fit_under) 
#>  beta_cov1  beta_cov2  beta_cov3   beta_d2c      sigma      gamma      rho_o 
#> -3.5437139  8.6136199  7.5113365 -0.3083557  4.3092006  0.5849726  0.0000000 
#>      tau_1      tau_2        psi 
#>  1.0000000  1.0000000  1.0000000

# confidence intervals for fixed effects
confint(fit_under) 
#>                2.5 %     97.5 %
#> beta_cov1 -5.9945650 -1.0928627
#> beta_cov2  4.1147592 13.1124807
#> beta_cov3  3.6940963 11.3285766
#> beta_d2c  -0.9436450  0.3269335
#> sigma      3.3292213  5.2891798
#> gamma      0.3012216  0.8687236
#> rho_o      0.0000000  0.0000000
#> tau_1      1.0000000  1.0000000
#> tau_2      1.0000000  1.0000000
#> psi        1.0000000  1.0000000

# confidence intervals for true locations
mu_ci <- confint(fit_under, type= "mu") 

head(mu_ci)
#>   id time_step     mu.x mu.x_2.5% mu.x_97.5%     mu.y mu.y_2.5% mu.y_97.5%
#> 1  1         1 1023.966  1023.817   1024.115 990.0594  989.9307   990.1882
#> 2  1         2 1024.019  1023.889   1024.150 990.1089  989.9969   990.2208
#> 3  1         3 1024.073  1023.959   1024.187 990.1574  990.0591   990.2557
#> 4  1         4 1024.126  1024.026   1024.226 990.2050  990.1168   990.2932
#> 5  1         5 1024.180  1024.091   1024.268 990.2515  990.1697   990.3333
#> 6  1         6 1024.233  1024.153   1024.313 990.2968  990.2182   990.3754

# AIC for comparing models with different fixed effects
AIC(fit_under) 
#> [1] 4135.94

# BIC for comparing models with different fixed effects
BIC(fit_under) 
#> [1] 4167.819
```

#### One-step-ahead residuals

``` r
# calculate one-step-ahead residuals for model diagnostics
res_under <- residuals(fit_under, data = exampleDat, spatialCovs = exampleCovs, ncores = 3)
res_under
#> 
#> === One-Step-Ahead (OSA) Residuals ===
#> Total observations: 1500 
#> Number of tracks:   3 
#> 
#> ---- Goodness-of-Fit Tests ----
#>  metric  statistic   p.value
#>    KS_x 0.01598432 0.8389506
#>    KS_y 0.02434843 0.3373330
#>  KS_mah 0.01965255 0.6097253
#>    LB_x 7.88961085 0.3424286
#>    LB_y 3.28849949 0.8570934
#>  LB_mah 5.37975199 0.6137248
#> -------------------------------
#> 
#> Residual Summary:
#>    residual.x           residual.y      
#>  Min.   :-3.0898455   Min.   :-3.57833  
#>  1st Qu.:-0.6737108   1st Qu.:-0.66399  
#>  Median :-0.0125678   Median : 0.01314  
#>  Mean   : 0.0003674   Mean   : 0.01898  
#>  3rd Qu.: 0.6900477   3rd Qu.: 0.71835  
#>  Max.   : 3.1302167   Max.   : 2.97038  
#>  NA's   :3            NA's   :3

# plot residuals to check model fit
p_under <- plot(res_under)
p_under$qq_x + p_under$qq_y + p_under$acf_x + p_under$acf_y + plot_layout(ncol=2)
```

![](man/figures/README-osa-1.png)<!-- -->

``` r

# can be used to compare "underdamped" vs "overdamped" models
res_over <- residuals(fit_over, data = exampleDat, spatialCovs = exampleCovs, ncores = 3)
res_over
#> 
#> === One-Step-Ahead (OSA) Residuals ===
#> Total observations: 1500 
#> Number of tracks:   3 
#> 
#> ---- Goodness-of-Fit Tests ----
#>  metric   statistic      p.value
#>    KS_x  0.02500516 3.065233e-01
#>    KS_y  0.02854070 1.744104e-01
#>  KS_mah  0.01735706 7.579011e-01
#>    LB_x 33.14761549 2.485147e-05
#>    LB_y 65.63052339 1.123146e-11
#>  LB_mah  6.54028837 4.782583e-01
#> -------------------------------
#> 
#> Residual Summary:
#>    residual.x         residual.y      
#>  Min.   :-3.11105   Min.   :-3.55955  
#>  1st Qu.:-0.67883   1st Qu.:-0.72711  
#>  Median :-0.02908   Median :-0.05798  
#>  Mean   :-0.03422   Mean   :-0.04454  
#>  3rd Qu.: 0.63190   3rd Qu.: 0.63178  
#>  Max.   : 2.89999   Max.   : 2.94656  
#>  NA's   :3          NA's   :3

p_over <- plot(res_over)
p_over$qq_x + p_over$qq_y + p_over$acf_x + p_over$acf_y + plot_layout(ncol=2)
```

![](man/figures/README-osa-2.png)<!-- -->

#### Bhattacharyya’s affinity

``` r
# calculate similarity of true and estimated UDs using Bhattacharyya's affinity
rasterOverlap(exp(UD), exp(trueUD))
#> [1] 0.9088081
```

#### Regional presence probability

``` r
# create a spatial mask for the region of interest
d2c <- exampleCovs$d2c < 2.5

reg_prob <- regionProb(fit_under,
                       spatialCovs = exampleCovs, 
                       mask = d2c, # region of interest
                       nSims = 1000, # number of Monte Carlo simulations
                       show_progress = FALSE)

reg_prob
#> Regional Probability Estimate
#> =============================
#> Point Estimate: 0.3099
#> 
#> Delta Method Approximation:
#>   Standard Error: 0.1976
#>   95% CI:         [0.0000, 0.6972]
#> 
#> Monte Carlo Simulation:
#>   Standard Error: 0.2172
#>   95% CI:         [0.0000, 0.7124]
#>   (Based on 1000 draws)

plot(reg_prob, log = TRUE)
```

![](man/figures/README-regionProb-1.png)<!-- -->

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
