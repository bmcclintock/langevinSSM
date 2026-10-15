library(langevinSSM)

n_sims <- 100
nbAnimals <- 3
obsPerAnimal <- 5000

obs_dist <- "Gaussian"
if(obs_dist == "Gaussian"){
  osaMethod <- "oneStepGaussianOffMode"
} else osaMethod <- "oneStepGeneric"

p_vals_ks_x <- p_vals_ks_y <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("underdamped", "overdamped")))
p_vals_lb_x <- p_vals_lb_y <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("underdamped", "overdamped")))
if(obs_dist=="Gaussian") p_vals_ks_mah <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("underdamped", "overdamped")))
if(obs_dist=="Gaussian") p_vals_lb_mah <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("underdamped", "overdamped")))
AIC_mat <- cAIC_mat <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("underdamped", "overdamped")))

true_model <- "underdamped"
true_par <- list(beta=c(-4,6,5,-0.1), sigma = 5, gamma = 0.5)

p_val <- 0.05

set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
for(i in 1:n_sims) {

  sim_data <- simLangevin(model = true_model, obs_dist = obs_dist, par = true_par, obsPerAnimal = obsPerAnimal, nbAnimals = nbAnimals, subSample = list(samplingRate = 10),
                         spatialCovs = exampleCovs, measurementError = list(smaj.sd=1.5,smin.sd=0.75))

  fit_under <- suppressMessages(fitLangevin(data = sim_data, obs_dist = obs_dist,
                          spatialCovs = exampleCovs, silent = TRUE, getJointPrecision = TRUE))

  fit_under$residuals <- residuals(fit_under, data = sim_data, spatialCovs = exampleCovs, ncores = nbAnimals, method = osaMethod, reverse=TRUE, splineApprox=TRUE)

  fit_over <- suppressMessages(fitLangevin(data = sim_data, obs_dist = obs_dist, model = "overdamped",
                     spatialCovs = exampleCovs, silent = TRUE, getJointPrecision = TRUE))

  fit_over$residuals <- residuals(fit_over, data = sim_data, spatialCovs = exampleCovs, ncores = nbAnimals, method = osaMethod, reverse=TRUE, splineApprox=TRUE)


  AIC_mat[i,] <- AIC(fit_under, fit_over)$AIC # not valid for comparing underdamped vs overdamped

  cAIC_mat[i,1] <- cAIC(fit_under,data=sim_data,spatialCovs=exampleCovs)$cAIC
  cAIC_mat[i,2] <- cAIC(fit_over,data=sim_data,spatialCovs=exampleCovs)$cAIC

  tests_under <- attr(fit_under$residuals, "tests")
  tests_over <- attr(fit_over$residuals, "tests")

  p_vals_ks_x[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "KS_x"]
  p_vals_ks_y[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "KS_y"]
  if(obs_dist=="Gaussian") p_vals_ks_mah[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "KS_mah"]
  p_vals_lb_x[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "LB_x"]
  p_vals_lb_y[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "LB_y"]
  if(obs_dist=="Gaussian") p_vals_lb_mah[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "LB_mah"]

  p_vals_ks_x[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "KS_x"]
  p_vals_ks_y[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "KS_y"]
  if(obs_dist=="Gaussian") p_vals_ks_mah[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "KS_mah"]
  p_vals_lb_x[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "LB_x"]
  p_vals_lb_y[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "LB_y"]
  if(obs_dist=="Gaussian") p_vals_lb_mah[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "LB_mah"]

  print(paste0("Simulation ", i, ": ",paste(apply(p_vals_ks_x[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for KS_x",colnames(p_vals_ks_x),collapse=";   ")))
  print(paste0("Simulation ", i, ": ",paste(apply(p_vals_ks_y[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for KS_y",colnames(p_vals_ks_y),collapse=";   ")))
  if(obs_dist=="Gaussian") print(paste0("Simulation ", i, ": ",paste(apply(p_vals_ks_mah[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for KS_mah",colnames(p_vals_ks_mah),collapse=";   ")))
  print(paste0("Simulation ", i, ": ",paste(apply(p_vals_lb_x[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for LB_x",colnames(p_vals_lb_x),collapse=";   ")))
  print(paste0("Simulation ", i, ": ",paste(apply(p_vals_lb_y[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for LB_y",colnames(p_vals_lb_y),collapse=";   ")))
  if(obs_dist=="Gaussian") print(paste0("Simulation ", i, ": ",paste(apply(p_vals_lb_mah[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for LB_mah",colnames(p_vals_lb_mah),collapse=";   ")))

}

apply(p_vals_ks_x < p_val, 2, mean)
apply(p_vals_ks_y < p_val, 2, mean)
if(obs_dist=="Gaussian") apply(p_vals_ks_mah < p_val, 2, mean)
apply(p_vals_lb_x < p_val,2, mean)
apply(p_vals_lb_y < p_val,2, mean)
if(obs_dist=="Gaussian") apply(p_vals_lb_mah < p_val,2, mean)
