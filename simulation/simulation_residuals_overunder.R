library(langevinSSM)


n_sims <- 100
p_vals_ks_x <- p_vals_ks_y <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("underdamped", "overdamped")))
p_vals_lb_x <- p_vals_lb_y <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("underdamped", "overdamped")))
p_vals_ks_mah <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("underdamped", "overdamped")))
p_vals_lb_mah <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("underdamped", "overdamped")))
AIC_mat <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("underdamped", "overdamped")))

true_model <- "underdamped"
true_par <- list(beta=c(-4,6,5,-0.1), sigma = 5, gamma = 0.5)

p_val <- 0.05

set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
for(i in 1:n_sims) {

  sim_data <- simLangevin(model = true_model, par = true_par, obsPerAnimal = 5000, nbAnimals = 3, subSample = list(samplingRate = 10),
                         spatialCovs = exampleCovs, measurementError = list(smaj.sd=1.5,smin.sd=0.75))

  fit_under <- fitLangevin(data = sim_data,
                          spatialCovs = exampleCovs,
                          calcResiduals = TRUE, silent = TRUE)

  fit_over <- fitLangevin(data = sim_data, model = "overdamped",
                     spatialCovs = exampleCovs,
                     calcResiduals = TRUE, silent = TRUE)

  AIC_mat[i,] <- AIC(fit_under, fit_over)$AIC

  tests_under <- attr(fit_under$residuals, "tests")
  tests_over <- attr(fit_over$residuals, "tests")

  p_vals_ks_x[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "KS_x"]
  p_vals_ks_y[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "KS_y"]
  p_vals_ks_mah[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "KS_mah"]
  p_vals_lb_x[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "LB_x"]
  p_vals_lb_y[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "LB_y"]
  p_vals_lb_mah[i,"underdamped"] <- tests_under$p.value[tests_under$metric == "LB_mah"]

  p_vals_ks_x[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "KS_x"]
  p_vals_ks_y[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "KS_y"]
  p_vals_ks_mah[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "KS_mah"]
  p_vals_lb_x[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "LB_x"]
  p_vals_lb_y[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "LB_y"]
  p_vals_lb_mah[i,"overdamped"] <- tests_over$p.value[tests_over$metric == "LB_mah"]

  print(paste0("Simulation ", i, ": ",paste(apply(p_vals_ks_x[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for KS_x",colnames(p_vals_ks_x),collapse=";   ")))
  print(paste0("Simulation ", i, ": ",paste(apply(p_vals_ks_y[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for KS_y",colnames(p_vals_ks_y),collapse=";   ")))
  print(paste0("Simulation ", i, ": ",paste(apply(p_vals_ks_mah[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for KS_mah",colnames(p_vals_ks_mah),collapse=";   ")))
  print(paste0("Simulation ", i, ": ",paste(apply(p_vals_lb_x[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for LB_x",colnames(p_vals_lb_x),collapse=";   ")))
  print(paste0("Simulation ", i, ": ",paste(apply(p_vals_lb_y[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for LB_y",colnames(p_vals_lb_y),collapse=";   ")))
  print(paste0("Simulation ", i, ": ",paste(apply(p_vals_lb_mah[1:i,,drop=FALSE]< p_val, 2, sum), "rejections for LB_mah",colnames(p_vals_lb_mah),collapse=";   ")))

}

apply(p_vals_ks_x < p_val, 2, mean)
apply(p_vals_ks_y < p_val, 2, mean)
apply(p_vals_ks_mah < p_val, 2, mean)
apply(p_vals_lb_x < p_val,2, mean)
apply(p_vals_lb_y < p_val,2, mean)
apply(p_vals_lb_mah < p_val,2, mean)
