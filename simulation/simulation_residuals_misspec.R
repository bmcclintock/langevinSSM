library(langevinSSM)


n_sims <- 100
p_vals_ks_x <- p_vals_ks_y <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("correct", "misspec")))
p_vals_lb_x <- p_vals_lb_y <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("correct", "misspec")))
p_vals_ks_mah <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("correct", "misspec")))
p_vals_lb_mah <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("correct", "misspec")))
AIC_mat <- matrix(0,n_sims,2,dimnames=list(1:n_sims, c("correct", "misspec")))

true_model <- "underdamped"
true_par <- list(beta=c(-4,6,5,-0.1), sigma = 5, gamma = 0.5)

p_val <- 0.05

set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
for(i in 1:n_sims) {

  sim_data <- simLangevin(model = true_model, par = true_par, obsPerAnimal = 5000, nbAnimals = 3, subSample = list(samplingRate = 10),
                          spatialCovs = exampleCovs, measurementError = list(smaj.sd=1.5,smin.sd=0.75))

  fit <- fitLangevin(data = sim_data, model = true_model,
                           spatialCovs = exampleCovs,
                           calcResiduals = TRUE, silent = TRUE)

  misCovs <- exampleCovs[1:3]
  for(j in 1:3){
    misCovs[[j]] <- terra::shift(simCov(sca = 200), dx=1000, dy=1000)
    names(misCovs[[j]]) <- paste0("cov",j)
  }
  misCovs$d2c <- exampleCovs$d2c

  misfit <- fitLangevin(data = sim_data,
                               spatialCovs = misCovs,
                               calcResiduals = TRUE, silent = TRUE)

  AIC_mat[i,] <- AIC(fit, misfit)$AIC

  tests_under <- attr(fit$residuals, "tests")
  tests_over <- attr(misfit$residuals, "tests")

  p_vals_ks_x[i,"correct"] <- tests_under$p.value[tests_under$metric == "KS_x"]
  p_vals_ks_y[i,"correct"] <- tests_under$p.value[tests_under$metric == "KS_y"]
  p_vals_ks_mah[i,"correct"] <- tests_under$p.value[tests_under$metric == "KS_mah"]
  p_vals_lb_x[i,"correct"] <- tests_under$p.value[tests_under$metric == "LB_x"]
  p_vals_lb_y[i,"correct"] <- tests_under$p.value[tests_under$metric == "LB_y"]
  p_vals_lb_mah[i,"correct"] <- tests_under$p.value[tests_under$metric == "LB_mah"]

  p_vals_ks_x[i,"misspec"] <- tests_over$p.value[tests_over$metric == "KS_x"]
  p_vals_ks_y[i,"misspec"] <- tests_over$p.value[tests_over$metric == "KS_y"]
  p_vals_ks_mah[i,"misspec"] <- tests_over$p.value[tests_over$metric == "KS_mah"]
  p_vals_lb_x[i,"misspec"] <- tests_over$p.value[tests_over$metric == "LB_x"]
  p_vals_lb_y[i,"misspec"] <- tests_over$p.value[tests_over$metric == "LB_y"]
  p_vals_lb_mah[i,"misspec"] <- tests_over$p.value[tests_over$metric == "LB_mah"]

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
