setwd("/group/diangelantonio/users/alessia_mapelli/Complete_projects/Prior-informed-conditional-GGMs/Computational_templates")
rm(list=ls(all=TRUE))
source("ggReg_main_functions.R")

load("/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/input_data_nodes.rda")
str(covariates)
covariates$x3 <- as.factor(covariates$x3)

########################################
## 1. Check the precision matrix estimation implementation
########################################
res <- GGReg_cov_estimation_sequential(
  Z0,
  known_ppi = known_ppi,
  covariates = covariates,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_sequential/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
plot_personalized_network(res$Dic_adj_matrics$Baseline)

res_2 <- GGReg_cov_estimation_sequential(
  Z0,
  known_ppi = known_ppi,
  covariates = covariates,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = FALSE,
  asparse_grid = 0.75,
  weight_grid = 1.0,
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_sequential/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
plot_personalized_network(res_2$Dic_adj_matrics$Baseline)


res_3 <- GGReg_cov_estimation_sequential(
  Z0,
  known_ppi = known_ppi,
  covariates = covariates,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = TRUE,
  p.rand.hyper = 0.1,
  K = 5,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_sequential/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
plot_personalized_network(res_3$Dic_adj_matrics$Baseline)


res_4 <- GGReg_cov_estimation_sequential(
  Z0,
  known_ppi = known_ppi,
  covariates = covariates,
  scr = TRUE,
  gamma = 0.1,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_sequential/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
plot_personalized_network(res_4$Dic_adj_matrics$Baseline)

res_5 <- GGReg_cov_estimation_sequential(
  Z0,
  known_ppi = known_ppi,
  covariates = NULL,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = FALSE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "./results_tests/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
plot_personalized_network(res_5$Dic_adj_matrics$Baseline)


res_6 <- GGReg_cov_estimation_sequential(
  Z0,
  known_ppi = NULL,
  covariates = covariates,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = FALSE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_sequential/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
plot_personalized_network(res_6$Dic_adj_matrics$Baseline)


res_7 <- GGReg_cov_estimation_sequential(
  Z0,
  known_ppi = NULL,
  covariates = NULL,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = FALSE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_sequential/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
plot_personalized_network(res_7$Dic_adj_matrics$Baseline)

res_8 <- GGReg_cov_estimation_sequential(
  Z0,
  known_ppi = known_ppi,
  covariates = covariates,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "1se",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_sequential/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
plot_personalized_network(res$Dic_adj_matrics$Baseline)
plot_personalized_network(res_8$Dic_adj_matrics$Baseline)

res_9 <- GGReg_cov_estimation_sequential(
  Z0,
  known_ppi = known_ppi,
  covariates = covariates,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_sequential/",
  name_output = "ggReg_result",
  symm_method ="AND",
  verbose = TRUE) 
plot_personalized_network(res$Dic_adj_matrics$Baseline)
plot_personalized_network(res_9$Dic_adj_matrics$Baseline)


pers_netw <- predict_personalized_network(
  res$Dic_adj_matrics, 
  training_covariates = covariates, 
  new_subject_covariates = data.frame(covariates[1,]),
  return_scaling_info = FALSE,
  verbose = TRUE) 

plot_personalized_network(pers_netw$personalized_network$Subject_1)



res_par <- GGReg_cov_estimation_parallel(
  Z0,
  known_ppi = known_ppi,
  covariates = covariates,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_par/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
plot_personalized_network(res_par$Dic_adj_matrics$Baseline)
plot_personalized_network(res$Dic_adj_matrics$Baseline)

res_slurm <- GGReg_cov_estimation_SLURM(
  Z0,
  known_ppi = known_ppi,
  covariates = covariates,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  slurm_script_path="./slurm_ggReg_node_v2.sbatch",
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_slurm/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
plot_personalized_network(res_slurm$Dic_adj_matrics$Baseline)
plot_personalized_network(res$Dic_adj_matrics$Baseline)

########################################
## 1. Check full implementation
########################################

source("ggReg_main_functions.R")
res <- GGReg_full_estimation(
  Z0,
  known_ppi = known_ppi,
  covariates = covariates,
  scr = FALSE,
  gamma = NULL,
  mean_estimation = FALSE,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  use_slurm = FALSE,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_sequential/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE) 
res$additional_info$dummy_params
new_subject <- data.frame(x1 =c(0.723269,0.22451) , x2 = c(0.22451,0.723269), x3 = c(1,0))
new_subject$x3 <- as.factor(new_subject$x3)
str(new_subject)
new_pers_networks <- predict_personalized_network(
  Dic_Delta_hat= res$results$Dic_adj_matrics,
  new_subject_covariates=new_subject,
  scaling_params=res$additional_info$scaling_params,
  dummy_params= res$additional_info$dummy_params,
  verbose = TRUE
)
new_pers_networks
compare_networks(new_pers_networks$Subject_1, new_pers_networks$Subject_2, method = "correlation")
plot_network_difference(new_pers_networks$Subject_1, new_pers_networks$Subject_2)

########################################
## 1. Usage example
########################################

set.seed(123)
n <- 1000  # number of samples
p <- 3   # number of features (nodes)
q <- 2    # number of covariates
# Generate synthetic expression data
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("Gene_", 1:p)
# Generate synthetic covariates
covariates <- data.frame(
  age = rnorm(n, mean = 50, sd = 10),
  sex = rbinom(n, 1, 0.5)
)
covariates$sex <- as.factor(ifelse(covariates$sex,"Male", "Female"))

# Optional: Prior knowledge network (PPI)
known_ppi <- matrix(0, p, p)
known_ppi[1:2, 1:2] <- 0.3
diag(known_ppi) <- 0

source("ggReg_main_functions.R")

results <- GGReg_full_estimation(
  x = X,
  known_ppi = known_ppi,
  covariates = covariates)

new_subject <- data.frame(age = 45, sex = as.factor("Female"))
pred_net <- predict_personalized_network(
  Dic_Delta_hat = results$results$Dic_adj_matrics,
  new_subject_covariates = new_subject,
  scaling_params = results$additional_info$scaling_params,
  dummy_params = results$additional_info$dummy_params
)

new_subjects <- data.frame(
 age = c(45, 30, 60),
 sex = as.factor(c("Female", "Male", "Female"))
)
pred_nets <- predict_personalized_network(
 Dic_Delta_hat = results$results$Dic_adj_matrics,
 new_subject_covariates = new_subjects,
 scaling_params = results$additional_info$scaling_params,
 dummy_params = results$additional_info$dummy_params
)

