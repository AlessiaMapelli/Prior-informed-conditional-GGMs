setwd("/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Code_git/4.Method_implementation")
rm(list=ls(all=TRUE))
source("ggReg_main_functions.R")

load("/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/input_data_nodes.rda")

########################################
## 1. Check the sequential implementation
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
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Example_data/n200_p10_q3_noisy_OR/rep1/results_sequential/",
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
