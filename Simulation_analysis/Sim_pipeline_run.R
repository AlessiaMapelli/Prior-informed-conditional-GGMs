# Simulation Pipeline for Covariate-Dependent Gaussian Graphical Models
# This script implements a comprehensive simulation study following the methodology
# described in the referenced paper
setwd("/group/diangelantonio/users/alessia_mapelli/Complete_projects/Prior-informed-conditional-GGMs/Simulation_analysis")
rm(list = ls(all = TRUE))

# Source required functions
source("Sim_functions.R")
source("/group/diangelantonio/users/alessia_mapelli/Complete_projects/Prior-informed-conditional-GGMs/Computational_templates/ggReg_main_functions.R")

#################################################
## SIMULATION CONFIGURATION
#################################################

# Define simulation grid
SIMULATION_GRID <- expand.grid(
  n_samples = c(500, 1000),
  n_nodes = c(50, 100),
  prior_type = c("perfect", "noisy", "none"),
  n_covariates = c(1, 3),  # 1=binary only, 3=2 continuous + 1 binary, 5=3 continuous (1 no effect) + 2 binary(1 no effect)
  symm_method = c("OR"), 
  stringsAsFactors = FALSE
)


# Number of replications per configuration
N_REPLICATIONS <- 10

# Fixed network parameters (following the paper)
NETWORK_CONFIG <- list(
  population_method = "SF",        # Scale-free for population network
  population_power_sf = 2.5,          # Power law exponent
  population_n_edges_sf = NULL,         # Will be computed based on density
  population_density_sf = 1,        # Target edge density for population network
  population_edge_prob_er = NULL,
  
  covariate_method = "ER",         # Erdős-Rényi for covariate effects
  covariate_power_sf = NULL,
  covariate_n_edges_sf = NULL,
  covariate_density_sf = NULL,
  covariate_edge_prob_er= 0.1,           # Edge probability for covariate networks
  
  
  coefficient_range = c(0.35, 0.5), # Magnitude range for coefficients
  sign_options = c(-1, 1)          # Allow both positive and negative edges
)

# Covariate configuration
base_covariate_config <- function(n_covariates) {
  if (n_covariates == 1) {
    list(n_continuous = 0, n_binary = 1, mean_sparsity = 0, effective_cov = 1, effective_mean=FALSE)
  } else if (n_covariates == 3) {
    list(n_continuous = 2, n_binary = 1, mean_sparsity = 0, effective_cov = 2, effective_mean=FALSE)
  } else if (n_covariates == 5) {
    list(n_continuous = 3, n_binary = 2, mean_sparsity = 0, effective_cov = 3, effective_mean=FALSE)
  }
}

# Prior knowledge noise parameters
PRIOR_NOISE_PARAMS <- list(
  sensitivity = 0.8,                    # Probability of detecting true edge
  specificity = 0.9,                    # Probability of correctly identifying non-edge
  confidence_range = c(0.7, 1)       # Range for edge confidence scores
)

get_method_params <- function(prior_type) {
  base_params <- list(
    use_slurm = TRUE,
    slurm_script_path = "./slurm_ggReg_node.sbatch",
    lambda_prec_type = "min",
    tune_hyperparams = TRUE,
    asparse_grid = c(0.5, 0.75, 0.9, 0.99),
    screening_procedure=FALSE,
    random_hyper_search = TRUE,
    p.rand.hyper = 0.5,
    K = 5,
    verbose = TRUE
  )
  
  # Adjust weight grid based on prior type
  if (prior_type == "none") {
    base_params$weight_grid = c(1.0)  # No prior weighting
  } else {
    base_params$weight_grid = c(0.8, 1.0, 1.1, 1.3)  # Include prior weighting
  }
  
  return(base_params)
}


TEST_SIMULATION_GRID <- expand.grid(
  n_samples = c(200),
  n_nodes = c(5),
  n_covariates = c(1),
  prior_type = c("perfect"),
  symm_method = c("OR"),
  stringsAsFactors = FALSE
)


#################################################
## SIMULATION EXECUTION FUNCTIONS
#################################################
simulation_output_folder = "Simulation_results"
dir.create(simulation_output_folder)
write.csv(SIMULATION_GRID, file=paste0(simulation_output_folder, "/input_computation_file.csv"))
simulation_output_file   = "simulation_results.RData"

#################################################
## SEQUENTIAL EXECUTION - Test simulation
#################################################

simulation_output_folder = "Simulation_results/test"
dir.create(simulation_output_folder)
write.csv(SIMULATION_GRID[16,], file=paste0(simulation_output_folder, "/input_computation_file.csv"))
simulation_output_file   = "simulation_results_test.RData"
source("Sim_functions.R")
all_results <- run_complete_simulation(SIMULATION_GRID = SIMULATION_GRID[16,],
                                       N_REPLICATIONS = 1,
                                       output_folder= simulation_output_folder,
                                       output_file = simulation_output_file)

#################################################
## SEQUENTIAL EXECUTION - FULL SIMULATION PIPELINE
#################################################

all_results <- run_complete_simulation(SIMULATION_GRID = SIMULATION_GRID,
                                    N_REPLICATIONS = N_REPLICATIONS,
                                    output_folder= simulation_output_folder,
                                    output_file = simulation_output_file)


#################################################
## PARALLEL EXECUTION - GENERATE INPUT DATASETS FOR ALL CONFIGURATIONS
#################################################

input_computation <- generate_input_datasets_simulation(
  SIMULATION_GRID = SIMULATION_GRID,
  N_REPLICATIONS = N_REPLICATIONS,
  output_folder= simulation_output_folder,
  output_file = simulation_output_file
  )

#################################################
## PARALLEL EXECUTION - ALGORITHM COMPUTATION IN PARALLEL ACROSS CONFIGURATIONS
#################################################
# Run from terminal
# cd /group/diangelantonio/users/alessia_mapelli/Complete_projects/Prior-informed-conditional-GGMs/Simulation_analysis
# Simulation_sequential.sh path_to_input_datasets

#################################################
## PARALLEL EXECUTION - EVALUATION OF THE OUTCOME
#################################################
all_results <- list()
for(i in 1:nrow(input_computation)){
  input_computation_row <- input_computation[i,]
  result <- collect_and_evaluate_resuts(p= input_computation_row$p, 
                                        output_path =input_computation_row$output_path,
                                        name_output =input_computation_row$name_output,
                                        symm_method = "OR")
  config <- result$config
  rep_id <- result
  result_id <- sprintf("n%d_p%d_q%d_%s_%s_rep%d", 
                       config$n_samples, config$n_nodes, config$n_covariates,
                       config$prior_type, config$symm_method, rep_id)
  all_results[[result_id]] <- result

  result <- collect_and_evaluate_resuts(p= input_computation_row$p, 
                                      output_path =input_computation_row$output_path,
                                      name_output =input_computation_row$name_output,
                                      symm_method = "AND")
  config <- result$config
  rep_id <- result
  result_id <- sprintf("n%d_p%d_q%d_%s_%s_rep%d", 
                       config$n_samples, config$n_nodes, config$n_covariates,
                       config$prior_type, config$symm_method, rep_id)
  all_results[[result_id]] <- result

  if (i %% N_REPLICATIONS == 0) {
    cat(paste0("Processing line ", i , " out of ", nrow(input_computation), "\n"))
    cat("Saving intermediate results...\n")
    save(all_results, SIMULATION_GRID, file = paste0(simulation_output_folder,"/", "temp_untill_row", i, "_", simulation_output_file))
  }
  }
cat("Saving final results...\n")
save(all_results, SIMULATION_GRID, file = paste(simulation_output_folder,simulation_output_file, sep="/"))


#################################################
## RESULTS ANALYSIS AND VISUALIZATION
#################################################

cat("\nProcessing and analyzing results...\n")

processed_data <- process_simulation_results(all_results)
#save(processed_data, file = paste(simulation_output_folder,simulation_output_file, sep="/"))

# Generate analysis
cat("Generating plots...\n")
plots <- generate_analysis_plots(processed_data)

library(forcats)
# 1. Prot performance by prior type
Prot_delta_data <- processed_data[processed_data$component == "Prot", ]
Prot_delta_data <- Prot_delta_data %>%
  mutate(prior_type = fct_relevel(prior_type, "none", "noisy", "perfect"))

# Aggregate by configuration
delta_summary_data <- Prot_delta_data %>%
    group_by(n_samples, n_nodes, n_covariates, prior_type, symm_method) %>%
    summarise(
      mean_TPR = mean(TPR, na.rm = TRUE),
      mean_FPR = mean(FPR, na.rm = TRUE),
      mean_F1 = mean(F1, na.rm = TRUE),
      mean_Magnitude_preserved= mean(Magnitude_preserved, na.rm = TRUE),
      sd_TPR = sd(TPR, na.rm = TRUE),
      sd_FPR = sd(FPR, na.rm = TRUE),
      sd_F1 = sd(F1, na.rm = TRUE),
      sd_Magnitude_preserved = sd(Magnitude_preserved, na.rm = TRUE),
      .groups = 'drop'
    )

new_plots <- list()

# TPR comparison for Delta 0
new_plots$delta_0_tpr_comparison <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_TPR, 
                                                               fill = prior_type, color = prior_type)) +
    geom_boxplot(position = position_dodge(0.8)) +
    geom_errorbar(aes(ymin = pmax(0, mean_TPR - sd_TPR), 
                      ymax = pmin(1, mean_TPR + sd_TPR), 
                      colour = prior_type),
                  position = position_dodge(0.8), width = 0.2) +
    facet_grid(n_nodes ~ n_covariates, labeller = label_both)+
    labs(title = "Performance in network reconstruction: True Positive Rate by Prior Type",
         subtitle = "Perfromance metrics measured on the baseline network for which the prior is provided",
         x = "Sample Size", y = "True Positive Rate",
         fill = "Prior Type") +
    theme_minimal() +
    guides(colour = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

new_plots$delta_0_fpr_comparison <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_FPR, 
                                                                 fill = prior_type, color = prior_type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = pmax(0, mean_FPR - sd_FPR), 
                    ymax = pmin(1, mean_FPR + sd_FPR), 
                    colour = prior_type),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(n_nodes ~ n_covariates, labeller = label_both) +
  labs(title = "Performance in network reconstruction: Frue Positive Rate by Prior Type",
       subtitle = "Perfromance metrics measured on the baseline network for which the prior is provided",
       x = "Sample Size", y = "Frue Positive Rate",
       fill = "Prior Type") +
  guides(colour = "none")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

  
# F1 Score comparison for Delta
new_plots$delta_0_f1_comparison <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_F1, 
                                                              fill = prior_type, color = prior_type)) +
    geom_boxplot(position = position_dodge(0.8)) +
    geom_errorbar(aes(ymin = pmax(0, mean_F1 - sd_F1), 
                      ymax = pmin(1, mean_F1 + sd_F1),
                      color = prior_type),
                  position = position_dodge(0.8), width = 0.2) +
    facet_grid(n_nodes ~ n_covariates, labeller = label_both) +
    labs(title = "Performance in network reconstruction: F1 Score by Prior Type",
         subtitle = "Perfromance metrics measured on the baseline network for which the prior is provided",
         x = "Sample Size", y = "F1 Score",
         fill = "Prior Type") +
    guides(colour = "none")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

new_plots$delta_0_error_comparison <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_Magnitude_preserved, 
                                                               fill = prior_type, color= prior_type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = pmax(0, mean_Magnitude_preserved - sd_Magnitude_preserved), 
                    ymax = mean_Magnitude_preserved + sd_Magnitude_preserved,
                    color= prior_type),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(n_nodes ~ n_covariates, labeller = label_both, scales = "free_y") +
  labs(title = "Performance in network reconstruction: Correlalation by Prior Type",
       subtitle = "Perfromance metrics measured on the baseline network for which the prior is provided",
       x = "Sample Size", y = "Spearman correlation",
       fill = "Prior Type") +
  guides(colour = "none")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


delta_component_data <- processed_data[processed_data$component_type == "delta_individual" & processed_data$n_covariates==1, ]
delta_component_data <- delta_component_data %>%
  mutate(prior_type = fct_relevel(prior_type, "none", "noisy", "perfect"))


delta_comp_summary <- delta_component_data %>%
    group_by(component,n_samples, prior_type, n_nodes) %>%
    summarise(
      mean_TPR = mean(TPR, na.rm = TRUE),
      mean_F1 = mean(F1, na.rm = TRUE),
      mean_Magnitude_preserved = mean(Magnitude_preserved, na.rm = TRUE),
      mean_accuracy = mean(Accuracy, na.rm = TRUE),
      sd_TPR = sd(TPR, na.rm = TRUE),
      sd_F1 = sd(F1, na.rm = TRUE),
      sd_Magnitude_preserved = sd(Magnitude_preserved, na.rm = TRUE),
      sd_accuracy = sd(Accuracy, na.rm = TRUE),
      .groups = 'drop'
    )
  
new_plots$delta_0_1_component_tpr <- ggplot(delta_comp_summary, aes(x = component, y = mean_TPR, 
                                                              fill = prior_type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = pmax(0, mean_TPR - sd_TPR), 
                      ymax = mean_TPR + sd_TPR),
                  position = position_dodge(0.8), width = 0.2) +
  facet_grid(n_nodes ~ n_samples, labeller = label_both) +
    labs(title = "Performance in network reconstruction: TPR by Component and Prior Type",
         subtitle = "Perfromance metrics measured on the baseline network and the single binary covariate",
         x = "Component", y = "Mean TPR",
         fill = "Prior Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
new_plots$delta_0_1_component_F1 <- ggplot(delta_comp_summary, aes(x = component, y = mean_F1, 
                                                             fill = prior_type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = pmax(0, mean_F1 - sd_F1), 
                      ymax = mean_F1 + sd_F1),
                  position = position_dodge(0.8), width = 0.2) +
  facet_grid(n_nodes ~ n_samples, labeller = label_both) +
    labs(title = "Performance in network reconstruction: F1 by Component and Prior Type",
         subtitle = "Perfromance metrics measured on the baseline network and the single binary covariate",
         x = "Component", y = "Mean F1",
         fill = "Prior Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
new_plots$delta_0_1_component_Acc <- ggplot(delta_comp_summary, aes(x = component, y = mean_accuracy, 
                                                              fill = prior_type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = pmax(0, mean_accuracy - sd_accuracy), 
                      ymax = mean_accuracy + sd_accuracy),
                  position = position_dodge(0.8), width = 0.2) +
  facet_grid(n_nodes ~ n_samples, labeller = label_both) +
    labs(title = "Performance in network reconstruction: Accuracy by Component and Prior Type",
         subtitle = "Perfromance metrics measured on the baseline network and the single binary covariate",
         x = "Component", y = "Mean accuracy",
         fill = "Prior Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
new_plots$delta_0_1_component_error <- ggplot(delta_comp_summary, aes(x = component, y = mean_Magnitude_preserved, 
                                                                fill = prior_type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = pmax(0, mean_Magnitude_preserved - sd_Magnitude_preserved), 
                      ymax = mean_Magnitude_preserved + sd_Magnitude_preserved),
                  position = position_dodge(0.8), width = 0.2) +
    facet_grid(n_nodes ~ n_samples, labeller = label_both) +
    labs(title = "Performance in network reconstruction: Correlation by Component and Prior Type",
         subtitle = "Perfromance metrics measured on the baseline network and the single binary covariate",
         x = "Component", y = "Spearman correlation",
         fill = "Prior Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

delta_component_data <- processed_data[processed_data$TP + processed_data$FN == 0 & processed_data$component_type=="delta_individual", ]
delta_component_data <- delta_component_data %>%
  mutate(prior_type = fct_relevel(prior_type, "none", "noisy", "perfect"))


delta_comp_summary <- delta_component_data %>%
  group_by(n_samples, prior_type, n_nodes) %>%
  summarise(
    mean_FPR = mean(FPR, na.rm = TRUE),
    mean_accuracy = mean(Accuracy, na.rm = TRUE),
    sd_FPR = sd(FPR, na.rm = TRUE),
    sd_accuracy = sd(Accuracy, na.rm = TRUE),
    .groups = 'drop'
  )

new_plots$delta_null_component_Acc <- ggplot(delta_comp_summary, aes(x = prior_type, y = mean_fpr, 
                                                                    fill = prior_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_accuracy - sd_accuracy), 
                    ymax = mean_accuracy + sd_accuracy),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(n_nodes ~ n_samples, labeller = label_both) +
  labs(title = "Performance in network reconstruction: Accuracy by Prior Type",
       subtitle = "Perfromance metrics measured on the network with zero prec matrix",
       x = "Component", y = "Mean accuracy",
       fill = "Prior Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

new_plots$delta_null_component_fpr <- ggplot(delta_comp_summary, aes(x = prior_type, y = mean_FPR, 
                                                                     fill = prior_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_FPR - sd_FPR), 
                    ymax = mean_FPR + sd_FPR),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(n_nodes ~ n_samples, labeller = label_both) +
  labs(title = "Performance in network reconstruction: FPR by Prior Type",
       subtitle = "Perfromance metrics measured on the network with zero prec matrix",
       x = "Component", y = "Mean FPR",
       fill = "Prior Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



#################################################
## PREPARE PLOTS FOR THE PRESENTATION
#################################################
library(dplyr)
library(ggplot2)
library(forcats)

Prot_delta_data <- processed_data[processed_data$component == "Prot" & processed_data$n_samples==5000 & processed_data$n_nodes==100, ]


Prot_delta_data <- Prot_delta_data %>%
  mutate(prior_type = fct_relevel(prior_type, "none", "noisy", "perfect"))

prior_cols <- c(
  "none"    = "#B7C9E2",
  "noisy"   = "#9BC2F9",   # orange
  "perfect" = "#6387FD"    # blue
)

# Aggregate by configuration
delta_summary_data <- Prot_delta_data %>%
  group_by(n_samples, n_nodes, n_covariates, prior_type, symm_method) %>%
  summarise(
    mean_TPR = mean(TPR, na.rm = TRUE),
    mean_FPR = mean(FPR, na.rm = TRUE),
    mean_F1 = mean(F1, na.rm = TRUE),
    mean_Magnitude_preserved= mean(Magnitude_preserved, na.rm = TRUE),
    sd_TPR = sd(TPR, na.rm = TRUE),
    sd_FPR = sd(FPR, na.rm = TRUE),
    sd_F1 = sd(F1, na.rm = TRUE),
    sd_Magnitude_preserved = sd(Magnitude_preserved, na.rm = TRUE),
    .groups = 'drop'
  )


ggplot(delta_summary_data, aes(x = prior_type, y = mean_F1, 
                               fill = prior_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_F1 - sd_F1), 
                    ymax = pmin(1, mean_F1 + sd_F1)),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(. ~ n_covariates, labeller = label_both) +
  scale_fill_manual(values = prior_cols, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Performance in baseline network reconstruction",
       subtitle = "n_nodes = 100 nodes and n_samples= 1000",
       x = "Prior Type", y = "F1 Score") +
  guides(colour = "none")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

delta_rep <- delta_summary_data %>%
  filter(n_covariates == 3)


ggplot(delta_rep, aes(x = prior_type, y = mean_F1, 
                               fill = prior_type)) +
  geom_col(width = 0.5) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_F1 - sd_F1), 
                    ymax = pmin(1, mean_F1 + sd_F1)),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = prior_cols, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Performance in baseline network reconstruction",
       subtitle = "(n_nodes = 100 nodes, n_samples= 1000, n_covariates=3)",
       x = "\n Prior Type", y = "F1 Score") +
  guides(colour = "none")+
  theme_minimal() +
  theme(axis.text.x = element_text(size=11, hjust = 0.5),
        plot.title = element_text(size=14,hjust = 0.5),
        plot.subtitle = element_text(size=13,hjust = 0.5),
        axis.text.y = element_text(size=11),
        axis.title =element_text(size=12) )

analysis_file <- file.path(simulation_output_folder, "simulation_analysis_no_weight_tuning_denser_x1_sparser_base.RData")
save(processed_data, plots, new_plots, file = analysis_file)

#################################################
## PREPARE TABLE FOR THE PAPER - TO DO
#################################################


cat("Generating summary tables...\n")
tables <- generate_summary_tables(processed_data)




build_param_dictionary <- function() {
  # Summarize the filtered simulation grid
  grid_summary <- list(
    n_rows = nrow(SIMULATION_GRID),
    n_samples = sort(unique(SIMULATION_GRID$n_samples)),
    n_nodes = sort(unique(SIMULATION_GRID$n_nodes)),
    n_covariates = sort(unique(SIMULATION_GRID$n_covariates)),
    prior_type = sort(unique(SIMULATION_GRID$prior_type)),
    symm_method = sort(unique(SIMULATION_GRID$symm_method)),
    filter_rules = c(
      #MANUALLY CHANGE IF NECESSARY
      "Removed: n_nodes >= 200 & n_samples <= 400",
      "Removed: n_nodes >= 100 & n_samples <= 200"
    )
  )
  
  # Covariate config rules (human-readable)
  #MANUALLY CHANGE IF NECESSARY
  covariate_config_rules <- list(
    `1` = list(n_continuous = 0, n_binary = 1, mean_sparsity = 0.3),
    `3` = list(n_continuous = 1, n_binary = 2, mean_sparsity = 0.2),
    `5` = list(n_continuous = 2, n_binary = 3, mean_sparsity = 0.15)
  )
  
  # Method params for each prior type as actually used
  method_params_by_prior <- lapply(sort(unique(SIMULATION_GRID$prior_type)), function(pt) {
    params <- get_method_params(pt)
    params
  })
  names(method_params_by_prior) <- sort(unique(SIMULATION_GRID$prior_type))
  
  # Execution subset info (if you used SIMULATION_GRID_TEST)
  exec_subset <- tryCatch({
    list(
      n_rows = nrow(SIMULATION_GRID_TEST),
      grid = SIMULATION_GRID_TEST
    )
  }, error = function(e) NULL)
  
  # Output paths (from your script)
  output <- list(
    simulation_output_folder = simulation_output_folder,
    simulation_output_file   = simulation_output_file,
    analysis_file            = analysis_file
  )
  
  list(
    created_at = as.character(Sys.time()),
    N_REPLICATIONS = N_REPLICATIONS,
    grid = grid_summary,
    network_config = NETWORK_CONFIG,
    prior_noise_params = PRIOR_NOISE_PARAMS,
    covariate_config_rules = covariate_config_rules,
    method_params_by_prior = method_params_by_prior,
    execution_subset = exec_subset,
    output_paths = output,
    session_info = capture.output(sessionInfo())
  )
}

param_dictionary <- build_param_dictionary()

# Save a human-readable Markdown summary (base R only)
md_lines <- c(
  "# Simulation Parameters",
  paste0("- **Created at:** ", param_dictionary$created_at),
  paste0("- **N_REPLICATIONS:** ", param_dictionary$N_REPLICATIONS),
  "",
  "## Grid (after filtering)",
  paste0("- **Rows:** ", param_dictionary$grid$n_rows),
  paste0("- **n_samples:** ", paste(param_dictionary$grid$n_samples, collapse = ", ")),
  paste0("- **n_nodes:** ", paste(param_dictionary$grid$n_nodes, collapse = ", ")),
  paste0("- **n_covariates:** ", paste(param_dictionary$grid$n_covariates, collapse = ", ")),
  paste0("- **prior_type:** ", paste(param_dictionary$grid$prior_type, collapse = ", ")),
  paste0("- **symm_method:** ", paste(param_dictionary$grid$symm_method, collapse = ", ")),
  "- **Filter rules:**",
  paste0("  - ", param_dictionary$grid$filter_rules),
  "",
  "## Network config",
  paste0("- **population_method:** ", NETWORK_CONFIG$population_method),
  paste0("- **population_power:** ", NETWORK_CONFIG$population_power),
  paste0("- **population_density:** ", NETWORK_CONFIG$population_density),
  paste0("- **covariate_method:** ", NETWORK_CONFIG$covariate_method),
  paste0("- **covariate_prob:** ", NETWORK_CONFIG$covariate_prob),
  paste0("- **coefficient_range:** [", paste(NETWORK_CONFIG$coefficient_range, collapse = ", "), "]"),
  paste0("- **sign_options:** [", paste(NETWORK_CONFIG$sign_options, collapse = ", "), "]"),
  "",
  "## Prior noise params",
  paste0("- **sensitivity:** ", PRIOR_NOISE_PARAMS$sensitivity),
  paste0("- **specificity:** ", PRIOR_NOISE_PARAMS$specificity),
  paste0("- **confidence_range:** [", paste(PRIOR_NOISE_PARAMS$confidence_range, collapse = ", "), "]"),
  "",
  "## Method params by prior",
  paste(
    vapply(names(param_dictionary$method_params_by_prior), function(nm) {
      p <- param_dictionary$method_params_by_prior[[nm]]
      paste0("- **", nm, ":** ",
             "use_slurm=", p$use_slurm, "; ",
             "tune_hyperparams=", p$tune_hyperparams, "; ",
             "asparse_grid=[", paste(p$asparse_grid, collapse = ", "), "]; ",
             "screening_procedure=", p$screening_procedure, "; ",
             "K=", p$K, "; ",
             "verbose=", p$verbose, "; ",
             "weight_grid=[", paste(p$weight_grid, collapse = ", "), "]")
    }, character(1L)),
    collapse = "\n"
  ),
  "",
  "## Execution subset (if used)",
  paste0("- **Rows:** ", if (!is.null(param_dictionary$execution_subset)) param_dictionary$execution_subset$n_rows else 0),
  paste0("- **n_samples:** ", if (!is.null(param_dictionary$execution_subset)) paste(sort(unique(param_dictionary$execution_subset$grid$n_samples)), collapse = ", ") else 0),
  paste0("- **n_nodes:** ", if (!is.null(param_dictionary$execution_subset)) paste(sort(unique(param_dictionary$execution_subset$grid$n_nodes)), collapse = ", ") else 0),
  paste0("- **n_covariates:** ", if (!is.null(param_dictionary$execution_subset)) paste(sort(unique(param_dictionary$execution_subset$grid$n_covariates)), collapse = ", ") else 0),
  paste0("- **prior_type:** ", if (!is.null(param_dictionary$execution_subset)) paste(sort(unique(param_dictionary$execution_subset$grid$prior_type)), collapse = ", ") else 0),
  paste0("- **symm_method:** ", if (!is.null(param_dictionary$execution_subset)) paste(sort(unique(param_dictionary$execution_subset$grid$symm_method)), collapse = ", ") else 0),
  "",
  "## Output paths",
  paste0("- **simulation_output_folder:** ", param_dictionary$output_paths$simulation_output_folder),
  paste0("- **simulation_output_file:** ", param_dictionary$output_paths$simulation_output_file),
  paste0("- **analysis_file:** ", param_dictionary$output_paths$analysis_file)
)
dict_md_path <- file.path(simulation_output_folder, "simulation_parameters_no_weight_tuning_denser_x1_sparser_base.md")
writeLines(md_lines, dict_md_path)
