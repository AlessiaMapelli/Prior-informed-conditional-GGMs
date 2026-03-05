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
  covariate_edge_prob_er= 0.05,           # Edge probability for covariate networks
  
  
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
    use_slurm = FALSE,
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
  n_samples = c(500),
  n_nodes = c(50),
  n_covariates = c(1),
  prior_type = c("perfect", "noisy", "none"),
  symm_method = c("OR"),
  stringsAsFactors = FALSE
)


#################################################
## SEQUENTIAL EXECUTION - Test simulation
#################################################

simulation_output_folder = "Simulation_results/test"
dir.create(simulation_output_folder)
write.csv(TEST_SIMULATION_GRID, file=paste0(simulation_output_folder, "/input_computation_file.csv"))
simulation_output_file   = "simulation_results_test.RData"
all_results <- run_complete_simulation(SIMULATION_GRID = TEST_SIMULATION_GRID,
                                       N_REPLICATIONS = N_REPLICATIONS,
                                       output_folder= simulation_output_folder,
                                       output_file = simulation_output_file)

#################################################
## SIMULATION EXECUTION FUNCTIONS
#################################################
simulation_output_folder = "Simulation_results"
dir.create(simulation_output_folder)
write.csv(SIMULATION_GRID, file=paste0(simulation_output_folder, "/simulation_grid.csv"))
simulation_output_file   = "simulation_results.RData"


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
## PARALLEL EXECUTION - CHECK RESULTS PRESENCE
#################################################
#input_computation <- read.csv("Simulation_results/input_computation_file.csv")
for(i in 1:nrow(input_computation)){
  input_computation_row <- input_computation[i,]
  output_path =input_computation_row$output_path
  name_output =input_computation_row$name_output
  first_result_file <- paste0(output_path, name_output, "_node_1.rda")
  if (!file.exists(first_result_file)) {
    cat("Result file for node 1 not found for row ", i, " \n")
  }
}

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
  rep_id <- result$rep_id
  result_id <- sprintf("n%d_p%d_q%d_%s_%s_rep%d", 
                       config$n_samples, config$n_nodes, config$n_covariates,
                       config$prior_type, "OR", rep_id)
  all_results[[result_id]] <- result

  result <- collect_and_evaluate_resuts(p= input_computation_row$p, 
                                      output_path =input_computation_row$output_path,
                                      name_output =input_computation_row$name_output,
                                      symm_method = "AND")
  config <- result$config
  rep_id <- result$rep_id
  result_id <- sprintf("n%d_p%d_q%d_%s_%s_rep%d", 
                       config$n_samples, config$n_nodes, config$n_covariates,
                       config$prior_type, "AND", rep_id)
  all_results[[result_id]] <- result

  if (i %% N_REPLICATIONS == 0) {
    cat(paste0("Processing line ", i , " out of ", nrow(input_computation), "\n"))
    # cat("Saving intermediate results...\n")
    # save(all_results, SIMULATION_GRID, file = paste0(simulation_output_folder,"/", "temp_", simulation_output_file))
  }
  }
cat("Saving final results...\n")
save(all_results, SIMULATION_GRID, file = paste(simulation_output_folder,simulation_output_file, sep="/"))



#################################################
## RESULTS ANALYSIS AND VISUALIZATION
#################################################

cat("\nProcessing and analyzing results...\n")
#load(paste(simulation_output_folder,simulation_output_file, sep="/"))
processed_data <- process_simulation_results(all_results)
save(processed_data, file = paste0(simulation_output_folder, "/", "simulation_results_processed.RData"))

# Generate analysis
cat("Generating plots...\n")
plots <- generate_analysis_plots(processed_data)

library(forcats)
# 1. Prot performance by prior type
#load(paste0(simulation_output_folder, "/", "simulation_results_processed.RData"))
str(processed_data)
Prot_delta_data <- processed_data %>%
  filter(component == "Baseline") %>%
  filter(symm_method=="OR") %>%
  mutate(prior_type = fct_relevel(prior_type, "none", "noisy", "perfect"))

# Aggregate by configuration
delta_summary_data <- Prot_delta_data %>%
  group_by(n_samples, n_nodes, n_covariates, prior_type, symm_method) %>%
  summarise(
    mean_TPR = mean(TPR, na.rm = TRUE),
    mean_FPR = mean(FPR, na.rm = TRUE),
    mean_F1 = mean(F1, na.rm = TRUE),
    mean_Magnitude_preserved= mean(Magnitude_preserved, na.rm = TRUE),
    mean_Accuracy = mean(Accuracy, na.rm = TRUE),
    sd_TPR = sd(TPR, na.rm = TRUE),
    sd_FPR = sd(FPR, na.rm = TRUE),
    sd_F1 = sd(F1, na.rm = TRUE),
    sd_Magnitude_preserved = sd(Magnitude_preserved, na.rm = TRUE),
    sd_Accuracy = sd(Accuracy, na.rm = TRUE),
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

new_plots$delta_0_accuracy_comparison <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_Accuracy, 
                                                                        fill = prior_type, color= prior_type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = pmax(0, mean_Accuracy - sd_Accuracy), 
                    ymax = mean_Accuracy + sd_Accuracy,
                    color= prior_type),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(n_nodes ~ n_covariates, labeller = label_both, scales = "free_y") +
  labs(title = "Performance in network reconstruction: Accuracy by Prior Type",
       subtitle = "Perfromance metrics measured on the baseline network for which the prior is provided",
       x = "Sample Size", y = "Accuracy",
       fill = "Prior Type") +
  guides(colour = "none")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

new_plots$delta_0_accuracy_comparison
new_plots$delta_0_f1_comparison
new_plots$delta_0_error_comparison



delta_component_data <- processed_data %>%
  filter(component_type == "delta_individual") %>%
  filter(n_covariates==1) %>%
  filter(symm_method == "OR")
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

new_plots$delta_0_1_component_Acc

delta_component_data <- processed_data %>%
  filter(TP + FN == 0) %>%
  filter(component_type=="delta_individual") %>%
  filter(n_covariates==3)%>%
  filter(symm_method=="OR")
delta_component_data <- delta_component_data %>%
  mutate(prior_type = fct_relevel(prior_type, "none", "noisy", "perfect"))


delta_comp_summary <- delta_component_data %>%
  group_by(n_samples, prior_type, n_nodes) %>%
  summarise(
    mean_fpr = mean(FPR, na.rm = TRUE),
    mean_accuracy = mean(Accuracy, na.rm = TRUE),
    sd_fpr = sd(FPR, na.rm = TRUE),
    sd_accuracy = sd(Accuracy, na.rm = TRUE),
    .groups = 'drop'
  )

new_plots$delta_null_component_Acc <- ggplot(delta_comp_summary, aes(x = prior_type, y = mean_accuracy, 
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

new_plots$delta_null_component_fpr <- ggplot(delta_comp_summary, aes(x = prior_type, y = mean_fpr, 
                                                                     fill = prior_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_fpr - sd_fpr), 
                    ymax = mean_fpr + sd_fpr),
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

new_plots$delta_null_component_fpr_rest <- ggplot(delta_comp_summary_rest, aes(x = prior_type, y = mean_fpr, 
                                                                               fill = prior_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_fpr - sd_fpr), 
                    ymax = mean_fpr + sd_fpr),
                position = position_dodge(0.8), width = 0.2) +
  labs(title = "Performance in network reconstruction: FPR by Prior Type",
       subtitle = "Perfromance metrics measured on the network with zero prec matrix",
       x = "Component", y = "Mean FPR",
       fill = "Prior Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


new_plots$delta_null_component_fpr

#################################################
## PREPARE PLOTS FOR THE PAPER
#################################################
library(dplyr)
library(ggplot2)
library(forcats)
library(gridExtra)
library(grid)

# Define consistent color scheme
prior_cols <- c(
  "none"    = "#E74C3C",  # Red for no prior
  "noisy"   = "#27AE60",  # Green for noisy prior  
  "perfect" = "#3498DB"   # Blue for perfect prior
)


Prot_delta_data <- processed_data %>%
  filter(component == "Baseline") %>%
  filter(symm_method=="OR") %>%
  mutate(prior_type = fct_relevel(prior_type, "none", "noisy", "perfect"))


# Aggregate by configuration
delta_summary_data <- Prot_delta_data %>%
  group_by(n_samples, n_nodes, n_covariates, prior_type, symm_method) %>%
  summarise(
    mean_TPR = mean(TPR, na.rm = TRUE),
    mean_FPR = mean(FPR, na.rm = TRUE),
    mean_F1 = mean(F1, na.rm = TRUE),
    mean_Magnitude_preserved= mean(Magnitude_preserved, na.rm = TRUE),
    mean_Accuracy = mean(Accuracy, na.rm = TRUE),
    sd_TPR = sd(TPR, na.rm = TRUE),
    sd_FPR = sd(FPR, na.rm = TRUE),
    sd_F1 = sd(F1, na.rm = TRUE),
    sd_Magnitude_preserved = sd(Magnitude_preserved, na.rm = TRUE),
    sd_Accuracy = sd(Accuracy, na.rm = TRUE),
    .groups = 'drop'
  )

# FIGURE 1
p <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_Accuracy,
                                    fill = prior_type, color = prior_type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_Accuracy - sd_Accuracy),
        ymax = mean_Accuracy + sd_Accuracy),
    position = position_dodge(0.8), width = 0.2
  ) +
  facet_grid(rows = vars(n_nodes),
             cols = vars(n_covariates),
             labeller = label_both,
             switch = "y") +
  labs(
    title = "Performance in network reconstruction: Baseline network",
    x = "Sample Size", y = "Accuracy",
    fill = "Prior Type", color = "Prior Type"
  ) +
  scale_fill_manual(values = prior_cols) +
  scale_color_manual(values = prior_cols) +
  guides(colour = "none") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )+
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(face = "bold")
  )  +
  coord_cartesian(ylim = c(0.85, 1))
p
ggsave(p, filename="Simulation_results/Figure_1.png", width=7, height=5, dpi=300)

# FIGURE 2
p2 <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_F1, 
                                    fill = prior_type, color = prior_type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = pmax(0, mean_F1 - sd_F1), 
                    ymax = pmin(1, mean_F1 + sd_F1),
                    color = prior_type),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(rows = vars(n_nodes),
             cols = vars(n_covariates),
             labeller = label_both,
             switch = "y") +
  labs(title = "Performance in network reconstruction: Baseline network",
       x = "Sample Size", y = "F1 Score",
       fill = "Prior Type") +
  scale_fill_manual(values = prior_cols) +
  scale_color_manual(values = prior_cols) +
  guides(colour = "none")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(face = "bold")
  )
p2
ggsave(p2, filename="Simulation_results/Figure_2.png", width=7, height=5, dpi=300)

# FIGURE 3
p3 <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_Magnitude_preserved, 
                                     fill = prior_type, color= prior_type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = pmax(0, mean_Magnitude_preserved - sd_Magnitude_preserved), 
                    ymax = mean_Magnitude_preserved + sd_Magnitude_preserved,
                    color= prior_type),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(rows = vars(n_nodes),
             cols = vars(n_covariates),
             labeller = label_both,
             switch = "y") +
  labs(title = "Performance in network reconstruction: Baseline network",
       x = "Sample Size", y = "Spearman correlation",
       fill = "Prior Type") +
  scale_fill_manual(values = prior_cols) +
  scale_color_manual(values = prior_cols) +
  guides(colour = "none")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(face = "bold")
  )
p3
ggsave(p3, filename="Simulation_results/Figure_3.png", width=7, height=5, dpi=300)

# FIGURE 4
delta_component_data <- processed_data %>%
  filter(component_type == "delta_individual") %>%
  filter(n_covariates==1) %>%
  filter(symm_method == "OR") %>%
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

p4 <- ggplot(delta_comp_summary, aes(x = component, y = mean_accuracy, 
                                     fill = prior_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_accuracy - sd_accuracy), 
                    ymax = mean_accuracy + sd_accuracy),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(n_nodes ~ n_samples, labeller = label_both) +
  labs(title = "Performance in network reconstruction: Baseline and covariance network",
       subtitle = "Setting q = 1",
       x = "Component", y = "Mean accuracy",
       fill = "Prior Type") +
  scale_fill_manual(values = prior_cols) +
  scale_color_manual(values = prior_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.3),
        plot.subtitle = element_text(hjust = 0.5))+
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(face = "bold")
  )

p4
ggsave(p4, filename="Simulation_results/Figure_4.png", width=7, height=5, dpi=300)

# FIGURE 5
delta_component_data <- processed_data %>%
  filter(TP + FN == 0) %>%
  filter(component_type=="delta_individual") %>%
  filter(n_covariates==3)%>%
  filter(symm_method=="OR") %>%
  mutate(prior_type = fct_relevel(prior_type, "none", "noisy", "perfect"))


delta_comp_summary <- delta_component_data %>%
  group_by(n_samples, prior_type, n_nodes) %>%
  summarise(
    mean_fpr = mean(FPR, na.rm = TRUE),
    mean_accuracy = mean(Accuracy, na.rm = TRUE),
    sd_fpr = sd(FPR, na.rm = TRUE),
    sd_accuracy = sd(Accuracy, na.rm = TRUE),
    .groups = 'drop'
  )

p5 <- ggplot(delta_comp_summary, aes(x = prior_type, y = mean_fpr, 
                                     fill = prior_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pmax(0, mean_fpr - sd_fpr), 
                    ymax = mean_fpr + sd_fpr),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(n_nodes ~ n_samples, labeller = label_both) +
  labs(title = "Performance in network reconstruction: Covariance network",
       subtitle = "Setting q = 3",
       x = "Component", y = "Mean FPR",
       fill = "Prior Type") +
  scale_fill_manual(values = prior_cols) +
  scale_color_manual(values = prior_cols) +
  scale_x_discrete(labels = c("", "x_null", ""))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(face = "bold")
  )
p5
ggsave(p5, filename="Simulation_results/Figure_5.png", width=7, height=5, dpi=300)


#################################################
## FIGURE 1: BASELINE NETWORK RECONSTRUCTION PERFORMANCE
#################################################

library(patchwork)
library(ggplot2)

# 1) Make each plot not repeat legend/title (keep legend only once)
pA <- p  + labs(title = NULL) + theme(legend.position = "none")
pB <- p2 + labs(title = NULL) + theme(legend.position = "none")
pC <- p3 + labs(title = NULL) + theme(legend.position = "none")

pA <- pA + theme(axis.title.x = element_blank())
pB <- pB + theme(axis.title.x = element_blank())

# 2) Combine with shared legend on the right
fig1 <- (pA / pB / pC) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# 3) Add panel tags and a global title
fig1 <- fig1 +
  plot_annotation(
    title = "Baseline Network Reconstruction Performance",
    tag_levels = "A"
  ) &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# draw
fig1

# save (SVG for paper)
ggsave("Simulation_results/Figure_1_combined.svg", fig1, width = 8, height = 12, units = "in")
ggsave("Simulation_results/Figure_1_combined.png", fig1, width = 8, height = 12, units = "in", dpi = 300)


pA <- p4  + labs(title = NULL) + theme(legend.position = "none")
pB <- p5  + labs(title = NULL) + theme(legend.position = "none")
# 2) Combine with shared legend on the right
fig2 <- (pA / pB) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

fig2 <- fig2 +
  plot_annotation(
    title = "Covariance Network Reconstruction Performance",
    tag_levels = "A"
  ) &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# draw
fig2

# save (SVG for paper)
ggsave("Simulation_results/Figure_2_combined.svg", fig2, width = 8, height = 12, units = "in")
ggsave("Simulation_results/Figure_2_combined.png", fig2, width = 8, height = 12, units = "in", dpi = 300)

#################################################
## PREPARE TABLE FOR THE PAPER
#################################################
cat("Generating summary tables...\n")
tables <- generate_summary_tables(processed_data)
save(tables, file = paste0(simulation_output_folder, "/", "simulation_results_summary_table.RData"))


tab <- tables$delta_component_complete_summary %>%
  filter(symm_method == "OR")
write.csv(tab, file=paste0(simulation_output_folder, "/Table_simulation_summary.csv"))
