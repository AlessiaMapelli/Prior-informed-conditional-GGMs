# Simulation Pipeline for Test Data in Prior-informed Covariate-Dependent Gaussian Graphical Models
setwd("/group/diangelantonio/users/alessia_mapelli/Complete_projects/Prior-informed-conditional-GGMs/")
rm(list = ls(all = TRUE))
library(fs)

# Source required functions
source("Simulation_analysis/Sim_functions.R")

#################################################
## SIMULATION CONFIGURATION
#################################################

# Define simulation grid
SIMULATION_GRID <- expand.grid(
  n_samples = c(1000),
  n_nodes = c(10),
  prior_type = c("perfect"),
  n_covariates = c(1),  # 1=binary only, 3=2 continuous + 1 binary, 5=3 continuous (1 no effect) + 2 binary(1 no effect)
  symm_method = c("OR"), 
  stringsAsFactors = FALSE
)


# Number of replications per configuration
N_REPLICATIONS <- 1

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


#################################################
## SIMULATION EXECUTION FUNCTIONS
#################################################
simulation_output_folder = "Computational_templates"
dir.create(simulation_output_folder)
simulation_output_file   = "Example_data.RData"


#################################################
## PARALLEL EXECUTION - GENERATE INPUT DATASETS FOR ALL CONFIGURATIONS
#################################################

input_computation <- generate_input_datasets_simulation(
  SIMULATION_GRID = SIMULATION_GRID,
  N_REPLICATIONS = N_REPLICATIONS,
  output_folder= simulation_output_folder,
  output_file = simulation_output_file
)

load("Computational_templates/n1000_p10_q1_perfect_OR/rep1/input_data_nodes.rda")
expression_data <- Z0
covariate_data <- covariates
covariate_data$x1 <- as.factor(covariate_data$x1)
ppi_weight_matrix <- known_ppi
save(expression_data, covariate_data,ppi_weight_matrix, file= "Computational_templates/Example_data.RData")

write.csv(expression_data,file= "Computational_templates/expression_data.csv")
write.csv(covariate_data,file= "Computational_templates/covariate_data.csv")
write.table(ppi_weight_matrix, file="Computational_templates/ppi_weight_matrix.txt")

file_delete("Computational_templates/n1000_p10_q1_perfect_OR")
file_delete("Computational_templates/input_computation_file.csv")

