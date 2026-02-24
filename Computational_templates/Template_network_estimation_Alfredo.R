rm(list=ls(all=TRUE))

# Set a base dir the one where you clone the git directory
base_dir = "/group/diangelantonio/users/alessia_mapelli/Complete_projects/Prior-informed-conditional-GGMs/Computational_templates/"
setwd(base_dir)

source(paste0(base_dir, "ggReg_main_functions.R"))

# Set you events embeddings dataframe path
data_path <- paste0(base_dir, "phrase_embeddings_df.csv")

# Check if file exists
if (!file.exists(data_path)) {
  stop("Data file not found. Please check the path: ", data_path)
}

# Load the data
cat("Loading UK Biobank disease embedding data...\n")
ukb_data <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE)


results <- GGReg_full_estimation(
  x = ukb_data,
  known_ppi = NULL,
  covariates = NULL,
  scr = FALSE,
  mean_estimation = FALSE,
  lambda_prec_type = "1se",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = TRUE,
  p.rand.hyper = 0.5,
  K = 5,
  use_slurm = TRUE,
  slurm_script_path = "./slurm_ggReg_node.sbatch",
  output_path = "./results/",
  name_output = "ggReg_result",
  symm_method ="AND",
  verbose = FALSE)
weighted_adj_matrix <- results$Dic_adj_matrics$Baseline
plot_personalized_network(weighted_adj_matrix)
save(weighted_adj_matrix, file = "adj_matrix_example.RData")
