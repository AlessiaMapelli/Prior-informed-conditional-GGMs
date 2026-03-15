#!/usr/bin/env Rscript

# ggReg_cov_node_job.R
# Individual node estimation script for SLURM array execution
# Processes a single node i with optional hyperparameter tuning

source("ggReg_main_functions.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript ggReg_node_job.R input_data_path node_index output_path name_output")
}

GGReg_cov_single_node_processing(
      input_data_path = args[1],
      node_index = as.numeric(args[2]),
      output_path = args[3],
      name_output = args[4],
      verbose = TRUE
)
