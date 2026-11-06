# Regenerate the simulation results from the saved node files, without
# refitting: re-collects every replication with the current
# ggReg_main_functions.R (Delta = -coef / sigma^2, OR/AND choice on the
# regression coefficients) and rebuilds simulation_results.RData, the
# processed results, the summary table and the figures, by running the
# "PARALLEL EXECUTION - EVALUATION OF THE OUTCOME" section of
# Sim_pipeline_run.R and everything after it.
# Run from the Simulation_analysis folder:
#   Rscript --vanilla Sim_regenerate_results.R
rm(list = ls(all = TRUE))

# matlib is loaded by Sim_functions.R but not used, and needs X11
src <- readLines("Sim_functions.R")
eval(parse(text = src[src != "library(matlib)"]))
source("ggReg_main_functions.R")

# The node files were written by the earlier code, which stored
# No_sim_Delta_hat_row = -sigma * coef. Recompute it from the stored
# regression coefficients and residual variance as the current node code
# does (-coef / sigma), without modifying the node files.
fun_src <- readLines("ggReg_main_functions.R")
old_line <- "No_sim_Delta_hat[i, ] <- node_result$No_sim_Delta_hat_row"
stopifnot(sum(grepl(old_line, fun_src, fixed = TRUE)) == 1)
fun_src <- sub(old_line,
               "No_sim_Delta_hat[i, ] <- -1/node_result$sigma * node_result$Prec_reg_matrix_row",
               fun_src, fixed = TRUE)
eval(parse(text = fun_src))

N_REPLICATIONS <- 10
simulation_output_folder <- "Simulation_results"
simulation_output_file <- "simulation_results.RData"
load("Simulation_results/simulation_results.RData")  # SIMULATION_GRID
input_computation <- read.csv("Simulation_results/input_computation_file.csv")

pipe <- readLines("Sim_pipeline_run.R")
start <- grep("^all_results <- list\\(\\)", pipe)
pdf(NULL)
eval(parse(text = pipe[start:length(pipe)]))
cat("Regeneration completed\n")
