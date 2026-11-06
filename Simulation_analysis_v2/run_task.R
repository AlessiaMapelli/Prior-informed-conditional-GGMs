#################################################
## Run one task of the simulation study: one dataset (n, p, q, replication)
## fitted under the three priors
#################################################
# Usage (from the Simulation_analysis_v2 folder):
#   Rscript run_task.R --n 5000 --p 100 --q 3 --rep 1 [--ncores 16]
#                      [--priors perfect,noisy,none] [--nodes all|1-16]
#                      [--selection ebic_relaxed] [--out results] [--tmp <dir>]
# The node-wise regressions use the template functions
# (ggReg_main_functions.R): GGReg_cov_single_node_processing() for each node,
# run in parallel on --ncores cores, and collect_node_results() with OR and
# AND symmetrisation. The input and node files are written to --tmp (default:
# the R temporary directory, i.e. the local disk of the compute node) and
# deleted at the end. With --nodes 1-16 only those nodes are fitted (timing
# test) and no evaluation is done.
# Outputs in <out>/n<n>_p<p>_q<q>/rep<rep>/:
#   evaluation.rds  metrics per prior, symmetrisation and component
#   nodes.rds       per prior and node: computing time, selected lambda,
#                   asparse, weight, number of selected predictors
#   config.rds      dataset configuration (seed, diagonal, edge counts) and
#                   data generation time

opt <- list(n = NA, p = NA, q = NA, rep = NA, ncores = NA, priors = "perfect,noisy,none",
            nodes = "all", selection = "ebic_relaxed", out = "results", tmp = tempdir())
args <- commandArgs(trailingOnly = TRUE)
for (k in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[k])]] <- args[k + 1]
for (v in c("n", "p", "q", "rep", "ncores")) opt[[v]] <- as.integer(opt[[v]])
if (any(is.na(unlist(opt[c("n", "p", "q", "rep")])))) stop("--n, --p, --q and --rep are required")
if (is.na(opt$ncores)) {
  opt$ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores()))
}
priors <- strsplit(opt$priors, ",")[[1]]
nodes <- if (opt$nodes == "all") seq_len(opt$p) else {
  r <- as.integer(strsplit(opt$nodes, "-")[[1]]); seq(r[1], min(r[2], opt$p))
}
full_run <- length(nodes) == opt$p

suppressPackageStartupMessages({
  library(sparsegl)
  library(glmnet)
})
source("ggReg_main_functions.R")
source("R/sim_design.R")

task <- sprintf("n%d_p%d_q%d", opt$n, opt$p, opt$q)
out_dir <- file.path(opt$out, task, paste0("rep", opt$rep))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (full_run && file.exists(file.path(out_dir, "evaluation.rds"))) {
  cat("Already evaluated:", out_dir, "\n")
  quit(save = "no")
}
cat(sprintf("Task %s rep %d: priors %s, %d nodes, %d cores (%s)\n", task, opt$rep,
            opt$priors, length(nodes), opt$ncores, Sys.info()[["nodename"]]))

# Data
seed <- dataset_seed(opt$n, opt$p, opt$q, opt$rep)
t0 <- Sys.time()
sim <- simulate_dataset(opt$n, opt$p, opt$q, seed, ncores = opt$ncores)
time_data <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Data generated in %.1f s (diagonal %.2f, baseline target edges %d)\n",
            time_data, sim$config$diagonal, sim$config$target_edges))

# Fixed settings of the template (see README)
settings <- list(scr = FALSE, gamma = NULL, scr_matrix = NULL, lambda_prec = NULL,
                 selection = opt$selection, refit = TRUE, selection_control = list(),
                 tune_hyperparams = TRUE, asparse_grid = c(0.5, 0.75, 0.9, 0.99),
                 weight_grid = c(0.8, 1, 1.1, 1.3), random_hyper_search = TRUE,
                 p.rand.hyper = 0.5, K = 5)

evaluation <- list(); node_info <- list(); wall <- c()
for (pr in priors) {
  work <- file.path(opt$tmp, sprintf("%s_rep%d_%s_%d", task, opt$rep, pr, Sys.getpid()))
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  work <- paste0(normalizePath(work), "/")
  input_data_path <- paste0(work, "input_data_nodes.rda")
  local({
    Z0 <- sim$x
    known_ppi <- make_prior(sim$T0, pr, seed + 7)
    covariates <- sim$covariates
    list2env(settings, environment())
    save(Z0, known_ppi, covariates, scr, gamma, scr_matrix, lambda_prec, selection, refit,
         selection_control, tune_hyperparams, asparse_grid, weight_grid,
         random_hyper_search, p.rand.hyper, K, file = input_data_path)
  })
  t0 <- Sys.time()
  status <- parallel::mclapply(nodes, function(i) {
    r <- tryCatch({
      invisible(capture.output(GGReg_cov_single_node_processing(
        input_data_path = input_data_path, node_index = i,
        output_path = work, name_output = "node")))
      "ok"
    }, error = function(e) conditionMessage(e))
    r
  }, mc.cores = opt$ncores, mc.preschedule = FALSE)
  wall[pr] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  bad <- unlist(status) != "ok"
  if (any(bad)) {
    cat("Failed nodes (", pr, "):", nodes[bad], "\n"); print(unique(unlist(status)[bad]))
    stop("Node failures")
  }
  # Per-node information
  node_info[[pr]] <- do.call(rbind, lapply(nodes, function(i) {
    load(paste0(work, "node_node_", i, ".rda"))
    op <- node_result$optimal_params
    data.frame(prior = pr, node = i,
               time_s = as.numeric(node_result$computation_time, units = "secs"),
               lambda = op$best_lambda_prec, asparse = op$asparse, weight = op$weight,
               n_selected = op$n_selected)
  }))
  cat(sprintf("Prior %-7s: %d nodes in %.1f s (%.1f s per node, mean node time %.1f s)\n",
              pr, length(nodes), wall[pr], wall[pr] * opt$ncores / length(nodes),
              mean(node_info[[pr]]$time_s)))
  # Collection and evaluation (full runs only)
  if (full_run) {
    for (symm in c("OR", "AND")) {
      res <- collect_node_results(opt$p, work, "node", symm)
      ev <- evaluate_networks(res$Dic_Delta_hat, sim)
      evaluation[[length(evaluation) + 1]] <- cbind(prior = pr, symm = symm, ev)
    }
  }
  unlink(work, recursive = TRUE)
}

node_info <- do.call(rbind, node_info)
node_info$n <- opt$n; node_info$p <- opt$p; node_info$q <- opt$q; node_info$rep <- opt$rep
saveRDS(node_info, file.path(out_dir, if (full_run) "nodes.rds" else "nodes_timing.rds"))
saveRDS(c(sim$config, list(time_data_s = time_data, wall_s = wall, ncores = opt$ncores,
                           selection = opt$selection)),
        file.path(out_dir, if (full_run) "config.rds" else "config_timing.rds"))
if (full_run) {
  evaluation <- do.call(rbind, evaluation)
  evaluation <- cbind(n = opt$n, p = opt$p, q = opt$q, rep = opt$rep, evaluation)
  saveRDS(evaluation, file.path(out_dir, "evaluation.rds"))
}
cat("Done:", out_dir, "\n")
