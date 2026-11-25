#################################################
## Reproducibility experiment: one task = one dataset drawn from a FIXED
## true network, fitted with several selection rules and priors.
#################################################
# Design (README.md, "Reproducibility experiment"):
# - The true networks (baseline and covariate networks) and the prior are fixed
#   for a (p, q) setting: graph seed repro_graph_seed(p, q), the same for all n
#   and all datasets; prior seed = graph seed + 7.
# - Datasets r = 1, ..., R differ only in the covariates and the subjects
#   (seed repro_data_seed(n, p, q, r)), so the variability of the estimated
#   networks across r is sampling variability only.
# - Each dataset is fitted with each rule in --rules and each prior in --priors,
#   with the main-study settings otherwise (grid, random hyperparameter search).
# Usage (from the Simulation_analysis_v2 folder):
#   Rscript run_repro_task.R --n 5000 --rep 1 [--p 100 --q 3] [--ncores 16]
#     [--rules bic_relaxed,ebic_relaxed,stability_bic] [--priors noisy,none]
#     [--nodes all|1-16] [--out results_reproducibility] [--tmp <dir>]
# Outputs in <out>/n<n>_p<p>_q<q>/rep<rep>/:
#   networks.rds    estimated edges (upper triangle, non-zero entries of the
#                   symmetrised matrices): rule, prior, symm, component, i, j, value
#   evaluation.rds  evaluate_networks() per rule, prior, symmetrisation, component
#   truth.rds       true edges on the scale of the fit (component, i, j, value)
#   nodes.rds       per rule, prior and node: time, lambda, asparse, weight, n_selected
#   config.rds      dataset configuration and wall times

opt <- list(n = NA, p = 100, q = 3, rep = NA, ncores = NA,
            rules = "bic_relaxed,ebic_relaxed,stability_bic", priors = "noisy,none",
            nodes = "all", out = "results_reproducibility", tmp = tempdir())
args <- commandArgs(trailingOnly = TRUE)
for (k in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[k])]] <- args[k + 1]
for (v in c("n", "p", "q", "rep", "ncores")) opt[[v]] <- as.integer(opt[[v]])
if (any(is.na(unlist(opt[c("n", "rep")])))) stop("--n and --rep are required")
if (is.na(opt$ncores)) {
  opt$ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores()))
}
rules <- strsplit(opt$rules, ",")[[1]]
priors <- strsplit(opt$priors, ",")[[1]]
nodes <- if (opt$nodes == "all") seq_len(opt$p) else {
  r <- as.integer(strsplit(opt$nodes, "-")[[1]]); seq(r[1], min(r[2], opt$p))
}
full_run <- length(nodes) == opt$p

suppressPackageStartupMessages({
  library(sparsegl)
  library(glmnet)
})
source("../Computational_templates/ggReg_main_functions.R")  # single template copy
source("R/sim_design.R")

#' Seeds of the reproducibility experiment (disjoint from dataset_seed(),
#' which is below 1.1e7)
repro_graph_seed <- function(p, q) as.integer(5e7 + p * 10 + q)
repro_data_seed <- function(n, p, q, rep) {
  as.integer(6e7 + match(n, c(1000, 5000, 10000, 20000)) * 1e6 + p * 1000 + q * 100 + rep)
}

task <- sprintf("n%d_p%d_q%d", opt$n, opt$p, opt$q)
out_dir <- file.path(opt$out, task, paste0("rep", opt$rep))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (full_run && file.exists(file.path(out_dir, "networks.rds"))) {
  cat("Already done:", out_dir, "\n")
  quit(save = "no")
}
cat(sprintf("Reproducibility task %s rep %d: rules %s, priors %s, %d nodes, %d cores (%s)\n",
            task, opt$rep, opt$rules, opt$priors, length(nodes), opt$ncores, Sys.info()[["nodename"]]))

gseed <- repro_graph_seed(opt$p, opt$q)
seed <- repro_data_seed(opt$n, opt$p, opt$q, opt$rep)
t0 <- Sys.time()
sim <- simulate_dataset(opt$n, opt$p, opt$q, seed, ncores = opt$ncores, graph_seed = gseed)
time_data <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Data generated in %.1f s (diagonal %.2f, baseline target edges %d)\n",
            time_data, sim$config$diagonal, sim$config$target_edges))

edges_of <- function(M) {
  M <- as.matrix(M); ut <- which(upper.tri(M) & M != 0, arr.ind = TRUE)
  data.frame(i = ut[, 1], j = ut[, 2], value = M[ut])
}
truth <- do.call(rbind, lapply(names(sim$truth), function(cn) {
  e <- edges_of(ifelse(abs(sim$truth[[cn]]) > 1e-8, sim$truth[[cn]], 0))
  if (nrow(e)) cbind(component = cn, e) else NULL
}))

settings <- list(scr = FALSE, gamma = NULL, scr_matrix = NULL, lambda_prec = NULL,
                 refit = TRUE, selection_control = list(),
                 tune_hyperparams = TRUE, asparse_grid = c(0.5, 0.75, 0.9, 0.99),
                 weight_grid = c(0.8, 1, 1.1, 1.3), random_hyper_search = TRUE,
                 p.rand.hyper = 0.5, K = 5)

networks <- list(); evaluation <- list(); node_info <- list(); wall <- c()
for (pr in priors) {
  known_ppi_pr <- make_prior(sim$T0, pr, gseed + 7)
  for (rule in rules) {
    work <- file.path(opt$tmp, sprintf("repro_%s_rep%d_%s_%s_%d", task, opt$rep, pr, rule, Sys.getpid()))
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    work <- paste0(normalizePath(work), "/")
    input_data_path <- paste0(work, "input_data_nodes.rda")
    local({
      Z0 <- sim$x; known_ppi <- known_ppi_pr; covariates <- sim$covariates
      selection <- rule
      list2env(settings, environment())
      save(Z0, known_ppi, covariates, scr, gamma, scr_matrix, lambda_prec, selection, refit,
           selection_control, tune_hyperparams, asparse_grid, weight_grid,
           random_hyper_search, p.rand.hyper, K, file = input_data_path)
    })
    t0 <- Sys.time()
    status <- parallel::mclapply(nodes, function(i) {
      tryCatch({
        invisible(capture.output(GGReg_cov_single_node_processing(
          input_data_path = input_data_path, node_index = i,
          output_path = work, name_output = "node")))
        "ok"
      }, error = function(e) conditionMessage(e))
    }, mc.cores = opt$ncores, mc.preschedule = FALSE)
    key <- paste(rule, pr, sep = "|")
    wall[key] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    bad <- unlist(status) != "ok"
    if (any(bad)) {
      cat("Failed nodes (", key, "):", nodes[bad], "\n"); print(unique(unlist(status)[bad]))
      stop("Node failures")
    }
    node_info[[key]] <- do.call(rbind, lapply(nodes, function(i) {
      load(paste0(work, "node_node_", i, ".rda"))
      op <- node_result$optimal_params
      data.frame(rule = rule, prior = pr, node = i,
                 time_s = as.numeric(node_result$computation_time, units = "secs"),
                 lambda = op$best_lambda_prec, asparse = op$asparse, weight = op$weight,
                 n_selected = op$n_selected)
    }))
    cat(sprintf("%-14s %-6s: %d nodes in %.1f s (mean node time %.1f s)\n", rule, pr,
                length(nodes), wall[key], mean(node_info[[key]]$time_s)))
    if (full_run) {
      for (symm in c("OR", "AND")) {
        res <- collect_node_results(opt$p, work, "node", symm)
        for (cn in names(sim$truth)) {
          e <- edges_of(res$Dic_Delta_hat[[cn]])
          if (nrow(e)) networks[[length(networks) + 1]] <-
            cbind(rule = rule, prior = pr, symm = symm, component = cn, e)
        }
        ev <- evaluate_networks(res$Dic_Delta_hat, sim)
        evaluation[[length(evaluation) + 1]] <- cbind(rule = rule, prior = pr, symm = symm, ev)
      }
    }
    unlink(work, recursive = TRUE)
  }
}

meta <- data.frame(n = opt$n, p = opt$p, q = opt$q, rep = opt$rep)
node_info <- cbind(meta, do.call(rbind, node_info), row.names = NULL)
saveRDS(node_info, file.path(out_dir, if (full_run) "nodes.rds" else "nodes_timing.rds"))
saveRDS(c(sim$config, list(time_data_s = time_data, wall_s = wall, ncores = opt$ncores)),
        file.path(out_dir, if (full_run) "config.rds" else "config_timing.rds"))
if (full_run) {
  saveRDS(cbind(meta, truth, row.names = NULL), file.path(out_dir, "truth.rds"))
  saveRDS(cbind(meta, do.call(rbind, evaluation), row.names = NULL), file.path(out_dir, "evaluation.rds"))
  # written last: its presence marks a completed task
  saveRDS(cbind(meta, do.call(rbind, networks), row.names = NULL), file.path(out_dir, "networks.rds"))
}
cat("Done:", out_dir, "\n")
