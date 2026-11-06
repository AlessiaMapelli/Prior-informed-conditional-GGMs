#################################################
## Run one simulation setting (sample size x replication)
#################################################
# Usage (from the Add_simulation_analysis folder):
#   Rscript run_setting.R --n 5000 --rep 1 [--p 100] [--prior noisy]
#                         [--ncores 16] [--nodes all|1-16] [--stab_B 50]
#                         [--out results]
# The node regressions run in parallel over --ncores cores (forked workers).
# Each node result is saved as soon as it is ready, so a job that hits the
# time limit can be resubmitted and continues where it stopped. Once all p
# nodes are available the setting is evaluated. With --nodes 1-16 only those
# nodes are fitted, which is used to measure the computing time.

opt <- list(n = NA, rep = NA, p = 100, prior = "noisy", ncores = NA,
            nodes = "all", stab_B = 50, out = "results", curve_nodes = "1-5")
args <- commandArgs(trailingOnly = TRUE)
for (k in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[k])]] <- args[k + 1]
for (v in c("n", "rep", "p", "stab_B", "ncores")) opt[[v]] <- as.integer(opt[[v]])
if (is.na(opt$n) || is.na(opt$rep)) stop("--n and --rep are required")
if (is.na(opt$ncores)) {
  opt$ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores()))
}
parse_range <- function(x, p) {
  if (x == "all") return(seq_len(p))
  r <- as.integer(strsplit(x, "-")[[1]])
  if (length(r) == 1) r else seq(r[1], min(r[2], p))
}

source("R/explore_functions.R")
load_original_functions("code")

setting_dir <- file.path(opt$out, sprintf("n%d_p%d_%s", opt$n, opt$p, opt$prior),
                         paste0("rep", opt$rep))
node_dir <- file.path(setting_dir, "nodes")
dir.create(node_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Setting n=%d p=%d prior=%s rep=%d, %d cores\n",
            opt$n, opt$p, opt$prior, opt$rep, opt$ncores))

# Data: the seed depends on the replication and the sample size
t0 <- Sys.time()
data_file <- file.path(setting_dir, "sim_dataset.rds")
if (file.exists(data_file)) {
  sim_dataset <- readRDS(data_file)
} else {
  sim_dataset <- simulate_setting(opt$n, opt$p, opt$prior, seed = 1000 * opt$rep + opt$n,
                                  out_dir = setting_dir)
  saveRDS(sim_dataset, data_file)
}
inp <- prepare_inputs(sim_dataset)
cat("Data ready in", format(Sys.time() - t0), "\n")

# Node regressions
nodes <- parse_range(opt$nodes, opt$p)
curve_nodes <- parse_range(opt$curve_nodes, opt$p)
todo <- nodes[!file.exists(file.path(node_dir, sprintf("node_%d.rds", nodes)))]
cat(length(todo), "of", length(nodes), "nodes to fit\n")
t0 <- Sys.time()
status <- parallel::mclapply(todo, function(i) {
  res <- tryCatch(explore_node(i, inp, stab_B = opt$stab_B,
                               keep_curves = i %in% curve_nodes),
                  error = function(e) e)
  if (inherits(res, "error")) return(conditionMessage(res))
  saveRDS(res, file.path(node_dir, sprintf("node_%d.rds", i)))
  "ok"
}, mc.cores = opt$ncores, mc.preschedule = FALSE)
failed <- todo[unlist(status) != "ok"]
if (length(failed) > 0) {
  cat("Failed nodes:", failed, "\n")
  print(unique(unlist(status)[unlist(status) != "ok"]))
}
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Fitted %d nodes in %.1f s (%.1f s per node per core)\n",
            length(todo), wall, wall * opt$ncores / max(length(todo), 1)))

# Timing summary of the available nodes
files <- list.files(node_dir, pattern = "^node_\\d+\\.rds$", full.names = TRUE)
res <- lapply(files, readRDS)
times <- do.call(rbind, lapply(res, function(r) data.frame(node = r$node, t(r$time))))
times$wall_total <- wall
saveRDS(times, file.path(setting_dir, "times.rds"))
cat("Mean time per node (s):\n"); print(colMeans(times[, c("path", "cv", "stab", "refit")]))

# Evaluation once all nodes are available
if (length(res) == opt$p) {
  cat("Evaluating...\n")
  evaluation <- evaluate_setting(res, inp, sim_dataset)
  evaluation$n <- opt$n; evaluation$p <- opt$p; evaluation$prior <- opt$prior
  evaluation$rep <- opt$rep
  saveRDS(evaluation, file.path(setting_dir, "evaluation.rds"))
  cv_diag <- do.call(rbind, lapply(res, `[[`, "cv_diag"))
  cv_diag$n <- opt$n; cv_diag$rep <- opt$rep
  saveRDS(cv_diag, file.path(setting_dir, "cv_diag.rds"))
  curves <- Filter(Negate(is.null), lapply(res, function(r) if (!is.null(r$curves)) list(node = r$node, curves = r$curves)))
  saveRDS(curves, file.path(setting_dir, "cv_curves.rds"))
  cat("Done:", setting_dir, "\n")
}
