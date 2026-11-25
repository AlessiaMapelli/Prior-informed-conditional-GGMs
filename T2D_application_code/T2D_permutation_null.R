#!/usr/bin/env Rscript
###############################################################################
# T2D case study - permutation null for the T2D-specific network component
#
# Question: are the T2D-specific edges (component Diab_Incident) more than
# expected by chance? (5-year endpoint: ~270 incident cases in the selection
# set; 10-year endpoint: 1,112.)
#
# Null: incident-T2D labels are permuted within sex x age-quintile strata of
# the selection set (preserving the T2D prevalence by sex and age); Sex and Age
# are untouched. For each permutation:
#   1. the mean step is re-run with the permuted design (same CV folds);
#   2. every node is refitted at the hyperparameters (lambda, alpha, theta)
#      it selected on the real data: the sparse-group lasso path is run from
#      lambda_max down to the selected lambda with the same log-spacing as the
#      original path (warm starts), followed by the OLS refit of the support.
#      The grid search is not repeated (~1/16 of the cost or less);
#   3. node-wise coefficients are symmetrised with AND and OR.
# Permutation 0 uses the real labels with the same procedure: it checks that
# the fixed-hyperparameter refit reproduces the original fit.
#
# Usage (from the private project folder, R with prCGGM >= 0.2.0):
#   Rscript --vanilla ../T2D_application_code/T2D_permutation_null.R prepare [B=100]
#   Rscript --vanilla T2D_permutation_null.R node <task_id> [nchunks]
#       (called by slurm_permutation_null.sbatch; task -> permutation x chunk)
#   Rscript --vanilla ../T2D_application_code/T2D_permutation_null.R collect
# The fit to test is chosen with T2D_RES (default results_bic_relaxed_10y),
# e.g. T2D_RES=results_bic_relaxed Rscript --vanilla ../T2D_application_code/T2D_permutation_null.R collect
###############################################################################

source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "config.R"))
use_prcggm()   # uses the internal prCGGM:::refit_ols() of version >= 0.2.0

# Results folder of the fit to test (environment variable T2D_RES):
#   results_bic_relaxed_10y (default; 10-year endpoint, primary)
#   results_bic_relaxed     (5-year endpoint, sensitivity analysis)
RES_DIR   <- Sys.getenv("T2D_RES", "results_bic_relaxed_10y")
BASE_DIR  <- file.path(RES_DIR, "prcggm_nodes")         # original fit
PERM_DIR  <- file.path(RES_DIR, "permutation_null")
NCHUNKS   <- 1L        # chunks of nodes per permutation (1 task ~ 3 h per permutation for 5y)
T2D_COL   <- "Diab_Incident"

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "help"
cat("Results folder:", RES_DIR, "\n")

base_input <- function() readRDS(file.path(BASE_DIR, "prcggm_input.rds"))

# ---------------------------------------------------------------------------
# Selected hyperparameters of the original fit, one row per node
# ---------------------------------------------------------------------------
read_original <- function(p) {
  files <- file.path(BASE_DIR, sprintf("prcggm_node_%d.rds", seq_len(p)))
  miss <- which(!file.exists(files))
  if (length(miss)) stop("Original fit incomplete; missing nodes: ",
                         paste(miss, collapse = ", "))
  res <- lapply(files, readRDS)
  list(tuning = do.call(rbind, lapply(res, `[[`, "tuning")),
       coef = do.call(rbind, lapply(res, `[[`, "coef_row")))
}

# ---------------------------------------------------------------------------
# Refit of node i at fixed (lambda, alpha, theta) - mirrors prCGGM:::fit_node
# ---------------------------------------------------------------------------
fit_node_fixed <- function(i, prep, lambda_sel, alpha, theta) {
  n <- prep$n; p <- prep$p; q <- prep$q
  y <- prep$z[, i]
  others <- seq_len(p)[-i]
  if (!is.null(prep$keep)) others <- others[prep$keep[i, others]]
  m <- length(others)
  zo <- prep$z[, others, drop = FALSE]
  X <- cbind(zo, do.call(cbind, lapply(seq_len(q), function(h) prep$U[, h] * zo)))
  pf_baseline <- if (is.null(prep$penalty)) rep(1, m) else prep$penalty[i, others]
  a <- list(x = X, y = y, group = rep(seq_len(q + 1), each = m),
            pf_group = c(0, rep(1, q)),
            pf_sparse = c(pf_baseline, rep(theta, q * m)),
            asparse = alpha, intercept = FALSE)
  # lambda_max of this design, then the original log-spacing
  # (default sparsegl path: nlambda = 100, lambda.factor = 1e-4 for n > P)
  # (sparsegl does not return lambda_max itself: from a 3-point log-spaced
  # path, lambda_max = lambda_2^2 / lambda_3)
  l3 <- do.call(sparsegl::sparsegl, c(a, list(nlambda = 3,
                                              lambda.factor = 0.5)))$lambda
  lmax <- l3[2]^2 / l3[3]
  coef_pen <- numeric(ncol(X))
  if (lambda_sel < lmax) {
    step <- log(1e4) / 99
    nl <- max(2L, ceiling(log(lmax / lambda_sel) / step) + 1L)
    lam <- exp(seq(log(lmax), log(lambda_sel), length.out = nl))
    path <- do.call(sparsegl::sparsegl, c(a, list(lambda = lam)))
    coef_pen <- as.matrix(path$beta)[, ncol(path$beta)]
  }
  support <- coef_pen != 0
  coef <- prCGGM:::refit_ols(X, y, support)
  coef_row <- numeric(1 + (q + 1) * p)
  coef_row[unlist(lapply(0:q, function(h) 1L + h * p + others))] <- coef
  coef_row
}

# Symmetrised edge sets (upper triangle) for component h
edges <- function(coef, p, h, rule) {
  D <- coef[, 1 + (h - 1) * p + seq_len(p)] != 0
  E <- if (rule == "AND") D & t(D) else D | t(D)
  E[upper.tri(E)]
}

# ---------------------------------------------------------------------------
if (mode == "prepare") {
  # Cheap: stores the original fit and the permuted labels; the mean step and
  # the node fits run in the array tasks (one task per permutation).
  B <- if (length(args) > 1) as.integer(args[2]) else 100L
  dir.create(PERM_DIR, recursive = TRUE, showWarnings = FALSE)
  prep <- base_input()
  stopifnot(T2D_COL %in% colnames(prep$U))
  orig <- read_original(prep$p)
  saveRDS(orig, file.path(PERM_DIR, "original.rds"))
  sex <- prep$U[, "Sex_Male"]
  age_q <- cut(prep$U[, "Age_c"], quantile(prep$U[, "Age_c"], 0:5 / 5),
               include.lowest = TRUE, labels = FALSE)
  strata <- interaction(sex, age_q, drop = TRUE)
  y0 <- prep$U[, T2D_COL]
  labels <- matrix(y0, length(y0), B + 1L)      # column 1 = real labels
  for (b in seq_len(B)) {
    set.seed(1000 + b)
    for (s in levels(strata)) {
      idx <- which(strata == s)
      labels[idx, b + 1L] <- y0[idx[sample.int(length(idx))]]
    }
  }
  overlap <- colSums(labels[, -1, drop = FALSE] == 1 & y0 == 1)
  saveRDS(list(labels = labels, strata = strata), file.path(PERM_DIR, "labels.rds"))
  cat("Cases:", sum(y0), "| strata:", nlevels(strata),
      "| cases kept in place by chance: mean", round(mean(overlap), 1),
      "max", max(overlap), "\n")
  cat(sprintf("Submit: sbatch --array=1-%d%%50 slurm_permutation_null.sbatch\n",
              (B + 1L) * NCHUNKS))

} else if (mode == "node") {
  task <- as.integer(args[2])
  nchunks <- if (length(args) > 2) as.integer(args[3]) else NCHUNKS
  b <- (task - 1L) %/% nchunks
  chunk <- (task - 1L) %% nchunks + 1L
  out <- file.path(PERM_DIR, sprintf("perm_%03d_chunk_%02d.rds", b, chunk))
  if (file.exists(out)) { cat("exists:", out, "\n"); quit(save = "no") }
  t0 <- Sys.time()
  prep <- base_input()
  orig <- readRDS(file.path(PERM_DIR, "original.rds"))
  lab <- readRDS(file.path(PERM_DIR, "labels.rds"))$labels
  # mean step with the permuted design (raw data reconstructed exactly from
  # residuals + fitted mean; same CV folds as the original fit)
  x <- prep$z + cbind(1, prep$U) %*% t(prep$mean_effects)
  prep$U[, T2D_COL] <- lab[, b + 1L]
  prep$z <- prCGGM:::fit_mean(x, prep$U, NULL, prep$settings$lambda_mean_type,
                              prep$foldid, FALSE)$z
  rm(x)
  cat(sprintf("perm %d: mean step done in %.1f min\n", b,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  nodes <- which(((seq_len(prep$p) - 1L) * nchunks) %/% prep$p + 1L == chunk)
  rows <- matrix(0, length(nodes), 1 + (prep$q + 1) * prep$p)
  for (k in seq_along(nodes)) {
    i <- nodes[k]; tu <- orig$tuning[i, ]
    rows[k, ] <- fit_node_fixed(i, prep, tu$lambda, tu$alpha, tu$theta)
  }
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("perm %d chunk %d: %d nodes in %.1f min\n", b, chunk,
              length(nodes), mins))
  saveRDS(list(b = b, nodes = nodes, coef = Matrix::Matrix(rows, sparse = TRUE),
               minutes = mins), out)

} else if (mode == "collect") {
  prep <- base_input()
  p <- prep$p; comps <- c("Baseline", colnames(prep$U))
  h_t2d <- match(T2D_COL, comps)
  orig <- readRDS(file.path(PERM_DIR, "original.rds"))
  files <- list.files(PERM_DIR, "^perm_\\d{3}_chunk_\\d{2}\\.rds$", full.names = TRUE)
  chunks <- lapply(files, readRDS)
  perms <- sort(unique(vapply(chunks, `[[`, integer(1), "b")))
  coef_by_perm <- list()
  for (b in perms) {
    cb <- chunks[vapply(chunks, `[[`, integer(1), "b") == b]
    C <- matrix(NA_real_, p, ncol(orig$coef))
    for (ch in cb) C[ch$nodes, ] <- as.matrix(ch$coef)
    if (anyNA(C)) { cat("permutation", b, "incomplete, skipped\n"); next }
    coef_by_perm[[as.character(b)]] <- C
  }
  stopifnot("0" %in% names(coef_by_perm))
  # check: permutation 0 (real labels) vs original fit
  same <- sapply(seq_along(comps), function(h)
    mean(edges(coef_by_perm[["0"]], p, h, "OR") == edges(orig$coef, p, h, "OR")))
  cat("Agreement of perm 0 with the original fit (OR edge sets):\n")
  print(setNames(round(same, 5), comps))
  # counts per permutation and component
  counts <- do.call(rbind, lapply(names(coef_by_perm), function(b) {
    do.call(rbind, lapply(seq_along(comps), function(h) data.frame(
      perm = as.integer(b), component = comps[h],
      AND = sum(edges(coef_by_perm[[b]], p, h, "AND")),
      OR = sum(edges(coef_by_perm[[b]], p, h, "OR")))))
  }))
  write.csv(counts, file.path(PERM_DIR, "edge_counts_by_permutation.csv"),
            row.names = FALSE)
  obs <- counts[counts$perm == 0, ]
  null <- counts[counts$perm > 0, ]
  B <- length(unique(null$perm))
  summ <- do.call(rbind, lapply(comps, function(cmp) {
    nl <- null[null$component == cmp, ]; ob <- obs[obs$component == cmp, ]
    data.frame(component = cmp, B = B,
               obs_AND = ob$AND, null_AND_mean = mean(nl$AND),
               null_AND_max = max(nl$AND),
               p_AND = (1 + sum(nl$AND >= ob$AND)) / (B + 1),
               obs_OR = ob$OR, null_OR_mean = mean(nl$OR),
               null_OR_max = max(nl$OR),
               p_OR = (1 + sum(nl$OR >= ob$OR)) / (B + 1))
  }))
  print(summ, row.names = FALSE)
  write.csv(summ, file.path(PERM_DIR, "permutation_summary.csv"), row.names = FALSE)
  # edge-wise: how often each observed T2D edge appears under the null
  ut <- which(upper.tri(diag(p)), arr.ind = TRUE)
  ew <- do.call(rbind, lapply(c("AND", "OR"), function(rule) {
    e0 <- edges(coef_by_perm[["0"]], p, h_t2d, rule)
    freq <- Reduce(`+`, lapply(names(coef_by_perm)[names(coef_by_perm) != "0"],
                               function(b) edges(coef_by_perm[[b]], p, h_t2d, rule)))
    k <- which(e0)
    if (!length(k)) return(NULL)
    data.frame(rule = rule, protein_1 = prep$nodes[ut[k, 1]],
               protein_2 = prep$nodes[ut[k, 2]], null_frequency = freq[k] / B)
  }))
  if (!is.null(ew)) write.csv(ew, file.path(PERM_DIR, "T2D_edges_null_frequency.csv"),
                              row.names = FALSE)
  # null frequency of every pair (does chance concentrate on some proteins?)
  hub <- Reduce(`+`, lapply(names(coef_by_perm)[names(coef_by_perm) != "0"], function(b) {
    D <- coef_by_perm[[b]][, 1 + (h_t2d - 1) * p + seq_len(p)] != 0
    rowSums(D | t(D))
  })) / B
  write.csv(data.frame(protein = prep$nodes, mean_null_T2D_degree_OR = hub),
            file.path(PERM_DIR, "T2D_null_degree_by_protein.csv"), row.names = FALSE)

} else {
  cat("Usage: T2D_permutation_null.R prepare [B] | node <task> [nchunks] | collect\n")
}
