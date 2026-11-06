#################################################
## Simulation study v2: data generation and evaluation
#################################################
# Design and rationale: README.md of this folder. The model fit uses the
# template functions (ggReg_main_functions.R, copy of
# Computational_templates/ggReg_main_functions.R).

suppressPackageStartupMessages({
  library(igraph)
})

#' Seed of a dataset: function of (n, p, q, replication), not of the prior,
#' so that the three priors are fitted on the same data
dataset_seed <- function(n, p, q, rep) {
  as.integer(rep * 100000 + match(n, c(1000, 5000, 10000, 20000)) * 10000 +
               match(p, c(50, 100, 250, 500)) * 1000 + q)
}

#' Covariates of the design
#' q = 1: x1 binary (Bernoulli(0.5)), with an effect.
#' q = 3: x1 continuous (U[0,1]) with an effect, x2 continuous (U[0,1])
#'        without effect, x3 binary (Bernoulli(0.5)) with an effect.
covariate_design <- function(q) {
  if (q == 1) {
    data.frame(name = "x1", type = "binary", effect = TRUE)
  } else if (q == 3) {
    data.frame(name = c("x1", "x2", "x3"), type = c("continuous", "continuous", "binary"),
               effect = c(TRUE, FALSE, TRUE))
  } else {
    stop("q must be 1 or 3")
  }
}

#' Simulate one dataset
#'
#' - Baseline network: power-law degree distribution with exponent 2.5
#'   (igraph::sample_fitness_pl), p edges, weights +-U[0.35, 0.5].
#' - Covariate networks (covariates with an effect): Erdos-Renyi with edge
#'   probability 1/(p - 1) (p/2 expected edges; 0.01 at p = 100 as in Zhang &
#'   Li, 2023), weights +-U[0.175, 0.25], zero diagonal.
#' - Common diagonal d: the smallest value keeping the minimum eigenvalue of
#'   Omega(u) = d I + B0 + sum_h u_h B_h at least eig_margin at the 2^q
#'   corners of the covariate box (exact, since the minimum eigenvalue is
#'   concave in u).
#' - Baseline target: as the template standardises continuous covariates and
#'   codes binary ones as 0/1 dummies, the estimated baseline is the network at
#'   the mean of the continuous covariates and at level 0 of the binary ones,
#'   T0 = Delta0 + sum_{h continuous} mean(x_h) B_h.
#' - Truth on the scale of the fit: the expression data are standardised
#'   (precision S Omega S), continuous covariates too (component x sd(x_h)).
#'
#' @param ncores Cores used to draw the subjects (one Cholesky per subject)
#' @return list(x = standardised expression matrix, covariates = data frame
#'   (binary covariates as factors with levels 0/1), T0 = baseline target on
#'   the generating scale (for the prior), truth = named list of true
#'   components on the scale of the fit, with the template names (Baseline,
#'   x1, x3_1, ...), info = data frame of the components, config)
simulate_dataset <- function(n, p, q, seed, eig_margin = 0.1, ncores = 1) {
  set.seed(seed)
  design <- covariate_design(q)
  weights <- function(adj, lo, hi) {
    w <- matrix(0, p, p)
    up <- which(upper.tri(adj) & adj)
    w[up] <- runif(length(up), lo, hi) * sample(c(-1, 1), length(up), replace = TRUE)
    w + t(w)
  }
  U <- sapply(seq_len(q), function(h) {
    if (design$type[h] == "binary") rbinom(n, 1, 0.5) else runif(n)
  })
  U <- matrix(U, n, q, dimnames = list(NULL, design$name))

  g <- sample_fitness_pl(p, no.of.edges = p, exponent.out = 2.5)
  B0 <- weights(as.matrix(as_adjacency_matrix(g)) > 0, 0.35, 0.5)
  Bh <- lapply(seq_len(q), function(h) {
    if (!design$effect[h]) return(matrix(0, p, p))
    adj <- matrix(runif(p * p) < 1 / (p - 1), p, p)
    adj[lower.tri(adj, diag = TRUE)] <- FALSE
    weights(adj | t(adj), 0.175, 0.25)
  })
  names(Bh) <- design$name

  corners <- as.matrix(expand.grid(rep(list(0:1), q)))
  min_eig <- min(apply(corners, 1, function(u) {
    min(eigen(B0 + Reduce(`+`, Map(`*`, Bh, u)), symmetric = TRUE, only.values = TRUE)$values)
  }))
  d <- eig_margin - min_eig
  Delta0 <- B0 + diag(d, p)

  # Subjects: x_i ~ N(0, Omega(u_i)^-1), drawn in chunks on ncores cores
  Zstd <- matrix(rnorm(n * p), n, p)
  chunks <- split(seq_len(n), cut(seq_len(n), max(1, ncores), labels = FALSE))
  draw <- function(idx) {
    t(sapply(idx, function(i) {
      omega <- Delta0 + Reduce(`+`, Map(`*`, Bh, U[i, ]))
      backsolve(chol(omega), Zstd[i, ])
    }))
  }
  x <- do.call(rbind, parallel::mclapply(chunks, draw, mc.cores = ncores))
  prot <- paste0("Prot", seq_len(p))
  colnames(x) <- prot

  continuous <- design$name[design$type == "continuous"]
  T0 <- Delta0
  for (h in continuous) T0 <- T0 + mean(U[, h]) * Bh[[h]]

  S <- diag(apply(x, 2, sd))
  on_fit_scale <- function(m) { m <- S %*% m %*% S; dimnames(m) <- list(prot, prot); m }
  truth <- list(Baseline = on_fit_scale(T0))
  info <- data.frame(component = "Baseline", covariate = NA, type = "baseline", effect = TRUE)
  for (h in seq_len(q)) {
    nm <- design$name[h]
    binary <- design$type[h] == "binary"
    comp <- if (binary) paste0(nm, "_1") else nm
    truth[[comp]] <- on_fit_scale(Bh[[nm]] * if (binary) 1 else sd(U[, nm]))
    info <- rbind(info, data.frame(component = comp, covariate = nm, type = design$type[h],
                                   effect = design$effect[h]))
  }

  covariates <- as.data.frame(U)
  for (h in design$name[design$type == "binary"]) {
    covariates[[h]] <- factor(covariates[[h]], levels = c(0, 1))
  }
  ut <- upper.tri(B0)
  list(x = scale(x), covariates = covariates, T0 = T0, truth = truth, info = info,
       config = list(n = n, p = p, q = q, seed = seed, diagonal = d,
                     baseline_edges = sum(B0[ut] != 0),
                     covariate_edges = sapply(Bh, function(m) sum(m[ut] != 0)),
                     target_edges = sum(abs(T0[ut]) > 1e-8)))
}

#' Prior knowledge from the support of the baseline target
#'
#' perfect: weight 1 on the edges of T0, 0 elsewhere.
#' noisy: edges of T0 kept with probability 0.8 (confidence U[0.7, 1]), non
#'   edges added with probability 0.1 (confidence U[0.1, 0.6]).
#' none: NULL (no prior; the template then uses penalty factor 1).
make_prior <- function(T0, type, seed, sensitivity = 0.8, specificity = 0.9,
                       confidence_range = c(0.7, 1)) {
  if (type == "none") return(NULL)
  p <- nrow(T0)
  edges <- abs(T0) > 1e-8
  diag(edges) <- FALSE
  W <- matrix(0, p, p, dimnames = dimnames(T0))
  if (type == "perfect") {
    W[edges] <- 1
    return(W)
  }
  set.seed(seed)
  ut <- which(upper.tri(W))
  true_e <- edges[ut]
  keep <- runif(length(ut)) < sensitivity
  false_pos <- runif(length(ut)) > specificity
  w <- numeric(length(ut))
  sel_t <- true_e & keep
  sel_f <- !true_e & false_pos
  w[sel_t] <- runif(sum(sel_t), confidence_range[1], confidence_range[2])
  w[sel_f] <- runif(sum(sel_f), 0.1, confidence_range[1] - 0.1)
  W[ut] <- w
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  W
}

#' Evaluate the estimated networks against the truth
#'
#' @param est Named list of estimated symmetric precision-scale matrices
#'   (Dic_Delta_hat of the template)
#' @param sim Output of simulate_dataset()
#' @return Data frame with one row per component
evaluate_networks <- function(est, sim) {
  out <- lapply(seq_len(nrow(sim$info)), function(r) {
    comp <- sim$info$component[r]
    Tm <- sim$truth[[comp]]; A <- as.matrix(est[[comp]])
    ut <- upper.tri(Tm)
    t_e <- abs(Tm[ut]) > 1e-8; e_e <- A[ut] != 0
    TP <- sum(t_e & e_e); FP <- sum(!t_e & e_e); FN <- sum(t_e & !e_e); TN <- sum(!t_e & !e_e)
    TPR <- if (TP + FN > 0) TP / (TP + FN) else NA
    PPV <- if (TP + FP > 0) TP / (TP + FP) else NA
    tp <- ut & A != 0 & abs(Tm) > 1e-8
    cor_tp <- function(method) if (sum(tp) > 2) suppressWarnings(cor(Tm[tp], A[tp], method = method)) else NA
    data.frame(
      component = comp, type = sim$info$type[r], effect = sim$info$effect[r],
      n_true_edges = sum(t_e), n_est_edges = sum(e_e), density = mean(e_e),
      TPR = TPR, FPR = FP / (FP + TN), Precision = PPV,
      F1 = if (is.na(TPR) || is.na(PPV) || TPR + PPV == 0) (if (sum(t_e) > 0) 0 else NA) else 2 * TPR * PPV / (TPR + PPV),
      Spearman_TP = cor_tp("spearman"), Pearson_TP = cor_tp("pearson"),
      Frobenius = sqrt(2 * sum((Tm[ut] - A[ut])^2)))
  })
  do.call(rbind, out)
}
