#################################################
## Exploratory analysis of the lambda selection rules
#################################################
# Node-wise regressions of the prior-informed conditional GGM, fitted once per
# node and hyperparameter combination and then compared under several
# selection rules:
#   min     cross-validation, lambda.min (current default)
#   1se     cross-validation, lambda.1se
#   bic     BIC along the lambda path (no cross-validation)
#   ebic    extended BIC along the lambda path (no cross-validation)
#   stab    stability selection (Meinshausen & Buhlmann, 2010) at the
#           hyperparameters chosen by "min"
# Each rule is evaluated after an OLS refit on the selected predictors.
# Post-selection thresholds (disparity-filter backbone, partial correlation
# threshold) are applied to the lambda.min networks in evaluate_setting().
# (Version 3, 100 replications; the v2 version, which also kept the penalised
# coefficients and the CV diagnostics, is in archive_v2_code/.)
#
# The design matrix, penalty factors, hyperparameter grid and random search
# mirror GGReg_cov_single_node_processing() in code/ggReg_main_functions.R.

suppressPackageStartupMessages({
  library(sparsegl)
  library(igraph)
  library(MASS)
  library(matrixcalc)
  library(reshape2)
  library(ggplot2)
})

#' Load the functions of the original code needed here, without running the
#' library() calls of the original files (some packages are not needed).
#'
#' @param code_dir Directory with the copies of Sim_functions.R and
#'   ggReg_main_functions.R
#' @return Invisibly, the names of the loaded functions
load_original_functions <- function(code_dir) {
  wanted <- list(
    "Sim_functions.R" = c("generate_network_structure", "generate_mean_effects",
                          "generate_prior_knowledge", "simulate_subject_data_robust",
                          "generate_simulation_dataset", "calculate_binary_metrics",
                          "calculate_frobenius_error", "calculate_magnitude_preserved",
                          "evaluate_estimation"),
    "ggReg_main_functions.R" = c("build_delta_dictionary"))
  for (f in names(wanted)) {
    for (e in parse(file.path(code_dir, f))) {
      if (is.call(e) && identical(e[[1]], as.name("<-")) &&
          as.character(e[[2]]) %in% wanted[[f]]) {
        eval(e, envir = globalenv())
      }
    }
  }
  invisible(unlist(wanted))
}

#################################################
## DATA
#################################################

#' Simulate one dataset from the conditional GGM
#'
#' Follows the design of Zhang & Li (2023, Section 6), with the precision
#' matrix positive definite for every covariate value by construction:
#' - covariates: x1, x2 ~ Uniform[0, 1] (continuous), x3 ~ Bernoulli(0.5)
#'   (binary); two of the three, chosen at random, have an effect on the
#'   network;
#' - population network: power-law degree distribution with exponent 2.5
#'   (igraph::sample_fitness_pl) and p edges, weights +-U[0.35, 0.5];
#' - covariate networks: Erdos-Renyi with edge probability 0.01, weights
#'   +-U[0.175, 0.25], zero diagonal;
#' - diagonal: one constant d for all nodes and all covariate values, the
#'   smallest that keeps Omega(u) = d I + B0 + sum_h u_h B_h positive definite
#'   (minimum eigenvalue >= eig_margin). Omega(u) is linear in u, so its minimum
#'   eigenvalue is concave in u and attains its minimum at a corner of the
#'   covariate box: checking the 2^q corners is exact.
#' The conditional variances 1/d are therefore free of the covariates and the
#' node-wise regression coefficients are linear in u, as the model assumes.
#'
#' Target of the baseline. As in the template and in the package, continuous
#' covariates are standardised before fitting and binary ones enter as 0/1
#' dummies (reference level 0). The estimated baseline is therefore the network
#' at the mean of the continuous covariates and at the reference level of the
#' binary ones:
#'   T0 = Delta0 + sum_{h continuous} mean(x_h) B_h.
#' The prior is generated from the support of T0, so that it describes the
#' network the baseline estimates.
#'
#' True components on the scale of the fit. The expression data are
#' standardised (Z = X S^-1, S = diag of the marginal standard deviations), so
#' the precision of Z is S Omega(u) S; continuous covariates are standardised
#' too, so their component is multiplied by sd(x_h). The returned `truth` is on
#' this scale and is directly comparable with the estimates; `truth_raw` keeps
#' the generating matrices.
#'
#' @param n Number of samples
#' @param p Number of nodes
#' @param prior_type "perfect", "noisy" or "none"
#' @param seed Random seed
#' @param out_dir Unused, kept for the interface of run_setting.R
#' @param n_edges_baseline Number of edges of the population network
#' @param edge_prob_cov Edge probability of the covariate networks
#' @param eig_margin Minimum eigenvalue of Omega(u) over the covariate range
#' @return A list with the data (pat_id, x1-x3, Prot1-Protp standardised),
#'   the true components on the fit scale (`truth`: Baseline = S T0 S,
#'   covariates), the generating matrices (`truth_raw`: Delta0, T0, B_h), the
#'   prior knowledge matrix, d and the configuration
simulate_setting <- function(n, p, prior_type, seed, out_dir = NULL,
                             n_edges_baseline = p, edge_prob_cov = 0.01,
                             eig_margin = 0.1) {
  set.seed(seed)
  weights <- function(adj, lo, hi) {
    w <- matrix(0, p, p)
    up <- which(upper.tri(adj) & adj)
    w[up] <- runif(length(up), lo, hi) * sample(c(-1, 1), length(up), replace = TRUE)
    w + t(w)
  }
  cov <- data.frame(x1 = runif(n), x2 = runif(n), x3 = rbinom(n, 1, 0.5))
  continuous <- c("x1", "x2"); binary <- "x3"
  q <- ncol(cov)

  g <- sample_fitness_pl(p, no.of.edges = n_edges_baseline, exponent.out = 2.5)
  B0 <- weights(as.matrix(as_adjacency_matrix(g)) > 0, 0.35, 0.5)
  effective <- sample(seq_len(q), 2)
  Bh <- lapply(seq_len(q), function(h) {
    if (!(h %in% effective)) return(matrix(0, p, p))
    adj <- matrix(runif(p * p) < edge_prob_cov, p, p)
    adj[lower.tri(adj, diag = TRUE)] <- FALSE
    weights(adj | t(adj), 0.175, 0.25)
  })
  names(Bh) <- names(cov)

  corners <- as.matrix(expand.grid(rep(list(0:1), q)))
  min_eig <- min(apply(corners, 1, function(u) {
    min(eigen(B0 + Reduce(`+`, Map(`*`, Bh, u)), symmetric = TRUE, only.values = TRUE)$values)
  }))
  d <- eig_margin - min_eig
  Delta0 <- B0 + diag(d, p)

  x <- matrix(0, n, p)
  U <- as.matrix(cov)
  for (i in seq_len(n)) {
    omega <- Delta0 + Reduce(`+`, Map(`*`, Bh, U[i, ]))
    x[i, ] <- backsolve(chol(omega), rnorm(p))
  }
  prot <- paste0("Prot", seq_len(p))
  colnames(x) <- prot
  data <- cbind(data.frame(pat_id = seq_len(n)), cov, as.data.frame(scale(x)))

  # Baseline target: mean of the continuous covariates, reference level of
  # the binary ones
  T0 <- Delta0
  for (h in continuous) T0 <- T0 + mean(cov[[h]]) * Bh[[h]]

  # True components on the scale of the fit
  S <- diag(apply(x, 2, sd))
  on_fit_scale <- function(m) { m <- S %*% m %*% S; dimnames(m) <- list(prot, prot); m }
  truth <- list(Baseline = on_fit_scale(T0))
  for (h in names(cov)) {
    truth[[h]] <- on_fit_scale(Bh[[h]] * if (h %in% continuous) sd(cov[[h]]) else 1)
  }
  truth_raw <- c(list(Delta0 = Delta0, T0 = T0), Bh)

  prior <- generate_prior_knowledge(
    true_adjacency = T0, prior_type = prior_type,
    noise_params = list(sensitivity = 0.8, specificity = 0.9, confidence_range = c(0.7, 1)))
  if (sum(prior > 0) == 0) prior <- NULL

  list(data = data, truth = truth, truth_raw = truth_raw, prior_knowledge = prior,
       diagonal = d,
       config = list(n_samples = n, n_nodes = p, prior_type = prior_type,
                     n_edges_baseline = n_edges_baseline, edge_prob_cov = edge_prob_cov,
                     continuous = continuous, binary = binary,
                     effective = names(cov)[effective], min_eig = eig_margin))
}

#' Build the inputs shared by all node regressions
#'
#' As in the template, continuous covariates are standardised and binary ones
#' are kept as 0/1 dummies.
#'
#' @param sim_dataset Output of simulate_setting()
#' @param nfolds Number of cross-validation folds
#' @param alpha_grid,theta_grid Hyperparameter grids
#' @param p_random Fraction of the grid explored per node (random search)
#' @param seed Seed for the fold assignment
#' @return A list with the data, the design blocks and the settings
prepare_inputs <- function(sim_dataset, nfolds = 5,
                           alpha_grid = c(0.5, 0.75, 0.9, 0.99),
                           theta_grid = c(0.8, 1.0, 1.1, 1.3),
                           p_random = 0.5, seed = 2024) {
  d <- sim_dataset$data
  Z <- as.matrix(d[, grepl("^Prot", colnames(d))])
  U <- as.matrix(d[, grepl("^x", colnames(d)), drop = FALSE])
  cont <- sim_dataset$config$continuous
  U[, cont] <- scale(U[, cont])
  W <- sim_dataset$prior_knowledge
  if (is.null(W)) theta_grid <- 1
  set.seed(seed)
  foldid <- sample(rep(seq_len(nfolds), length.out = nrow(Z)))
  list(Z = Z, U = U, W = W, n = nrow(Z), p = ncol(Z), q = ncol(U),
       cov_names = colnames(U), foldid = foldid, nfolds = nfolds,
       alpha_grid = alpha_grid, theta_grid = theta_grid, p_random = p_random)
}


#################################################
## NODE-WISE EXPLORATION
#################################################

#' Columns of the template layout (intercept, then one block of p columns per
#' component) holding the predictors of node i
node_layout_columns <- function(i, p, q) {
  others <- setdiff(seq_len(p), i)
  unlist(lapply(0:q, function(b) 1 + b * p + others))
}

#' OLS refit on a set of predictors (no intercept, as in the penalised fit)
refit_ols <- function(X, y, support) {
  coef <- numeric(ncol(X))
  if (any(support)) {
    fit <- lm.fit(X[, support, drop = FALSE], y)
    b <- fit$coefficients
    b[is.na(b)] <- 0
    coef[support] <- b
    rss <- sum(fit$residuals^2)
  } else {
    rss <- sum(y^2)
  }
  list(coef = coef, sigma2 = rss / max(length(y) - sum(support), 1))
}


#' Fit node i under all selection rules
#'
#' @param i Node index
#' @param inp Output of prepare_inputs()
#' @param stab_B Number of subsamples for stability selection (0 to skip)
#' @param stab_frac Subsample fraction
#' @param ebic_gamma EBIC parameter
#' @param cv Run the cross-validation (rules min and 1se). With FALSE only the
#'   lambda paths are fitted (rules bic and ebic)
#' @param stab_at Rule whose hyperparameters (alpha, theta) are used for
#'   stability selection: "min" (needs cv = TRUE), "bic" or "bic_relaxed"
#'   (needs relaxed = TRUE)
#' @param relaxed Also compute the relaxed criteria (rules bic_relaxed and
#'   ebic_relaxed): BIC / EBIC on the RSS of the OLS refit of each support on
#'   the path, from the Gram matrix X'X computed once per node
#' @return A list with, per rule, the OLS-refitted coefficients and residual
#'   variance on the selected predictors, the stability selection counts and
#'   the computing times
explore_node <- function(i, inp, stab_B = 50, stab_frac = 0.5, ebic_gamma = 0.5,
                         cv = TRUE, stab_at = "min", relaxed = FALSE) {
  stopifnot(stab_at %in% c("min", "bic", "bic_relaxed"), cv || stab_at != "min",
            relaxed || stab_at != "bic_relaxed")
  Z <- inp$Z; U <- inp$U; n <- inp$n; p <- inp$p; q <- inp$q
  others <- setdiff(seq_len(p), i)
  y <- Z[, i]
  Zo <- Z[, others, drop = FALSE]
  X <- do.call(cbind, c(list(Zo), lapply(seq_len(q), function(h) U[, h] * Zo)))
  group <- rep(seq_len(q + 1), each = p - 1)
  P <- ncol(X)
  pf_base <- if (is.null(inp$W)) rep(1, p - 1) else 1 - inp$W[i, others]

  combos <- expand.grid(alpha = inp$alpha_grid, theta = inp$theta_grid)
  if (nrow(combos) > 1 && inp$p_random < 1) {
    set.seed(234 + 10 * i)  # same per-node random search as the template
    combos <- combos[sample(seq_len(nrow(combos)),
                            size = floor(inp$p_random * nrow(combos))), , drop = FALSE]
  }

  fits <- vector("list", nrow(combos))
  time_path <- time_cv <- time_relaxed <- 0
  if (relaxed) {
    # Relaxed RSS of a support s: y'y - b'X_s'y with b = (X_s'X_s)^-1 X_s'y
    # (QR fallback if X_s'X_s is singular)
    t0 <- Sys.time()
    XtX <- crossprod(X); Xty <- drop(crossprod(X, y)); yty <- sum(y^2)
    relaxed_rss <- function(s) {
      if (!any(s)) return(yty)
      R <- tryCatch(chol(XtX[s, s, drop = FALSE]), error = function(e) NULL)
      if (is.null(R)) return(sum(lm.fit(X[, s, drop = FALSE], y)$residuals^2))
      b <- backsolve(R, forwardsolve(t(R), Xty[s]))
      max(yty - sum(Xty[s] * b), 0)
    }
    time_relaxed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  for (k in seq_len(nrow(combos))) {
    pf_sparse <- c(pf_base, rep(combos$theta[k], q * (p - 1)))
    args <- list(x = X, y = y, group = group, pf_group = c(0, rep(1, q)),
                 pf_sparse = pf_sparse, asparse = combos$alpha[k], intercept = FALSE)
    t0 <- Sys.time()
    path <- do.call(sparsegl, args)
    time_path <- time_path + as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cvfit <- NULL
    if (cv) {
      t0 <- Sys.time()
      cvfit <- do.call(cv.sparsegl, c(args, list(lambda = path$lambda, foldid = inp$foldid)))
      time_cv <- time_cv + as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    B <- as.matrix(path$beta)
    rss <- colSums((y - X %*% B)^2)
    df <- colSums(B != 0)
    bic <- n * log(rss / n) + log(n) * df
    fits[[k]] <- list(
      args = args, lambda = path$lambda, B = B, bic = bic,
      ebic = bic + 2 * ebic_gamma * lchoose(P, df),
      idx_min = if (cv) which(path$lambda == cvfit$lambda.min)[1] else NA,
      idx_1se = if (cv) which(path$lambda == cvfit$lambda.1se)[1] else NA,
      L = length(path$lambda))
    if (relaxed) {
      t0 <- Sys.time()
      S <- B != 0
      sk <- apply(S, 2, function(s) paste(which(s), collapse = ","))
      uk <- unique(sk)
      ru <- vapply(uk, function(u) relaxed_rss(S[, match(u, sk)]), numeric(1))
      bic_rel <- n * log(unname(ru[match(sk, uk)]) / n) + log(n) * df
      fits[[k]]$bic_relaxed <- bic_rel
      fits[[k]]$ebic_relaxed <- bic_rel + 2 * ebic_gamma * lchoose(P, df)
      time_relaxed <- time_relaxed + as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
  }

  pick <- function(k, idx) {
    f <- fits[[k]]
    list(combo = k, idx = idx, coef = f$B[, idx], alpha = combos$alpha[k],
         theta = combos$theta[k], lambda = f$lambda[idx])
  }
  # min / 1se: lambda by CV within each combination, combination by BIC
  # (as in the template); bic / ebic: joint minimum over combinations and path
  by_cv <- function(field) {
    bics <- sapply(fits, function(f) f$bic[f[[field]]])
    k <- which.min(bics)
    pick(k, fits[[k]][[field]])
  }
  by_ic <- function(field) {
    best <- sapply(fits, function(f) min(f[[field]]))
    k <- which.min(best)
    pick(k, which.min(fits[[k]][[field]]))
  }
  sel <- list(bic = by_ic("bic"), ebic = by_ic("ebic"))
  if (cv) sel <- c(list(min = by_cv("idx_min"), `1se` = by_cv("idx_1se")), sel)
  if (relaxed) sel <- c(sel, list(bic_relaxed = by_ic("bic_relaxed"),
                                  ebic_relaxed = by_ic("ebic_relaxed")))

  # Stability selection at the hyperparameters (alpha, theta) chosen by the
  # rule stab_at, on the lambda path of the full-data fit
  stab <- NULL
  time_stab <- 0
  if (stab_B > 0) {
    f <- fits[[sel[[stab_at]]$combo]]
    counts <- matrix(0L, P, f$L)
    m <- floor(stab_frac * n)
    set.seed(5000 + i)
    t0 <- Sys.time()
    for (b in seq_len(stab_B)) {
      s <- sample.int(n, m)
      a <- f$args; a$x <- X[s, , drop = FALSE]; a$y <- y[s]; a$lambda <- f$lambda
      fb <- do.call(sparsegl, a)
      Lb <- ncol(fb$beta)
      counts[, seq_len(Lb)] <- counts[, seq_len(Lb)] + (as.matrix(fb$beta) != 0)
    }
    time_stab <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    stab <- list(counts = counts, B = stab_B, lambda = f$lambda, at = stab_at)
  }

  # OLS refit on the predictors selected by each rule
  t0 <- Sys.time()
  rules <- lapply(sel, function(s) {
    c(s[c("alpha", "theta", "lambda", "idx")],
      list(refit = refit_ols(X, y, s$coef != 0)))
  })
  time_refit <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / length(rules)

  list(node = i, P = P, rules = rules, stab = stab,
       time = c(path = time_path, cv = time_cv, stab = time_stab, refit = time_refit,
                relaxed = time_relaxed))
}

#################################################
## COLLECTION AND EVALUATION
#################################################

#' Stability selection decision from the stored counts
#'
#' @param stab The stab element of a node result
#' @param pi_thr Selection-probability threshold
#' @param pfer Per-node bound on the expected number of false selections,
#'   which sets the lambda region (Meinshausen & Buhlmann, 2010, Theorem 1)
#' @return Logical vector of selected predictors
stability_support <- function(stab, pi_thr = 0.75, pfer = 1) {
  P <- nrow(stab$counts)
  freq <- stab$counts / stab$B
  q_lambda <- sqrt(pfer * (2 * pi_thr - 1) * P)
  # Lambda region: from the largest lambda down to the last value at which the
  # average number of selected predictors is still within the bound
  avg_sel <- colSums(freq)
  first_over <- which(avg_sel > q_lambda)[1]
  region <- if (is.na(first_over)) seq_along(avg_sel) else seq_len(first_over - 1)
  if (length(region) == 0) return(rep(FALSE, P))
  apply(freq[, region, drop = FALSE], 1, max) >= pi_thr
}

#' Assemble node results into template-layout matrices and symmetrise them
#'
#' @param nodes List of explore_node() results, ordered by node
#' @param inp Output of prepare_inputs()
#' @param coef_fun Function(node_result) returning list(coef, sigma2)
#' @param symm "OR" or "AND"
#' @return Named list (Baseline and covariate names) of symmetric matrices of
#'   precision-scale coefficients (delta), plus the regression-scale ones
assemble_networks <- function(nodes, inp, coef_fun, symm) {
  p <- inp$p; q <- inp$q
  ncol_all <- 1 + (q + 1) * p
  Bm <- Dm <- matrix(0, p, ncol_all)
  for (r in nodes) {
    cs <- coef_fun(r)
    cols <- node_layout_columns(r$node, p, q)
    Bm[r$node, cols] <- cs$coef
    Dm[r$node, cols] <- -cs$coef / cs$sigma2
  }
  groups <- c(0, rep(seq_len(q + 1), each = p))
  U_names <- c("Intercept", "Baseline", inp$cov_names)
  list(delta = build_delta_dictionary(Dm, groups, U_names, symm, ref = Bm),
       coef = build_delta_dictionary(Bm, groups, U_names, symm))
}

#' Disparity-filter backbone (Serrano et al., 2009) of a weighted network, as
#' in the T2D report (alpha = 0.05, Benjamini-Hochberg). An edge is kept when
#' it is significant for at least one of its two endpoints.
disparity_backbone <- function(A, alpha = 0.05) {
  W <- abs(A); diag(W) <- 0
  k <- rowSums(W > 0); s <- rowSums(W)
  pv <- matrix(1, nrow(W), ncol(W))
  for (i in which(k > 1)) {
    j <- which(W[i, ] > 0)
    pv[i, j] <- (1 - W[i, j] / s[i])^(k[i] - 1)
  }
  pmin_edge <- pmin(pv, t(pv))
  ut <- upper.tri(W) & W > 0
  keep <- matrix(FALSE, nrow(W), ncol(W))
  if (any(ut)) keep[ut] <- p.adjust(pmin_edge[ut], method = "BH") < alpha
  keep <- keep | t(keep)
  A * keep
}

#' Partial-correlation scale of a symmetric delta matrix, given the node
#' residual variances: rho_jk = -delta_jk * sigma_j * sigma_k
partial_correlation <- function(D, sigma2) {
  -D * sqrt(outer(sigma2, sigma2))
}

#' Evaluate all rules on one simulated dataset
#'
#' All rules use the OLS-refitted coefficients. The post-selection thresholds
#' (backbone, partial correlation) are applied to the lambda.min networks.
#' Edges are compared with the true components on the scale of the fit
#' (sim_dataset$truth; the baseline target is the network at the mean of the
#' continuous covariates and the reference level of the binary ones).
#'
#' @param nodes List of explore_node() results
#' @param inp Output of prepare_inputs()
#' @param sim_dataset Output of simulate_setting()
#' @param pcor_thr Threshold on the absolute partial correlation
#' @return A data frame with one row per rule, post-processing,
#'   symmetrisation and component
evaluate_setting <- function(nodes, inp, sim_dataset, pcor_thr = 0.05) {
  nodes <- nodes[order(sapply(nodes, `[[`, "node"))]
  truth <- sim_dataset$truth
  rule_names <- names(nodes[[1]]$rules)

  variants <- list()
  for (rn in rule_names) {
    variants[[rn]] <- local({ rn <- rn; function(r) r$rules[[rn]]$refit })
  }
  if (!is.null(nodes[[1]]$stab)) {
    # OLS refit on the stable set needs the data, so it is done here
    stab_refit <- lapply(nodes, function(r) {
      i <- r$node; p <- inp$p; q <- inp$q
      Zo <- inp$Z[, setdiff(seq_len(p), i), drop = FALSE]
      X <- do.call(cbind, c(list(Zo), lapply(seq_len(q), function(h) inp$U[, h] * Zo)))
      refit_ols(X, inp$Z[, i], stability_support(r$stab))
    })
    names(stab_refit) <- sapply(nodes, `[[`, "node")
    stab_name <- switch(nodes[[1]]$stab$at, min = "stab", bic = "stab_bic",
                        bic_relaxed = "stab_bic_relaxed")
    variants[[stab_name]] <- function(r) stab_refit[[as.character(r$node)]]
  }

  if ("min" %in% rule_names) sigma2_min <- sapply(nodes, function(r) r$rules$min$refit$sigma2)
  out <- list()
  for (vn in names(variants)) for (symm in c("OR", "AND")) {
    nets <- assemble_networks(nodes, inp, variants[[vn]], symm)$delta
    post <- list(none = nets)
    if (vn == "min") {
      post$backbone <- lapply(nets, disparity_backbone)
      post$pcor <- lapply(nets, function(D) D * (abs(partial_correlation(D, sigma2_min)) > pcor_thr))
    }
    for (pp in names(post)) for (comp in names(truth)) {
      A <- post[[pp]][[comp]]; Tm <- truth[[comp]]
      ut <- upper.tri(Tm)
      m <- calculate_binary_metrics(Tm, A)
      # Magnitude on the true-positive edges
      tp <- ut & A != 0 & abs(Tm) > 1e-8
      cor_tp <- function(method) if (sum(tp) > 2) suppressWarnings(cor(Tm[tp], A[tp], method = method)) else NA
      out[[length(out) + 1]] <- data.frame(
        rule = vn, post = pp, symm = symm, component = comp,
        component_type = if (comp == "Baseline") "baseline" else
          if (any(abs(Tm[ut]) > 1e-8)) "effect" else "null",
        n_true_edges = sum(abs(Tm[ut]) > 1e-8), n_est_edges = sum(A[ut] != 0),
        density = mean(A[ut] != 0), TPR = m$TPR, FPR = m$FPR,
        Precision = m$Precision, F1 = m$F1,
        Spearman_TP = cor_tp("spearman"), Pearson_TP = cor_tp("pearson"),
        # off-diagonal entries only
        Frobenius = sqrt(2 * sum((Tm[ut] - A[ut])^2)))
    }
  }
  do.call(rbind, out)
}
