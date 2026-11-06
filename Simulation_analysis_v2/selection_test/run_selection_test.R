#################################################
## Selection-rule test: BIC on the penalised vs the relaxed (OLS refit) RSS,
## current vs extended hyperparameter grid, stability selection at the
## hyperparameters of either BIC
#################################################
# Usage (from the Simulation_analysis_v2 folder):
#   Rscript selection_test/run_selection_test.R --n 20000 --p 100 --q 3 --rep 1 [--ncores 16]
#           [--stab 1] [--out selection_test/results]
# --stab 0 skips the two stability-selection variants (Part 2 of the
# sensitivity analysis: BIC / EBIC relaxed, current vs extended grid).
# Same datasets as the main study (dataset_seed), fitted under the three
# priors. For each node the lambda path of every hyperparameter combination is
# fitted once (union of the current and the extended grid) and all rules are
# applied to the same paths. All rules report the OLS refit on the selected
# support. Variants (rule x grid):
#   bic           BIC on the RSS of the penalised fit (current rule)
#   bic_relaxed   BIC on the RSS of the OLS refit of each support on the path
#   ebic_relaxed  EBIC (gamma 0.5) on the relaxed RSS
#   stab_bic          stability selection at the combination chosen by bic
#   stab_bic_relaxed  stability selection at the combination chosen by bic_relaxed
# Grids (50% random search, the same sampling as the template, seed 234 + 10 i):
#   current   alpha {0.5, 0.75, 0.9, 0.99}       x theta {0.8, 1, 1.1, 1.3}
#   extended  alpha {0.5, 0.75, 0.9, 0.95, 0.99, 1} x theta {0.4, 0.6, 0.8, 1, 1.3}
# Stability selection is run on the current grid only.
# Check: bic on the current grid must reproduce results/<setting>/rep<r>/evaluation.rds.
# Outputs in selection_test/results/n<n>_p<p>_q<q>/rep<rep>/:
#   evaluation.rds  metrics per prior, variant, symmetrisation and component
#   timing.rds      per prior and node: seconds for each path, each relaxed
#                   RSS computation and each stability run
#   choices.rds     per prior, node and variant: alpha, theta, lambda, df

opt <- list(n = NA, p = NA, q = NA, rep = NA, ncores = NA, out = "selection_test/results", stab = 1)
args <- commandArgs(trailingOnly = TRUE)
for (k in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[k])]] <- args[k + 1]
for (v in c("n", "p", "q", "rep", "ncores", "stab")) opt[[v]] <- as.integer(opt[[v]])
if (is.na(opt$ncores)) opt$ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores()))

suppressPackageStartupMessages({library(sparsegl); library(glmnet)})
source("ggReg_main_functions.R")
source("R/sim_design.R")

task <- sprintf("n%d_p%d_q%d", opt$n, opt$p, opt$q)
out_dir <- file.path(opt$out, task, paste0("rep", opt$rep))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out_dir, "evaluation.rds"))) {
  cat("Already done:", out_dir, "\n"); quit(save = "no")
}
seed <- dataset_seed(opt$n, opt$p, opt$q, opt$rep)
sim <- simulate_dataset(opt$n, opt$p, opt$q, seed, ncores = opt$ncores)
n <- opt$n; p <- opt$p
Z0 <- sim$x

# Covariates as in the template: continuous scaled, factors one-hot (reference 0)
cov <- sim$covariates
num <- sapply(cov, is.numeric)
cov[, num] <- scale(cov[, num])
C <- model.matrix(~ ., cov)
cn <- colnames(C)
for (v in names(cov)) if (is.factor(cov[[v]]))
  for (l in levels(cov[[v]])) cn <- sub(paste0("^", v, l, "$"), paste0(v, "_", l), cn)
colnames(C) <- cn
q <- ncol(C) - 1
groups_all <- c(0, rep(seq_len(q + 1), each = p))
U_names <- c("Intercept", "Baseline", cn[-1])

grid_current <- expand.grid(asparse = c(0.5, 0.75, 0.9, 0.99), weight = c(0.8, 1, 1.1, 1.3))
grid_extended <- expand.grid(asparse = c(0.5, 0.75, 0.9, 0.95, 0.99, 1),
                             weight = c(0.4, 0.6, 0.8, 1, 1.3))
sample_grid <- function(g, i) {
  set.seed(234 + 10 * i)
  g[sample(seq_len(nrow(g)), size = max(1, floor(0.5 * nrow(g)))), , drop = FALSE]
}
ctrl <- list(ebic_gamma = 0.5, stab_B = 50, stab_fraction = 0.5, stab_pi = 0.75, stab_pfer = 1)

ols_rss <- function(X, y, s) {
  if (!any(s)) return(list(rss = sum(y^2), coef = numeric(length(s))))
  fit <- lm.fit(X[, s, drop = FALSE], y)
  b <- fit$coefficients; b[is.na(b)] <- 0
  coef <- numeric(length(s)); coef[s] <- b
  list(rss = sum(fit$residuals^2), coef = coef)
}

stability <- function(f, X, y, i) {
  P <- ncol(X); counts <- matrix(0, P, length(f$lambda))
  m <- floor(ctrl$stab_fraction * n)
  set.seed(5000 + i)
  for (b in seq_len(ctrl$stab_B)) {
    s <- sample.int(n, m)
    a <- f$args; a$x <- X[s, , drop = FALSE]; a$y <- y[s]; a$lambda <- f$lambda
    fb <- do.call(sparsegl, a)
    Lb <- ncol(fb$beta)
    counts[, seq_len(Lb)] <- counts[, seq_len(Lb)] + (as.matrix(fb$beta) != 0)
  }
  freq <- counts / ctrl$stab_B
  q_lambda <- sqrt(ctrl$stab_pfer * (2 * ctrl$stab_pi - 1) * P)
  first_over <- which(colSums(freq) > q_lambda)[1]
  region <- if (is.na(first_over)) seq_along(f$lambda) else seq_len(first_over - 1)
  if (length(region) == 0) rep(FALSE, P) else apply(freq[, region, drop = FALSE], 1, max) >= ctrl$stab_pi
}

fit_node <- function(i, pf_baseline_all) {
  y <- Z0[, i]
  X <- do.call(cbind, lapply(seq_len(q + 1), function(j) Z0[, -i] * C[, j]))
  grp <- rep(seq_len(q + 1), each = p - 1)
  pf_group <- c(0, rep(1, q))
  pf_baseline <- pf_baseline_all[i, -i]
  make_pf <- function(w) { v <- rep(w, ncol(X)); v[grp == 1] <- pf_baseline; v }
  P <- ncol(X)
  t_gram <- system.time({XtX <- crossprod(X); Xty <- drop(crossprod(X, y)); yty <- sum(y^2)})[["elapsed"]]
  relaxed_rss <- function(s) {
    if (!any(s)) return(yty)
    R <- tryCatch(chol(XtX[s, s, drop = FALSE]), error = function(e) NULL)
    if (is.null(R)) return(ols_rss(X, y, s)$rss)  # rank deficient: QR fallback
    b <- backsolve(R, forwardsolve(t(R), Xty[s]))
    max(yty - sum(Xty[s] * b), 0)
  }
  cur <- sample_grid(grid_current, i); ext <- sample_grid(grid_extended, i)
  key <- function(g) paste(g$asparse, g$weight)
  all <- unique(rbind(cur, ext))
  fits <- lapply(seq_len(nrow(all)), function(k) {
    args <- list(x = X, y = y, group = grp, pf_group = pf_group,
                 pf_sparse = make_pf(all$weight[k]), asparse = all$asparse[k], intercept = FALSE)
    t_path <- system.time(path <- do.call(sparsegl, args))[["elapsed"]]
    B <- as.matrix(path$beta); S <- B != 0; df <- colSums(S)
    rss_pen <- colSums((y - X %*% B)^2)
    # Relaxed RSS of each distinct support from the Gram matrix (computed once
    # per node, outside this timing): RSS = y'y - b' X_s'y, b = (X_s'X_s)^-1 X_s'y
    t_rel <- system.time({
      sk <- apply(S, 2, function(s) paste(which(s), collapse = ","))
      uk <- unique(sk)
      ru <- sapply(uk, function(u) relaxed_rss(S[, match(u, sk)]))
      rss_rel <- unname(ru[match(sk, uk)])
    })[["elapsed"]]
    bic <- function(rss) n * log(rss / n) + log(n) * df
    list(args = args, lambda = path$lambda, S = S, df = df,
         bic = bic(rss_pen), bic_relaxed = bic(rss_rel),
         ebic_relaxed = bic(rss_rel) + 2 * ctrl$ebic_gamma * lchoose(P, df),
         t_path = t_path, t_rel = t_rel, n_supports = length(uk))
  })
  names(fits) <- key(all)
  pick <- function(grid, field) {
    fs <- fits[key(grid)]
    k <- which.min(sapply(fs, function(f) min(f[[field]])))
    f <- fs[[k]]; idx <- which.min(f[[field]])
    list(f = f, combo = grid[k, ], support = f$S[, idx], lambda = f$lambda[idx])
  }
  variants <- list(
    bic_current = pick(cur, "bic"), bic_relaxed_current = pick(cur, "bic_relaxed"),
    ebic_relaxed_current = pick(cur, "ebic_relaxed"),
    bic_extended = pick(ext, "bic"), bic_relaxed_extended = pick(ext, "bic_relaxed"),
    ebic_relaxed_extended = pick(ext, "ebic_relaxed"))
  t_stab <- c()
  for (base in if (opt$stab == 1) c("bic_current", "bic_relaxed_current") else character(0)) {
    v <- variants[[base]]
    nm <- sub("_current", "", paste0("stab_", base))
    same <- base == "bic_relaxed_current" &&
      identical(key(v$combo), key(variants$bic_current$combo))
    if (same) {
      variants[[nm]] <- variants$stab_bic; t_stab[nm] <- t_stab[["stab_bic"]]; next
    }
    t_stab[nm] <- system.time(s <- stability(v$f, X, y, i))[["elapsed"]]
    variants[[nm]] <- list(f = v$f, combo = v$combo, support = s, lambda = NA)
  }
  rows <- lapply(variants, function(v) {
    r <- ols_rss(X, y, v$support)
    sigma2 <- r$rss / max(n - sum(v$support), 1)
    list(coef = r$coef, sigma2 = sigma2)
  })
  timing <- data.frame(node = i, combo = names(fits),
                       in_current = names(fits) %in% key(cur), in_extended = names(fits) %in% key(ext),
                       t_path = sapply(fits, `[[`, "t_path"), t_relaxed = sapply(fits, `[[`, "t_rel"),
                       n_supports = sapply(fits, `[[`, "n_supports"), t_gram = t_gram)
  choices <- do.call(rbind, lapply(names(variants), function(nm) {
    v <- variants[[nm]]
    data.frame(node = i, variant = nm, asparse = v$combo$asparse, weight = v$combo$weight,
               lambda = v$lambda, df = sum(v$support))
  }))
  list(rows = rows, timing = timing,
       t_stab = data.frame(node = rep(i, length(t_stab)), variant = names(t_stab), t = unname(t_stab)),
       choices = choices)
}

evaluation <- list(); timing <- list(); choices <- list(); stab_t <- list()
for (pr in c("perfect", "noisy", "none")) {
  known_ppi <- make_prior(sim$T0, pr, seed + 7)
  pf_all <- if (is.null(known_ppi)) matrix(1, p, p) else 1 - known_ppi
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(p), fit_node, pf_baseline_all = pf_all,
                            mc.cores = opt$ncores, mc.preschedule = FALSE)
  bad <- sapply(res, inherits, "try-error")
  if (any(bad)) { print(res[bad][[1]]); stop("Node failures") }
  cat(sprintf("Prior %-7s: %.1f min\n", pr, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  for (vn in names(res[[1]]$rows)) {
    Prec <- matrix(0, p, 1 + (q + 1) * p); Delta <- Prec
    for (i in seq_len(p)) {
      cols <- setdiff(seq_len(ncol(Prec)), 1 + (0:q) * p + i)[-1]
      Prec[i, cols] <- res[[i]]$rows[[vn]]$coef
      Delta[i, cols] <- -res[[i]]$rows[[vn]]$coef / res[[i]]$rows[[vn]]$sigma2
    }
    for (symm in c("OR", "AND")) {
      est <- build_delta_dictionary(Delta, groups_all, U_names, symm, ref = Prec)
      evaluation[[length(evaluation) + 1]] <-
        cbind(prior = pr, variant = vn, symm = symm, evaluate_networks(est, sim))
    }
  }
  timing[[pr]] <- cbind(prior = pr, do.call(rbind, lapply(res, `[[`, "timing")))
  st <- do.call(rbind, lapply(res, `[[`, "t_stab"))
  if (NROW(st) > 0) stab_t[[pr]] <- cbind(prior = pr, st)
  choices[[pr]] <- cbind(prior = pr, do.call(rbind, lapply(res, `[[`, "choices")))
}
meta <- data.frame(n = n, p = p, q = opt$q, rep = opt$rep)
evaluation <- cbind(meta, do.call(rbind, evaluation))
saveRDS(evaluation, file.path(out_dir, "evaluation.rds"))
saveRDS(list(paths = cbind(meta, do.call(rbind, timing)),
             stability = if (length(stab_t) > 0) cbind(meta, do.call(rbind, stab_t)) else NULL),
        file.path(out_dir, "timing.rds"))
saveRDS(cbind(meta, do.call(rbind, choices)), file.path(out_dir, "choices.rds"))

# Check against the main study
main_file <- file.path("results", task, paste0("rep", opt$rep), "evaluation.rds")
if (file.exists(main_file)) {
  main <- readRDS(main_file)
  mine <- evaluation[evaluation$variant == "bic_current", ]
  m <- merge(main, mine, by = c("prior", "symm", "component"))
  cat("bic_current reproduces the main study (F1, ", nrow(m), " rows): ",
      isTRUE(all.equal(m$F1.x, m$F1.y)), "\n", sep = "")
}
cat("Done:", out_dir, "\n")
