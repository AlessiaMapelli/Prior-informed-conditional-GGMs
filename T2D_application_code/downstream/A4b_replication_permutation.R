#################################################
## A4.1b Permutation calibration of the held-out replication.
## The binomial sign test assumes independent edges; edges share nodes (e.g. TIMP1).
## Here the T2D labels of the held-out data are permuted within sex x age-quintile
## strata, the residuals and designs are rebuilt, the same supports refitted, and
## the edge-level statistics recomputed. b = 0 is the observed labelling.
## Run from the private project folder (see README.md):
##   Rscript --vanilla downstream/A4b_replication_permutation.R [B]
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))
suppressPackageStartupMessages(use_prcggm())
args <- commandArgs(trailingOnly = TRUE); B <- if (length(args)) as.integer(args[1]) else 100

ho <- rbind(train_data, test_data)
strata <- interaction(ho$Sex, cut(ho$Age, quantile(ho$Age, 0:5 / 5), include.lowest = TRUE))
edges <- edge_list(A_or)
dir_terms <- which(B_t2d != 0, arr.ind = TRUE)
nodes <- sort(unique(dir_terms[, 1]))
supports <- lapply(nodes, function(j) { cf <- fit$coefficients[j, ]; setdiff(names(cf)[cf != 0], "(Intercept)") })
# For each edge, the directed term that defines the OR weight
edge_dir <- t(mapply(function(a, b) {
  i <- match(a, prot); k <- match(b, prot)
  if (abs(B_t2d[i, k]) >= abs(B_t2d[k, i])) c(i, k) else c(k, i)
}, edges$protein_1, edges$protein_2))
sign_sel <- sign(B_t2d[edge_dir])

stats_for <- function(diab) {
  d <- ho; d$Diab <- diab
  U <- prCGGM:::covariate_design(code_covariates(d), fit$covariates)
  R <- as.matrix(d[, prot]) - cbind(1, U) %*% t(fit$mean_effects[prot, ])
  est <- matrix(NA_real_, p, p); pv <- matrix(NA_real_, p, p)
  for (a in seq_along(nodes)) {
    j <- nodes[a]; cols <- supports[[a]]
    X <- sapply(cols, function(cn) { s <- strsplit(cn, ":", fixed = TRUE)[[1]]
      if (length(s) == 1) R[, s] else U[, s[1]] * R[, s[2]] })
    X <- cbind(1, X, U); y <- R[, j]
    qx <- qr(X); cf <- qr.coef(qx, y); res <- y - X %*% ifelse(is.na(cf), 0, cf)
    df <- nrow(X) - qx$rank; s2 <- sum(res^2) / df
    XtXi <- chol2inv(qr.R(qx)); se <- sqrt(diag(XtXi) * s2)
    piv <- qx$pivot[seq_len(qx$rank)]; names(se) <- colnames(X)[piv]
    ks <- dir_terms[dir_terms[, 1] == j, 2]
    for (k in ks) { tn <- paste0("Diab_Incident:", prot[k]); b <- cf[tn]
      est[j, k] <- b; pv[j, k] <- 2 * pt(-abs(b / se[tn]), df) }
  }
  e <- est[edge_dir]; pe <- pv[edge_dir]
  c(concordant = sum(sign(e) == sign_sel), p05 = sum(pe < 0.05), p05_conc = sum(pe < 0.05 & sign(e) == sign_sel))
}

set.seed(2026)
res <- matrix(NA, B + 1, 3, dimnames = list(NULL, c("concordant", "p05", "p05_conc")))
t0 <- Sys.time()
for (b in 0:B) {
  diab <- if (b == 0) ho$Diab else { x <- ho$Diab; for (s in levels(strata)) { i <- which(strata == s); x[i] <- x[i][sample.int(length(i))] }; x }
  res[b + 1, ] <- stats_for(diab)
  if (b %in% c(0, 1, 10) || b %% 25 == 0) cat("b =", b, ":", res[b + 1, ], "|", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
}
res <- data.frame(b = 0:B, res); write.csv(res, file.path(OUT_DIR, "A4_replication_permutation.csv"), row.names = FALSE)
obs <- res[1, ]; nul <- res[-1, ]
for (v in c("concordant", "p05", "p05_conc"))
  cat(sprintf("%-11s observed %3d | null mean %.1f, 95%% %.1f, max %d | permutation p = %.3f\n", v, obs[[v]],
              mean(nul[[v]]), quantile(nul[[v]], 0.95), max(nul[[v]]), (1 + sum(nul[[v]] >= obs[[v]])) / (B + 1)))
