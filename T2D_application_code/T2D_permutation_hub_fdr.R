# FDR of the hub classes defined by the hub test (BH on the marginal p-values,
# permutation FDR by degree threshold). Run after T2D_permutation_hub_degree.R:
#   Rscript --vanilla ../T2D_application_code/T2D_permutation_hub_fdr.R
D <- "results_bic_relaxed_10y/permutation_null"; o <- readRDS(file.path(D, "original.rds")); p <- nrow(o$coef)
nm <- read.csv(file.path(D, "T2D_null_degree_by_protein.csv"))$protein
deg <- function(C) { B <- C[, 1 + 3 * p + seq_len(p)]; diag(B) <- 0; rowSums((B != 0) | (t(B) != 0)) }
fl <- sort(list.files(D, "^perm_[0-9]+_chunk")); Dn <- NULL
for (f in fl) { x <- readRDS(file.path(D, f)); C <- matrix(0, p, ncol(o$coef)); C[x$nodes, ] <- as.matrix(x$coef)
  d <- deg(C); if (x$b == 0) d0 <- d else Dn <- rbind(Dn, d) }
saveRDS(list(obs = setNames(d0, nm), null = `colnames<-`(Dn, nm)), file.path(D, "T2D_degree_obs_null.rds"))
B <- nrow(Dn)
# 1. BH on marginal p
pm <- (1 + colSums(Dn >= rep(d0, each = B))) / (B + 1)
qbh <- p.adjust(pm, "BH")
cat("Marginal p at floor (0.0099):", sum(pm <= 1/(B+1) + 1e-12), "| BH q <= 0.05:", sum(qbh <= .05), "| <= 0.10:", sum(qbh <= .10), "| <= 0.20:", sum(qbh <= .20), "\n")
cat("min BH q:", round(min(qbh), 3), "\n")
# 2. permutation FDR by degree threshold (pooled over proteins)
cat("\nDegree threshold: observed proteins, expected under null, FDR\n")
for (k in 3:10) { obs <- sum(d0 >= k); exp <- mean(rowSums(Dn >= k)); cat(sprintf("deg >= %2d: obs %3d  null %.2f  FDR %.3f\n", k, obs, exp, min(1, exp / obs))) }
# 3. empirical FDR of the rule (marginal p <= a & degree >= k), null p from the other permutations
rule_fdr <- function(a, k) {
  nf <- sapply(seq_len(B), function(b) { pb <- (1 + colSums(Dn[-b, , drop = FALSE] >= rep(Dn[b, ], each = B - 1))) / B
    sum(pb <= a & Dn[b, ] >= k) })
  obs <- sum(pm <= a & d0 >= k); c(obs = obs, null_mean = mean(nf), null_q95 = unname(quantile(nf, .95)), FDR = min(1, mean(nf) / obs)) }
cat("\nRule (marginal p <= a & degree >= k): observed, mean null, FDR\n")
for (a in c(0.05, 0.01)) for (k in c(3, 4, 5, 6)) { r <- rule_fdr(a, k); cat(sprintf("p <= %.2f, deg >= %d: obs %2d  null mean %.2f (95%% %g)  FDR %.3f\n", a, k, r["obs"], r["null_mean"], r["null_q95"], r["FDR"])) }
