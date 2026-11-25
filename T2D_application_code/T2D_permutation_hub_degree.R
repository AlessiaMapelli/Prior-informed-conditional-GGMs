# Hub test: degree of each protein in the T2D OR graph against its permutation
# null (permutation_null/ of T2D_permutation_null.R). Run from the private
# project folder after T2D_permutation_null.R collect:
#   Rscript --vanilla ../T2D_application_code/T2D_permutation_hub_degree.R
D <- "results_bic_relaxed_10y/permutation_null"; o <- readRDS(file.path(D, "original.rds")); p <- nrow(o$coef)
nm <- read.csv(file.path(D, "T2D_null_degree_by_protein.csv"))$protein
deg <- function(C) { B <- C[, 1 + 3 * p + seq_len(p)]; diag(B) <- 0; A <- (B != 0) | (t(B) != 0); rowSums(A) }
fl <- sort(list.files(D, "^perm_[0-9]+_chunk")); Dn <- NULL
for (f in fl) { x <- readRDS(file.path(D, f)); C <- matrix(0, p, ncol(o$coef)); C[x$nodes, ] <- as.matrix(x$coef)
  d <- deg(C); if (x$b == 0) d0 <- d else Dn <- rbind(Dn, d) }
maxd <- apply(Dn, 1, max)
cat("null max OR T2D degree quantiles:", quantile(maxd, c(.5,.9,.95,1)), "\n")
r <- data.frame(protein = nm, obs = d0, null_mean = colMeans(Dn), null_q95 = apply(Dn, 2, quantile, .95),
  p_marg = (1 + colSums(Dn >= rep(d0, each = nrow(Dn)))) / (nrow(Dn) + 1),
  p_maxT = sapply(d0, function(v) (1 + sum(maxd >= v)) / (length(maxd) + 1)))
r <- r[order(-r$obs), ]; print(head(r, 12), row.names = FALSE, digits = 3)
write.csv(r, file.path(D, "T2D_hub_degree_test.csv"), row.names = FALSE)
