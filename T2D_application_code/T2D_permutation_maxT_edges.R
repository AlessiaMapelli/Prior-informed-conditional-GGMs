# Edge-level maxT test of the T2D OR edges against the permutation null
# (family-wise error over all edges). Run after T2D_permutation_null.R collect:
#   Rscript --vanilla ../T2D_application_code/T2D_permutation_maxT_edges.R
suppressMessages(library(Matrix))
D <- "results_bic_relaxed_10y/permutation_null"; o <- readRDS(file.path(D, "original.rds")); p <- nrow(o$coef)
nm <- read.csv(file.path(D, "T2D_null_degree_by_protein.csv"))$protein
wmat <- function(C) { B <- C[, 1 + 3 * p + seq_len(p)]; diag(B) <- 0; A <- abs(B); At <- t(A)
  W <- ifelse(A <= At, B, t(B)); W[!(B != 0 & t(B) != 0)] <- 0; W }
mins <- numeric(0); maxw <- numeric(0); fl <- sort(list.files(D, "^perm_[0-9]+_chunk"))
for (f in fl) { x <- readRDS(file.path(D, f)); C <- matrix(0, p, ncol(o$coef)); C[x$nodes, ] <- as.matrix(x$coef)
  W <- wmat(C); if (x$b == 0) { W0 <- W; m0 <- x$minutes } else { maxw <- c(maxw, max(abs(W[upper.tri(W)]))); mins <- c(mins, x$minutes) } }
idx <- which(W0 != 0 & upper.tri(W0), arr.ind = TRUE)
cat("minutes per permutation: perm0", round(m0), "| null median", round(median(mins)), "max", round(max(mins)), "\n")
cat("null max |w| quantiles:", round(quantile(maxw, c(.5, .9, .95, 1)), 3), "| perms with no AND edge:", sum(maxw == 0), "\n")
for (k in seq_len(nrow(idx))) { w <- W0[idx[k,1], idx[k,2]]
  cat(sprintf("%s-%s  w=%.3f  FWER p (maxT) = %.3f\n", nm[idx[k,1]], nm[idx[k,2]], w, (1 + sum(maxw >= abs(w))) / (length(maxw) + 1))) }
# OR rule: edge weight = larger |B| of the two directions; maxT over the union graph
wor <- function(C) { B <- C[, 1 + 3 * p + seq_len(p)]; diag(B) <- 0; pmax(abs(B), abs(t(B))) }
maxo <- numeric(0)
for (f in fl) { x <- readRDS(file.path(D, f)); C <- matrix(0, p, ncol(o$coef)); C[x$nodes, ] <- as.matrix(x$coef)
  W <- wor(C); if (x$b == 0) Wo0 <- W else maxo <- c(maxo, max(W[upper.tri(W)])) }
io <- which(Wo0 != 0 & upper.tri(Wo0), arr.ind = TRUE)
po <- sapply(seq_len(nrow(io)), function(k) (1 + sum(maxo >= Wo0[io[k,1], io[k,2]])) / (length(maxo) + 1))
ro <- data.frame(p1 = nm[io[,1]], p2 = nm[io[,2]], w = Wo0[io], p_maxT = po); ro <- ro[order(ro$p_maxT, -ro$w), ]
cat("\nOR edges:", nrow(ro), "| null max |w| quantiles:", round(quantile(maxo, c(.5,.9,.95,1)),3), "\n")
cat("OR edges with FWER p <= 0.05:", sum(ro$p_maxT <= .05), "| <= 0.10:", sum(ro$p_maxT <= .10), "\n"); print(head(ro, 15), row.names = FALSE)
write.csv(ro, file.path(D, "T2D_OR_edges_maxT.csv"), row.names = FALSE)
