#################################################
## A4.2b STRING physical interactions vs the permutation null.
## Physical subnetwork of STRING v12.5 (9606.protein.physical.links.v12.5.txt.gz,
## downloaded by the user from string-db.org into results_bic_relaxed_10y/downstream/external/).
## ENSP ids mapped to the Olink target names with STRING v12.0 protein.info (all v12.5
## physical ids are present there). The prior of the fit (STRING_W, v12.0 combined score,
## all channels) is not used here.
## Run from the private project folder: Rscript --vanilla ../T2D_application_code/downstream/A4d_STRING_physical.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))

info <- read.delim(file.path(STRING_DIR, "9606.protein.info.v12.0.txt"),
                   header = FALSE, quote = "")[, 1:2]
names(info) <- c("ensp", "name"); info <- info[info$name %in% prot, ]
ph <- read.table(gzfile(file.path(OUT_DIR, "external", "9606.protein.physical.links.v12.5.txt.gz")), header = TRUE)
ph <- ph[ph$protein1 %in% info$ensp & ph$protein2 %in% info$ensp, ]
ph$a <- info$name[match(ph$protein1, info$ensp)]; ph$b <- info$name[match(ph$protein2, info$ensp)]
Sp <- matrix(0, p, p, dimnames = list(prot, prot))
for (i in seq_len(nrow(ph))) Sp[ph$a[i], ph$b[i]] <- max(Sp[ph$a[i], ph$b[i]], ph$combined_score[i])
Sp <- pmax(Sp, t(Sp)); diag(Sp) <- 0
mapped <- prot %in% info$name
cat("Panel proteins mapped to STRING:", sum(mapped), "of", p, "| unmapped:", prot[!mapped], "\n")
up <- upper.tri(Sp) & outer(mapped, mapped)
thr <- c(any = 150, medium = 400, high = 700)
for (t in names(thr)) cat(sprintf("Panel pairs with physical score >= %d (%s): %d of %d (%.2f%%)\n", thr[t], t,
                                  sum(Sp[up] >= thr[t]), sum(up), 100 * mean(Sp[up] >= thr[t])))

stats <- function(Bt) {
  diag(Bt) <- 0; A <- ((Bt != 0) | (t(Bt) != 0)) & up
  c(n = sum(A), sapply(thr, function(t) mean(Sp[A] >= t)), mean_score = mean(Sp[A]))
}
o <- readRDS(file.path(PERM_DIR, "original.rds"))
fl <- sort(list.files(PERM_DIR, "^perm_[0-9]+_chunk"))
st <- do.call(rbind, lapply(fl, function(f) {
  x <- readRDS(file.path(PERM_DIR, f)); C <- matrix(0, p, ncol(o$coef)); C[x$nodes, ] <- as.matrix(x$coef)
  c(b = x$b, stats(C[, 1 + 3 * p + seq_len(p)]))
}))
st <- as.data.frame(st); obs <- st[st$b == 0, ]; nul <- st[st$b > 0, ]
res <- do.call(rbind, lapply(c(names(thr), "mean_score"), function(v)
  data.frame(statistic = if (v == "mean_score") "mean physical score" else sprintf("share with physical score >= %d (%s)", thr[v], v),
             observed = obs[[v]], observed_n = if (v == "mean_score") NA else round(obs[[v]] * obs$n),
             null_mean = mean(nul[[v]]), null_q95 = quantile(nul[[v]], 0.95),
             panel_pairs = if (v == "mean_score") mean(Sp[up]) else mean(Sp[up] >= thr[v]),
             p_perm = (1 + sum(nul[[v]] >= obs[[v]])) / (nrow(nul) + 1))))
rownames(res) <- NULL
cat("\nT2D OR edges among mapped proteins:", obs$n, "\n"); print(res, row.names = FALSE, digits = 3)
write.csv(res, file.path(OUT_DIR, "A4d_STRING_physical_vs_null.csv"), row.names = FALSE)
write.csv(st, file.path(OUT_DIR, "A4d_STRING_physical_by_permutation.csv"), row.names = FALSE)

# Edges with physical support, and TIMP1 against its own physical degree
edges <- edge_list(A_or)
edges$physical_score <- Sp[cbind(edges$protein_1, edges$protein_2)]
cat("\nT2D OR edges with physical support (score >= 150):\n")
print(edges[edges$physical_score > 0, c("protein_1", "protein_2", "weight", "in_AND", "physical_score")][order(-edges$physical_score[edges$physical_score > 0]), ],
      row.names = FALSE, digits = 3)
tp <- Sp["TIMP1", setdiff(prot[mapped], "TIMP1")]
part <- setdiff(unlist(edges[edges$protein_1 == "TIMP1" | edges$protein_2 == "TIMP1", 1:2]), "TIMP1")
for (t in names(thr)) { K <- sum(tp >= thr[t]); k <- sum(tp[part] >= thr[t])
  cat(sprintf("TIMP1 (%s): partners with physical support %d/%d; TIMP1 physical degree in panel %d/%d; hypergeometric p = %.3g\n",
              t, k, length(part), K, length(tp), phyper(k - 1, K, length(tp) - K, length(part), lower.tail = FALSE))) }
write.csv(edges[, c("protein_1", "protein_2", "physical_score")], file.path(OUT_DIR, "A4d_edges_physical_score.csv"), row.names = FALSE)
