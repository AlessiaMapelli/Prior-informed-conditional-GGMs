#################################################
## A4.2c Human Protein Atlas interaction consensus vs the permutation null.
## File: results_bic_relaxed_10y/downstream/external/interaction_consensus.tsv
## (proteinatlas.org/humanproteome/interaction; consensus of BioGRID, IntAct,
##  OpenCell, BioPlex and RoseTTAFold2-PPI predictions; Ensembl gene ids).
## Independent of the STRING prior. Symbols mapped to Ensembl genes via msigdbr.
## Two references: all consensus pairs; pairs with >= 2 experimental sources
## (RoseTTAFold2-PPI predictions not counted).
## Run from the private project folder: Rscript --vanilla ../T2D_application_code/downstream/A4e_HPA_interactions.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))
suppressPackageStartupMessages(library(msigdbr))

h <- read.delim(file.path(OUT_DIR, "external", "interaction_consensus.tsv"))
src <- strsplit(h$datasets, "_")
h$n_exp <- sapply(src, function(s) sum(s != "rosettafold2ppi"))
map <- unique(rbind(msigdbr(species = "Homo sapiens", collection = "C5")[, c("gene_symbol", "ensembl_gene")],
                    msigdbr(species = "Homo sapiens", collection = "C2")[, c("gene_symbol", "ensembl_gene")]))
map <- map[map$gene_symbol %in% prot & !is.na(map$ensembl_gene), ]
cat("Panel proteins with an Ensembl gene id:", length(unique(map$gene_symbol)), "of", p,
    "| unmapped:", setdiff(prot, map$gene_symbol), "\n")
h <- h[h$ensembl_gene_id_1 %in% map$ensembl_gene & h$ensembl_gene_id_2 %in% map$ensembl_gene, ]
build <- function(hh) {
  M <- matrix(FALSE, p, p, dimnames = list(prot, prot))
  for (i in seq_len(nrow(hh))) {
    a <- unique(map$gene_symbol[map$ensembl_gene == hh$ensembl_gene_id_1[i]])
    b <- unique(map$gene_symbol[map$ensembl_gene == hh$ensembl_gene_id_2[i]])
    M[a, b] <- TRUE; M[b, a] <- TRUE
  }
  diag(M) <- FALSE; M
}
refs <- list(`HPA consensus (all)` = build(h), `HPA consensus, >= 2 experimental sources` = build(h[h$n_exp >= 2, ]))
mapped <- prot %in% map$gene_symbol; up <- upper.tri(diag(p)) & outer(mapped, mapped)
for (r in names(refs)) cat(sprintf("%s: %d of %d mapped panel pairs (%.2f%%)\n", r, sum(refs[[r]][up]), sum(up), 100 * mean(refs[[r]][up])))

o <- readRDS(file.path(PERM_DIR, "original.rds"))
fl <- sort(list.files(PERM_DIR, "^perm_[0-9]+_chunk"))
st <- do.call(rbind, lapply(fl, function(f) {
  x <- readRDS(file.path(PERM_DIR, f)); C <- matrix(0, p, ncol(o$coef)); C[x$nodes, ] <- as.matrix(x$coef)
  Bt <- C[, 1 + 3 * p + seq_len(p)]; diag(Bt) <- 0; A <- ((Bt != 0) | (t(Bt) != 0)) & up
  c(b = x$b, n = sum(A), sapply(refs, function(R) sum(R[A])))
}))
st <- as.data.frame(st, check.names = FALSE); obs <- st[st$b == 0, ]; nul <- st[st$b > 0, ]
res <- do.call(rbind, lapply(names(refs), function(r) data.frame(
  reference = r, observed = obs[[r]], observed_share = obs[[r]] / obs$n,
  null_mean_share = mean(nul[[r]] / nul$n), null_q95_share = quantile(nul[[r]] / nul$n, 0.95),
  panel_share = mean(refs[[r]][up]),
  p_perm = (1 + sum(nul[[r]] / nul$n >= obs[[r]] / obs$n)) / (nrow(nul) + 1))))
rownames(res) <- NULL
cat("\nT2D OR edges among mapped proteins:", obs$n, "\n"); print(res, row.names = FALSE, digits = 3)
write.csv(res, file.path(OUT_DIR, "A4e_HPA_interactions_vs_null.csv"), row.names = FALSE)

edges <- edge_list(A_or)
edges$HPA <- refs[[1]][cbind(edges$protein_1, edges$protein_2)]
edges$HPA_exp2 <- refs[[2]][cbind(edges$protein_1, edges$protein_2)]
cat("\nT2D OR edges in the HPA consensus:\n"); print(edges[edges$HPA, c("protein_1", "protein_2", "weight", "in_AND", "HPA_exp2")], row.names = FALSE, digits = 3)
write.csv(edges[, c("protein_1", "protein_2", "HPA", "HPA_exp2")], file.path(OUT_DIR, "A4e_edges_HPA.csv"), row.names = FALSE)
