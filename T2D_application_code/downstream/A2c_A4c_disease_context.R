#################################################
## A2 (DisGeNET part) and A4.3: disease context of the T2D OR graph.
##  - DisGeNET curated T2D genes (UMLS C0011860; academic tier = CURATED sources only),
##    downloaded with disgenet2r to results_bic_relaxed_10y/downstream/external/.
##  - Diabetes / insulin pathways from msigdbr (KEGG Legacy + Reactome).
## Protein level: enrichment of T2D genes among candidate hubs / rewired proteins
##   (hypergeometric vs the 366 panel) and degree-label association vs the
##   permutation null degrees.
## Edge level (A4.3): share of T2D edges touching / joining T2D genes, and share of
##   edges whose two proteins share a diabetes pathway, observed vs permutation null.
## Run from the private project folder: Rscript --vanilla ../T2D_application_code/downstream/A2c_A4c_disease_context.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))
suppressPackageStartupMessages(library(msigdbr))

dg <- read.csv(file.path(OUT_DIR, "external", "DisGeNET_T2D_C0011860_curated.csv"))
dg <- dg[dg$gene_symbol %in% prot, ]
labels <- list(`DisGeNET curated` = unique(dg$gene_symbol),
               `DisGeNET curated, score >= 0.6` = unique(dg$gene_symbol[dg$normalized_score >= 0.6]))
for (l in names(labels)) cat(sprintf("%s: %d panel proteins: %s\n", l, length(labels[[l]]), paste(sort(labels[[l]]), collapse = ", ")))

# ---- Protein level -----------------------------------------------------------
dn <- readRDS(file.path(PERM_DIR, "T2D_degree_obs_null.rds"))
prot_res <- do.call(rbind, lapply(names(labels), function(l) {
  L <- labels[[l]]
  hg <- function(set) { k <- sum(set %in% L)
    c(k = k, n = length(set), p = phyper(k - 1, length(L), length(prot) - length(L), length(set), lower.tail = FALSE)) }
  h1 <- hg(candidate_hubs); h2 <- hg(rewired); h3 <- hg(names(dn$obs)[dn$obs > 0])
  # degree-label association: mean T2D degree of labelled proteins minus others, vs null
  stat <- function(d) mean(d[names(d) %in% L]) - mean(d[!names(d) %in% L])
  obs <- stat(dn$obs); nul <- apply(dn$null, 1, function(d) stat(setNames(d, colnames(dn$null))))
  data.frame(label = l, n_label = length(L),
             hubs = sprintf("%d/%d", h1["k"], h1["n"]), p_hubs = h1["p"],
             rewired = sprintf("%d/%d", h2["k"], h2["n"]), p_rewired = h2["p"],
             connected = sprintf("%d/%d", h3["k"], h3["n"]), p_connected = h3["p"],
             degree_diff_obs = obs, degree_diff_null_mean = mean(nul),
             p_degree_perm = (1 + sum(nul >= obs)) / (length(nul) + 1),
             hubs_in_label = paste(intersect(rewired, L), collapse = ","))
}))
rownames(prot_res) <- NULL
cat("\nProtein level (panel share of label:", round(length(labels[[1]]) / length(prot), 3), "):\n")
print(prot_res, row.names = FALSE, digits = 3)
write.csv(prot_res, file.path(OUT_DIR, "A2c_DisGeNET_protein_level.csv"), row.names = FALSE)

# ---- Edge level (A4.3) ------------------------------------------------------------
pw <- rbind(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"),
            msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
pw <- pw[grepl("DIABETES|INSULIN|ADIPOCYTOKINE|PPAR|GLUCAGON|MATURITY_ONSET", pw$gs_name) &
         !grepl("INSULIN_LIKE_GROWTH|IGF", pw$gs_name) & pw$gene_symbol %in% prot, c("gs_name", "gene_symbol")]
pw <- unique(pw); pw_sizes <- table(pw$gs_name); pw <- pw[pw$gs_name %in% names(pw_sizes)[pw_sizes >= 2], ]
cat("\nDiabetes / insulin pathways with >= 2 panel proteins:\n"); print(table(pw$gs_name))
P <- matrix(FALSE, p, p, dimnames = list(prot, prot))
for (s in split(pw$gene_symbol, pw$gs_name)) P[s, s] <- TRUE
diag(P) <- FALSE
cat("Panel pairs sharing a diabetes pathway:", sum(P[upper.tri(P)]), "of", sum(upper.tri(P)), "\n")

L1 <- prot %in% labels[[1]]; L2 <- prot %in% labels[[2]]
edge_stats <- function(A) { up <- A & upper.tri(A); idx <- which(up, arr.ind = TRUE)
  c(n = nrow(idx),
    touch_DGN = mean(L1[idx[, 1]] | L1[idx[, 2]]), both_DGN = mean(L1[idx[, 1]] & L1[idx[, 2]]),
    touch_DGN06 = mean(L2[idx[, 1]] | L2[idx[, 2]]), both_DGN06 = mean(L2[idx[, 1]] & L2[idx[, 2]]),
    pathway = mean(P[idx])) }
o <- readRDS(file.path(PERM_DIR, "original.rds"))
fl <- sort(list.files(PERM_DIR, "^perm_[0-9]+_chunk"))
st <- do.call(rbind, lapply(fl, function(f) {
  x <- readRDS(file.path(PERM_DIR, f)); C <- matrix(0, p, ncol(o$coef)); C[x$nodes, ] <- as.matrix(x$coef)
  Bt <- C[, 1 + 3 * p + seq_len(p)]; diag(Bt) <- 0
  c(b = x$b, edge_stats((Bt != 0) | (t(Bt) != 0)))
}))
st <- as.data.frame(st); obs <- st[st$b == 0, ]; nul <- st[st$b > 0, ]
stopifnot(obs$n == 170)
panel <- edge_stats(matrix(TRUE, p, p))
edge_res <- do.call(rbind, lapply(c("touch_DGN", "both_DGN", "touch_DGN06", "both_DGN06", "pathway"), function(v)
  data.frame(statistic = v, observed = obs[[v]], observed_n = round(obs[[v]] * obs$n), null_mean = mean(nul[[v]]),
             null_q95 = quantile(nul[[v]], 0.95), panel_pairs = panel[[v]],
             p_perm = (1 + sum(nul[[v]] >= obs[[v]])) / (nrow(nul) + 1))))
rownames(edge_res) <- NULL
cat("\nEdge level (A4.3), share of the 170 T2D OR edges:\n"); print(edge_res, row.names = FALSE, digits = 3)
write.csv(edge_res, file.path(OUT_DIR, "A4c_disease_context_edges.csv"), row.names = FALSE)

edges <- edge_list(A_or)
dg_edges <- edges[edges$protein_1 %in% labels[[1]] & edges$protein_2 %in% labels[[1]] | P[cbind(edges$protein_1, edges$protein_2)], ]
cat("\nT2D edges joining two DisGeNET T2D genes or sharing a diabetes pathway:\n"); print(dg_edges[, 1:5], row.names = FALSE, digits = 3)
