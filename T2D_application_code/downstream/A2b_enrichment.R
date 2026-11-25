#################################################
## A2b. Over-representation analysis (ORA) of the T2D OR communities and protein
## sets, with the 366 panel proteins as background (fixes bug 12).
## Gene sets from msigdbr (local, nothing sent out): KEGG Legacy, Reactome, GO BP.
## Run from the private project folder: Rscript --vanilla ../T2D_application_code/downstream/A2b_enrichment.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))
suppressPackageStartupMessages(library(msigdbr))

nodes_tab <- read.csv(file.path(OUT_DIR, "A2_node_table.csv"))
dbs <- list(KEGG = c("C2", "CP:KEGG_LEGACY"), Reactome = c("C2", "CP:REACTOME"), GO_BP = c("C5", "GO:BP"))
gs <- do.call(rbind, lapply(names(dbs), function(d) {
  x <- msigdbr(species = "Homo sapiens", collection = dbs[[d]][1], subcollection = dbs[[d]][2])
  data.frame(db = d, set = x$gs_name, gene = x$gene_symbol)
}))
universe <- intersect(prot, unique(gs$gene))
cat("Panel proteins found in the gene-set collections:", length(universe), "of", length(prot), "\n")
cat("Not found:", setdiff(prot, universe), "\n")
gs <- unique(gs[gs$gene %in% universe, ])
size <- table(paste(gs$db, gs$set))
gs <- gs[paste(gs$db, gs$set) %in% names(size)[size >= 3], ]       # >= 3 panel proteins in the set

ora <- function(query, label) {
  q <- intersect(query, universe); if (length(q) < 3) return(NULL)
  out <- do.call(rbind, lapply(split(gs, paste(gs$db, gs$set)), function(s) {
    k <- length(intersect(q, s$gene)); if (k < 2) return(NULL)
    M <- length(unique(s$gene))
    data.frame(query = label, query_size = length(q), db = s$db[1], set = s$set[1], set_size_in_panel = M,
               overlap = k, expected = length(q) * M / length(universe),
               p = phyper(k - 1, M, length(universe) - M, length(q), lower.tail = FALSE),
               genes = paste(sort(intersect(q, s$gene)), collapse = ","))
  }))
  if (is.null(out)) return(NULL)
  out$q_BH <- ave(out$p, out$db, FUN = function(p) p.adjust(p, "BH"))   # BH within query x database
  out
}
comm <- nodes_tab[!is.na(nodes_tab$community), ]
csize <- table(comm$community)
queries <- c(setNames(lapply(names(csize)[csize >= 5], function(c) comm$protein[comm$community == c]),
                      paste0("Community ", names(csize)[csize >= 5])),
             list(`Candidate hubs` = candidate_hubs, `Rewired proteins` = rewired,
                  `TIMP1 + partners` = c("TIMP1", nodes_tab$protein[nodes_tab$protein %in%
                     unlist(edge_list(A_or)[edge_list(A_or)$protein_1 == "TIMP1" | edge_list(A_or)$protein_2 == "TIMP1", 1:2])]),
                  `All connected` = comm$protein))
res <- do.call(rbind, Map(ora, queries, names(queries)))
res <- res[order(res$query, res$q_BH, res$p), ]
write.csv(res, file.path(OUT_DIR, "A2b_enrichment_ORA.csv"), row.names = FALSE)
cat("\nCommunities (>= 5 proteins):\n")
for (c in names(csize)[csize >= 5]) cat(sprintf("  Community %s (%d): %s\n", c, csize[[c]], paste(sort(comm$protein[comm$community == c]), collapse = ", ")))
cat("\nSignificant terms (BH q < 0.05 within query x database):\n")
sig <- res[res$q_BH < 0.05, ]
if (nrow(sig)) print(sig[, c("query", "db", "set", "overlap", "set_size_in_panel", "expected", "q_BH", "genes")], row.names = FALSE, digits = 2) else cat("  none\n")
cat("\nTop 3 per query (nominal):\n")
top <- do.call(rbind, lapply(split(res, res$query), function(x) head(x[order(x$p), ], 3)))
print(top[, c("query", "db", "set", "overlap", "set_size_in_panel", "p", "q_BH", "genes")], row.names = FALSE, digits = 2)
# T2D-specific pathways, whatever their significance
t2d_sets <- unique(gs$set[grepl("DIABETES|INSULIN|ADIPOCYTOKINE|PPAR|AGE_RAGE|GLUCOSE", gs$set)])
cat("\nDiabetes / insulin related sets touched by the queries:\n")
print(res[res$set %in% t2d_sets, c("query", "db", "set", "overlap", "set_size_in_panel", "p", "q_BH", "genes")], row.names = FALSE, digits = 2)
