#################################################
## A5. Edge table: all 170 T2D OR edges with weight, sign, direction of change,
## maxT p-value, STRING support (prior W and physical score), community,
## held-out replication, hub classes and DisGeNET flags. Descriptive only.
## Needs the outputs of A2, A2c, A3/A4 and A4d.
## Run from the private project folder: Rscript --vanilla ../T2D_application_code/downstream/A5_edge_table.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))

e <- read.csv(file.path(OUT_DIR, "A4_edges_with_replication_STRING.csv"))
ph <- read.csv(file.path(OUT_DIR, "A4d_edges_physical_score.csv"))
nodes <- read.csv(file.path(OUT_DIR, "A2_node_table.csv"))
dg <- unique(read.csv(file.path(OUT_DIR, "external", "DisGeNET_T2D_C0011860_curated.csv"))$gene_symbol)
stopifnot(nrow(e) == 170, identical(paste(e$protein_1, e$protein_2), paste(ph$protein_1, ph$protein_2)))

com <- setNames(nodes$community, nodes$protein); cls <- setNames(nodes$hub_class, nodes$protein)
tab <- data.frame(
  protein_1 = e$protein_1, protein_2 = e$protein_2,
  in_AND = e$in_AND,
  T2D_coefficient = e$weight, baseline_coefficient = e$baseline,
  dependence_in_cases = e$in_cases,
  p_maxT = e$p_maxT,
  STRING_prior_W = round(e$STRING_w, 4), STRING_combined_score = ifelse(e$STRING_w > 0, round(150 + 849 * e$STRING_w), 0),
  STRING_physical_score = ph$physical_score,
  community_1 = com[e$protein_1], community_2 = com[e$protein_2],
  same_community = com[e$protein_1] == com[e$protein_2],
  class_1 = cls[e$protein_1], class_2 = cls[e$protein_2],
  DisGeNET_T2D_1 = e$protein_1 %in% dg, DisGeNET_T2D_2 = e$protein_2 %in% dg,
  heldout_coefficient = e$beta_heldout, heldout_p = e$p_heldout, heldout_sign_concordant = e$concordant,
  replicated = e$concordant & e$p_heldout < 0.05)
tab <- tab[order(!tab$in_AND, !tab$replicated, tab$p_maxT, -abs(tab$T2D_coefficient)), ]
write.csv(tab, file.path(OUT_DIR, "A5_T2D_OR_edge_table.csv"), row.names = FALSE)

cat("Edges:", nrow(tab), "| AND:", sum(tab$in_AND), "| replicated (concordant, p < 0.05):", sum(tab$replicated),
    "| STRING (prior) support:", sum(tab$STRING_prior_W > 0), "| physical:", sum(tab$STRING_physical_score > 0),
    "| within community:", sum(tab$same_community), "\n")
cat("Direction of change: ", paste(names(table(tab$dependence_in_cases)), table(tab$dependence_in_cases), sep = " ", collapse = "; "), "\n")
cat("\nReplicated edges:\n")
print(tab[tab$replicated, c("protein_1", "protein_2", "in_AND", "T2D_coefficient", "dependence_in_cases", "heldout_coefficient", "heldout_p",
                            "STRING_combined_score", "class_1", "class_2")], row.names = FALSE, digits = 3)
cat("\nAND edges:\n")
print(tab[tab$in_AND, c("protein_1", "protein_2", "T2D_coefficient", "dependence_in_cases", "p_maxT", "heldout_coefficient", "heldout_p", "STRING_combined_score")],
      row.names = FALSE, digits = 3)
