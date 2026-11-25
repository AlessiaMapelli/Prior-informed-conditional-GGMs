#################################################
## Draft main and supplementary figures of the 10-year T2D case study.
## Reads the outputs of the downstream scripts (results_bic_relaxed_10y/downstream/)
## and writes figures to results_bic_relaxed_10y/downstream/figures/.
## Run from the private project folder: Rscript --vanilla ../T2D_application_code/downstream/Figures.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))
suppressPackageStartupMessages({ library(igraph); library(ggraph); library(ggrepel); library(cowplot) })
FIG_DIR <- file.path(OUT_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)
save_f <- function(g, name, w, h) { for (ext in c("pdf", "png")) ggsave(file.path(FIG_DIR, paste0(name, ".", ext)), g, width = w, height = h, dpi = 300, bg = "white"); cat("Written", name, "\n") }
theme_set(theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")))
col_obs <- "#d95f02"; col_null <- "grey70"
cls_cols <- c(`Confirmed hub` = "#d95f02", `Candidate hub` = "#7570b3", Rewired = "#1b9e77", Other = "grey75")
chg_cols <- c(stronger = "#b2182b", weaker = "#2166ac", `case-specific` = "grey30")

null_panel <- function(nul, obs, xlab, title, lab = NULL) {
  p <- (1 + sum(nul >= obs)) / (length(nul) + 1)
  ggplot(data.frame(x = nul), aes(x)) + geom_histogram(bins = 25, fill = col_null, colour = "white") +
    geom_vline(xintercept = obs, colour = col_obs, linewidth = 0.9) +
    annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = 2.4, label.size = 0,
             label = paste0(if (is.null(lab)) sprintf("observed %s", format(round(obs, 3))) else lab, "\npermutation p = ", formatC(p, format = "g", digits = 2))) +
    labs(x = xlab, y = "Permutations", title = title) + theme(plot.title = element_text(size = 8, face = "bold"))
}

# ---- Figure 1 (main): validation of the T2D component ---------------------------
ec <- read.csv(file.path(PERM_DIR, "edge_counts_by_permutation.csv")); ec <- ec[ec$component == "Diab_Incident", ]
dd <- read.csv(file.path(OUT_DIR, "A2_degree_distribution.csv"))
rp <- read.csv(file.path(OUT_DIR, "A4_replication_permutation.csv"))
sv <- read.csv(file.path(OUT_DIR, "A4_STRING_vs_null.csv"))
dn <- readRDS(file.path(PERM_DIR, "T2D_degree_obs_null.rds"))
dgl <- unique(read.csv(file.path(OUT_DIR, "external", "DisGeNET_T2D_C0011860_curated.csv"))$gene_symbol)
dstat <- function(d) mean(d[names(d) %in% dgl]) - mean(d[!names(d) %in% dgl])
d_null <- apply(dn$null, 1, function(d) dstat(setNames(d, colnames(dn$null))))

pa <- null_panel(ec$OR[ec$perm > 0], ec$OR[ec$perm == 0], "T2D-specific edges (OR)", "a  Global test: number of T2D edges", sprintf("observed %d (null mean %.0f)", ec$OR[ec$perm == 0], mean(ec$OR[ec$perm > 0])))
tail_df <- rbind(data.frame(k = dd$degree, n = dd$obs_ge, what = "Observed"), data.frame(k = dd$degree, n = dd$null_ge, what = "Permutation null (mean)"))
tail_df <- tail_df[tail_df$k >= 1 & tail_df$n > 0, ]
pb <- ggplot(tail_df, aes(k, n, colour = what)) + geom_step(linewidth = 0.7) + scale_y_log10() +
  scale_colour_manual(values = c(Observed = col_obs, `Permutation null (mean)` = "grey40"), name = NULL) +
  annotate("text", x = 24, y = 1.4, label = "TIMP1", size = 2.4, hjust = 1) +
  labs(x = "T2D degree k", y = "Proteins with degree >= k", title = "b  Degree distribution") +
  theme(legend.position = c(0.68, 0.85), legend.background = element_blank(), plot.title = element_text(size = 8, face = "bold"))
pc <- null_panel(rp$concordant[-1], rp$concordant[1], "Sign-concordant edges (of 170)", "c  Held-out replication: sign concordance", sprintf("observed %d / 170", rp$concordant[1]))
pd <- null_panel(rp$p05_conc[-1], rp$p05_conc[1], "Concordant edges with p < 0.05", "d  Held-out replication: replicated edges", sprintf("observed %d", rp$p05_conc[1]))
pe <- null_panel(sv$OR_meanW[sv$b > 0], sv$OR_meanW[sv$b == 0], "Mean STRING weight of T2D edges", "e  STRING support vs null")
pf <- null_panel(d_null, dstat(dn$obs), "Mean T2D degree: DisGeNET T2D genes - others", "f  T2D genes are more rewired", sprintf("observed %+.2f", dstat(dn$obs)))
save_f(plot_grid(pa, pb, pc, pd, pe, pf, ncol = 3, align = "hv"), "Fig1_T2D_validation", 9, 5.6)

# ---- Figure 2 (main): the T2D OR network ------------------------------------------
et <- read.csv(file.path(OUT_DIR, "A5_T2D_OR_edge_table.csv"))
nt <- read.csv(file.path(OUT_DIR, "A2_node_table.csv")); nt <- nt[nt$degree > 0, ]
csz <- table(nt$community); big <- names(sort(csz, decreasing = TRUE))[1:min(9, sum(csz >= 5))]
big <- as.character(sort(as.integer(big)))
nt$comm <- factor(ifelse(as.character(nt$community) %in% big, paste("C", nt$community), "other (< 5)"),
                  levels = c(paste("C", big), "other (< 5)"))
nt$cls <- ifelse(nt$hub_class == "", "Other", nt$hub_class)
g <- graph_from_data_frame(data.frame(from = et$protein_1, to = et$protein_2, w = abs(et$T2D_coefficient),
                                      change = et$dependence_in_cases, rep = et$replicated, and = et$in_AND),
                           vertices = nt[, c("protein", "degree", "comm", "cls")], directed = FALSE)
set.seed(7)
lay <- create_layout(g, layout = "fr", niter = 2000)
comm_cols <- setNames(c(RColorBrewer::brewer.pal(12, "Set3")[c(1:8, 10, 11)][seq_along(big)], "grey90"), c(paste("C", big), "other (< 5)"))
p2 <- ggraph(lay) +
  geom_edge_link(aes(edge_colour = change, edge_width = rep, edge_alpha = rep)) +
  scale_edge_colour_manual(values = chg_cols, name = "Dependence in cases") +
  scale_edge_width_manual(values = c(`FALSE` = 0.3, `TRUE` = 1.1), labels = c("not replicated", "replicated (held-out)"), name = NULL) +
  scale_edge_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 1), guide = "none") +
  geom_node_point(aes(size = degree, fill = comm, colour = cls), shape = 21, stroke = 0.9) +
  scale_fill_manual(values = comm_cols, name = "Community") +
  scale_colour_manual(values = c(`Confirmed hub` = "#d95f02", `Candidate hub` = "#7570b3", Rewired = "#1b9e77", Other = "grey55"), name = "Protein class") +
  scale_size_continuous(range = c(1.2, 7), breaks = c(1, 5, 10, 24), name = "T2D degree") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3.5, colour = "grey40")),
         colour = guide_legend(override.aes = list(shape = 21, size = 3.5, fill = "white", stroke = 1.2))) +
  geom_node_text(aes(label = ifelse(cls != "Other" | degree >= 6, name, "")), repel = TRUE, size = 2.4, max.overlaps = Inf) +
  theme_void(base_size = 8) + theme(plot.background = element_rect(fill = "white", colour = NA), legend.position = "right") +
  labs(title = "T2D-specific protein network (OR rule, 170 edges among 144 proteins)")
save_f(p2, "Fig2_T2D_network", 10, 7.5)

# ---- Figure 3 (main): hubs vs differential expression, and TIMP1 --------------------
an <- read.csv(file.path(OUT_DIR, "A2_node_table.csv")); an$cls <- ifelse(an$hub_class == "", "Other", an$hub_class)
p3a <- ggplot(an, aes(abs(logOR_perSD), degree, colour = factor(cls, levels = names(cls_cols)))) +
  geom_point(alpha = 0.85, size = 1.5) +
  geom_text_repel(data = an[an$cls != "Other" | an$rank_abs_effect <= 8 | an$degree >= 7, ], aes(label = protein),
                  size = 2.3, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
  scale_colour_manual(values = cls_cols, name = NULL) +
  labs(x = "|log OR per SD| for incident T2D (DE, selection set)", y = "T2D-specific degree",
       title = "a  Rewiring vs differential expression") +
  theme(legend.position = c(0.8, 0.8), legend.background = element_blank(), plot.title = element_text(size = 8, face = "bold"))
tn <- read.csv(file.path(OUT_DIR, "A3_TIMP1_neighbourhood.csv"))
tn$rep <- tn$concordant & tn$p_heldout < 0.05
gt <- graph_from_data_frame(data.frame(from = "TIMP1", to = tn$partner, w = abs(tn$weight), change = tn$in_cases, string = tn$STRING_w > 0, rep = tn$rep),
                            vertices = data.frame(name = c("TIMP1", tn$partner), cls = c("Confirmed hub", tn$partner_class)), directed = FALSE)
p3b <- ggraph(gt, layout = "star", center = 1) +
  geom_edge_link(aes(edge_colour = change, edge_width = w, edge_linetype = string)) +
  scale_edge_colour_manual(values = chg_cols, name = "Dependence\nin cases") +
  scale_edge_width(range = c(0.3, 2.2), name = "|T2D coef.|") +
  scale_edge_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), labels = c("no", "yes"), name = "STRING") +
  geom_node_point(aes(colour = cls), size = 3.5) + scale_colour_manual(values = cls_cols, name = NULL) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2.4) +
  theme_void(base_size = 8) + labs(title = "b  TIMP1 neighbourhood (FWER p = 0.03)") +
  theme(plot.title = element_text(size = 8, face = "bold"), plot.background = element_rect(fill = "white", colour = NA))
save_f(plot_grid(p3a, p3b, ncol = 2, rel_widths = c(1, 1.3)), "Fig3_hubs_DE_TIMP1", 12, 5.2)

# ---- Figure 4 (main): prediction -----------------------------------------------------
pr <- read.csv(file.path(OUT_DIR, "B_model_performance.csv"))
lab <- c(M0 = "Age + sex", M1 = "AND proteins (7)", M2 = "AND proteins + 5 products", M3a = "TIMP1", M3b = "Candidate hubs (12)",
         M3c = "Rewired proteins (21)", M3d = "Hubs + TIMP1 products", M3e = "Top 10% central (15)", M4 = "DE, q < 0.05 (303)",
         M4b_k7 = "Top-7 DE", M4b_k12 = "Top-12 DE", M8 = "Lasso (173)", M8b_k7 = "Lasso path, first 7", M8b_k12 = "Lasso path, first 12",
         M4b_k15 = "Top-15 DE", M4b_k21 = "Top-21 DE", M8b_k15 = "Lasso path, first 15", M8b_k21 = "Lasso path, first 21",
         M5 = "QDiabetes variables", M6 = "QDiabetes + HbA1c", M7 = "QDiabetes + HbA1c + PRS")
pr$label <- ifelse(pr$model %in% names(lab), lab[pr$model],
                   paste(ifelse(pr$base == "M5", "QDiabetes", "QDiabetes + HbA1c"), "+", lab[sub("^M[56]\\+", "", pr$model)]))
pr$panel <- ifelse(pr$base == "M0", "a  Protein sets (with age and sex)", ifelse(pr$model %in% c("M5", "M6", "M7"), "b  Clinical benchmark and added proteins", "b  Clinical benchmark and added proteins"))
pr$type <- ifelse(grepl("^M[123]|\\+M[23]", pr$model), "Network-based", ifelse(pr$model %in% c("M0", "M5", "M6", "M7"), "Reference", "Standard selection (DE / lasso)"))
pr$label <- factor(pr$label, levels = rev(unique(pr$label)))
p4 <- ggplot(pr, aes(AUC_test, label, colour = type)) +
  geom_pointrange(aes(xmin = AUC_test_lo, xmax = AUC_test_hi), size = 0.25) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(`Network-based` = col_obs, `Standard selection (DE / lasso)` = "#7570b3", Reference = "grey35"), name = NULL) +
  labs(x = "Test AUC, held-out 30% (95% CI, DeLong)", y = NULL, title = "10-year incident T2D: discrimination") +
  theme(legend.position = "top")
save_f(p4, "Fig4_prediction", 10, 7)

# ---- Supplementary ---------------------------------------------------------------------
# S1: permutation null, 10-year vs 5-year, AND and OR (copy of A1 in the figure style)
file.copy(file.path(OUT_DIR, c("FigS_A1_permutation_null.pdf", "FigS_A1_permutation_null.png")),
          file.path(FIG_DIR, c("FigS1_permutation_null_5y_10y.pdf", "FigS1_permutation_null_5y_10y.png")), overwrite = TRUE)
# S2: enrichment of communities and protein sets
en <- read.csv(file.path(OUT_DIR, "A2b_enrichment_ORA.csv"))
en <- en[en$q_BH < 0.05, ]
en <- do.call(rbind, lapply(split(en, en$query), function(x) { x <- x[order(x$q_BH, -x$overlap), ]
  unique(rbind(head(x, 4), x[grepl("(^|,)TIMP1(,|$)", x$genes), ])) }))   # top 4 per query + terms containing TIMP1
en$term <- sub("^(KEGG|REACTOME|GOBP)_", "", en$set); en$term <- substr(gsub("_", " ", en$term), 1, 55)
en$term <- paste0(en$term, " [", en$db, "]")
pS2 <- ggplot(en, aes(overlap / query_size, reorder(term, -log10(q_BH)), size = overlap, colour = -log10(q_BH))) +
  geom_point() + facet_grid(query ~ ., scales = "free_y", space = "free_y") +
  scale_colour_viridis_c(name = "-log10 q") + scale_size_continuous(range = c(1.5, 4), name = "Overlap") +
  labs(x = "Share of the query proteins in the term", y = NULL, title = "Over-representation (panel background, BH within query x database)") +
  theme(strip.text.y = element_text(angle = 0, size = 6.5), axis.text.y = element_text(size = 6))
save_f(pS2, "FigS2_enrichment", 9, 11)
# S3: hub test: observed degree vs protein-specific null 95th percentile
hb <- read.csv(file.path(PERM_DIR, "T2D_hub_degree_test.csv")); hb <- hb[hb$obs >= 3, ]
hb$cls <- ifelse(hb$protein %in% confirmed_hubs, "Confirmed hub", ifelse(hb$protein %in% candidate_hubs, "Candidate hub", ifelse(hb$protein %in% rewired, "Rewired", "Other")))
pS3 <- ggplot(hb, aes(reorder(protein, obs))) +
  geom_col(aes(y = obs, fill = factor(cls, levels = names(cls_cols)))) +
  geom_point(aes(y = null_q95), shape = 4, size = 1.6) + geom_point(aes(y = null_mean), shape = 1, size = 1.4) +
  scale_fill_manual(values = cls_cols, name = NULL) + coord_flip() +
  labs(x = NULL, y = "T2D degree (bar); null mean (o) and 95th percentile (x)",
       title = "Hub test: proteins with T2D degree >= 3", subtitle = "Candidate hubs: p <= 0.01 and degree >= 5 (FDR 5.7%); rewired: p <= 0.01 and degree >= 3 (FDR 4.9%)") +
  theme(legend.position = "top")
save_f(pS3, "FigS3_hub_test", 6, 7)
# S4: external references vs the permutation null (STRING fraction, physical, DisGeNET edges)
ph <- read.csv(file.path(OUT_DIR, "A4d_STRING_physical_by_permutation.csv"))
dc <- read.csv(file.path(OUT_DIR, "A4c_disease_context_edges.csv"))
pS4a <- null_panel(sv$OR_frac[sv$b > 0], sv$OR_frac[sv$b == 0], "Share of T2D edges with STRING support (prior W > 0)", "a  STRING (combined score, v12.0)")
pS4b <- null_panel(ph$any[ph$b > 0], ph$any[ph$b == 0], "Share of T2D edges with a physical interaction", "b  STRING physical (v12.5, score >= 150)")
pS4c <- null_panel(rp$p05[-1], rp$p05[1], "Edges with held-out p < 0.05 (any sign)", "c  Held-out replication: p < 0.05")
save_f(plot_grid(pS4a, pS4b, pS4c, ncol = 3), "FigS4_external_references", 9, 2.9)
# S5: prediction with the first design (age- and sex-matched training set), for the record
pv <- read.csv(file.path(OUT_DIR, "B_model_performance_v1_matched.csv"))
pv$label <- ifelse(pv$model %in% names(lab), lab[pv$model], paste(ifelse(pv$base == "M5", "QDiabetes", "QDiabetes + HbA1c"), "+", lab[sub("^M[56]\\+", "", pv$model)]))
pv$label <- factor(pv$label, levels = rev(unique(pv$label)))
pS5 <- ggplot(pv, aes(AUC_test, label)) + geom_pointrange(aes(xmin = AUC_test_lo, xmax = AUC_test_hi), size = 0.2) +
  labs(x = "Test AUC (95% CI)", y = NULL, title = "First design: fit on the age- and sex-matched training set (214 cases)",
       subtitle = "Superseded: age and sex carry no information in training (M0 AUC = 0.5)")
save_f(pS5, "FigS5_prediction_matched_design", 6, 6)
