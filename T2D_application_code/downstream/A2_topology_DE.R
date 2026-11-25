#################################################
## A2. Topology of the T2D OR graph (degree vs null, communities, centrality,
## hub classes) and differential expression (DE) on the selection set.
## DisGeNET and community enrichment are in A2b (external resources).
## Run from the private project folder: Rscript --vanilla ../T2D_application_code/downstream/A2_topology_DE.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))
suppressPackageStartupMessages(library(igraph))

# ---- Degree distribution: observed vs null ----------------------------------
dn <- readRDS(file.path(PERM_DIR, "T2D_degree_obs_null.rds"))
stopifnot(identical(unname(dn$obs), unname(rowSums(A_or))))
kmax <- max(dn$obs, dn$null)
cnt <- function(d) tabulate(d + 1, nbins = kmax + 1)          # proteins with degree 0..kmax
null_cnt <- t(apply(dn$null, 1, cnt))
deg_tab <- data.frame(degree = 0:kmax, observed = cnt(dn$obs),
                      null_mean = colMeans(null_cnt),
                      null_lo = apply(null_cnt, 2, quantile, 0.025),
                      null_hi = apply(null_cnt, 2, quantile, 0.975))
# Tail: proteins with degree >= k
deg_tab$obs_ge  <- rev(cumsum(rev(deg_tab$observed)))
deg_tab$null_ge <- rev(cumsum(rev(deg_tab$null_mean)))
write.csv(deg_tab, file.path(OUT_DIR, "A2_degree_distribution.csv"), row.names = FALSE)
cat("Connected proteins: observed", sum(dn$obs > 0), "| null mean", round(mean(rowSums(dn$null > 0)), 1), "\n")
cat("Proteins with degree >= 3 / 5 / 10: observed", sum(dn$obs >= 3), sum(dn$obs >= 5), sum(dn$obs >= 10),
    "| null mean", round(mean(rowSums(dn$null >= 3)), 2), round(mean(rowSums(dn$null >= 5)), 2),
    round(mean(rowSums(dn$null >= 10)), 2), "\n")

tail_df <- rbind(data.frame(k = deg_tab$degree, n = deg_tab$obs_ge, what = "Observed"),
                 data.frame(k = deg_tab$degree, n = deg_tab$null_ge, what = "Null mean (B = 100)"))
tail_df <- tail_df[tail_df$k >= 1 & tail_df$n > 0, ]
g <- ggplot(tail_df, aes(k, n, colour = what)) + geom_step(linewidth = 0.7) +
  scale_y_log10() + scale_colour_manual(values = c(Observed = "#d95f02", `Null mean (B = 100)` = "grey40"), name = NULL) +
  labs(x = "T2D degree k (OR graph)", y = "Proteins with degree >= k (log scale)",
       title = "T2D-specific degree: observed vs permutation null") +
  theme(legend.position = c(0.75, 0.85))
save_fig(g, "FigS_A2_degree_vs_null", 5, 4)

# ---- Graph, communities, centrality (connected proteins only) ---------------
Wabs <- abs(W_or)
g_all <- graph_from_adjacency_matrix(Wabs, mode = "undirected", weighted = TRUE, diag = FALSE)
g_con <- induced_subgraph(g_all, V(g_all)[degree(g_all) > 0])
comp <- components(g_con)
cat("Connected proteins:", vcount(g_con), "| edges:", ecount(g_con), "| components:", comp$no,
    "| largest component:", max(comp$csize), "\n")

set.seed(123)
cfg <- cluster_fast_greedy(g_con)
cat("Communities (fast greedy):", length(cfg), "| modularity:", round(modularity(cfg), 3),
    "| sizes >= 5:", sum(sizes(cfg) >= 5), "\n")

w <- E(g_con)$weight
node <- data.frame(protein = V(g_con)$name, community = as.integer(membership(cfg)),
                   component = comp$membership,
                   degree = degree(g_con),
                   eigen_centrality = eigen_centrality(g_con)$vector,
                   page_rank = page_rank(g_con, directed = FALSE, weights = w)$vector,
                   betweenness = betweenness(g_con, directed = FALSE, weights = 1 / w),
                   hub_score = hub_score(g_con)$vector, stringsAsFactors = FALSE)
z <- function(x) as.numeric(scale(x))
node$av_centrality <- rowMeans(cbind(z(node$degree), z(node$eigen_centrality),
                                     z(node$page_rank), z(node$betweenness)))
node$top10_central <- node$av_centrality > quantile(node$av_centrality, 0.90)
cat("Top 10% central proteins (av_centrality, connected only):", sum(node$top10_central), "\n")
print(node[order(-node$av_centrality), c("protein", "degree", "betweenness", "av_centrality")][1:15, ],
      row.names = FALSE, digits = 3)

# ---- Hub classes -------------------------------------------------------------
node <- merge(node, setNames(hub[, c("protein", "null_mean", "null_q95", "p_marg", "p_maxT")],
                             c("protein", "null_mean_degree", "null_q95_degree", "p_hub_marginal", "p_hub_FWER")),
              by = "protein", all.x = TRUE)
node$hub_class <- ifelse(node$protein %in% confirmed_hubs, "Confirmed hub",
                  ifelse(node$protein %in% candidate_hubs, "Candidate hub",
                  ifelse(node$protein %in% rewired, "Rewired", "")))

# ---- Differential expression on the selection set ----------------------------
sd_ <- selection_data
y <- as.integer(sd_$Diab == "Incident")
cat("\nDE on the selection set:", nrow(sd_), "participants,", sum(y), "cases\n")
de <- do.call(rbind, lapply(prot, function(j) {
  x <- as.numeric(scale(sd_[[j]]))
  m <- glm(y ~ x + Age + Sex, family = binomial, data = data.frame(y, x, Age = sd_$Age, Sex = sd_$Sex))
  s <- summary(m)$coefficients["x", ]
  data.frame(protein = j, logOR_perSD = s[1], se = s[2], p = s[4],
             p_wilcoxon = wilcox.test(sd_[[j]][y == 1], sd_[[j]][y == 0])$p.value,
             median_diff = median(sd_[[j]][y == 1]) - median(sd_[[j]][y == 0]))
}))
de$q_BH <- p.adjust(de$p, "BH"); de$q_BH_wilcoxon <- p.adjust(de$p_wilcoxon, "BH")
de$OR_perSD <- exp(de$logOR_perSD)
de$rank_abs_effect <- rank(-abs(de$logOR_perSD), ties.method = "first")
de$mean_effect_T2D <- fit$mean_effects[de$protein, "Diab_Incident"]
cat("DE proteins, BH q < 0.05: logistic", sum(de$q_BH < 0.05), "| Wilcoxon", sum(de$q_BH_wilcoxon < 0.05), "\n")
cat("Spearman(logOR, prCGGM mean effect):", round(cor(de$logOR_perSD, de$mean_effect_T2D, method = "spearman"), 2), "\n")
print(head(de[order(de$rank_abs_effect), c("protein", "OR_perSD", "q_BH", "rank_abs_effect")], 21), row.names = FALSE, digits = 3)

# ---- Node table (all 366 proteins) -------------------------------------------
all_nodes <- merge(data.frame(protein = prot), node, by = "protein", all.x = TRUE)
all_nodes$degree[is.na(all_nodes$degree)] <- 0
all_nodes <- merge(all_nodes, de, by = "protein")
all_nodes$hub_class[is.na(all_nodes$hub_class)] <- ""
all_nodes$top10_central[is.na(all_nodes$top10_central)] <- FALSE
write.csv(all_nodes[order(-all_nodes$degree, all_nodes$protein), ],
          file.path(OUT_DIR, "A2_node_table.csv"), row.names = FALSE)
write.csv(de[order(de$rank_abs_effect), ], file.path(OUT_DIR, "A2_DE_selection_set.csv"), row.names = FALSE)

# Network classes vs DE ranks
cmp <- all_nodes[all_nodes$hub_class != "" | all_nodes$top10_central,
                 c("protein", "degree", "hub_class", "top10_central", "OR_perSD", "q_BH", "rank_abs_effect")]
cat("\nHubs / rewired / central proteins and their DE rank (out of 366):\n")
print(cmp[order(-cmp$degree), ], row.names = FALSE, digits = 3)
for (s in list(list("candidate hubs", candidate_hubs), list("rewired", rewired),
               list("top 10% central", all_nodes$protein[all_nodes$top10_central]))) {
  k <- length(s[[2]]); topk <- de$protein[de$rank_abs_effect <= k]
  cat(sprintf("Overlap %s (k = %d) with top-%d DE: %d | median DE rank %.0f\n", s[[1]], k, k,
              length(intersect(s[[2]], topk)), median(de$rank_abs_effect[de$protein %in% s[[2]]])))
}
cat("Spearman(T2D degree, |logOR|) over 366 proteins:",
    round(cor(all_nodes$degree, abs(all_nodes$logOR_perSD), method = "spearman"), 2), "\n")

g <- ggplot(all_nodes, aes(abs(logOR_perSD), degree, colour = factor(hub_class, levels = c("Confirmed hub", "Candidate hub", "Rewired", "")))) +
  geom_point(alpha = 0.8, size = 1.6) +
  ggrepel::geom_text_repel(data = all_nodes[all_nodes$hub_class != "" | all_nodes$rank_abs_effect <= 12 | all_nodes$degree >= 7, ],
            aes(label = protein), size = 2.3, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
  scale_colour_manual(values = c(`Confirmed hub` = "#d95f02", `Candidate hub` = "#7570b3", Rewired = "#1b9e77", "grey70"),
                      labels = c("Confirmed hub", "Candidate hub", "Rewired", "Other"), name = NULL) +
  labs(x = "|log OR per SD| for incident T2D (selection set, adjusted for age and sex)",
       y = "T2D-specific degree (OR graph)", title = "Network rewiring vs differential expression") +
  theme(legend.position = "top")
save_fig(g, "Fig_A2_degree_vs_DE", 7, 5)
