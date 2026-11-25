#################################################
## Reproducibility experiment: summary of the networks estimated from
## independent datasets drawn from the same true network
#################################################
# Usage (from the Simulation_analysis_v2 folder):
#   Rscript --vanilla analyse_repro.R [results_reproducibility]
# For each (n, p, q, rule, prior, symmetrisation, component), over the R datasets:
# - Jaccard similarity |E_r & E_s| / |E_r | E_s| for all pairs r < s of datasets
#   (pairs where both networks are empty are skipped);
# - edge-selection frequency f_e = share of datasets selecting edge e, summarised
#   separately for true edges and non-edges: mean f; share of true edges with
#   f >= 0.8 (stably recovered); number of non-edges with f >= 0.2 and >= 0.5
#   (recurrent false positives); number of edges ever selected; mean 4 f (1 - f)
#   over the edges ever selected (instability index: 0 if every edge is selected
#   always or never, 1 if every edge is selected in half of the datasets);
# - network density: mean, SD and coefficient of variation (CV) of the number of
#   estimated edges, and mean F1 / TPR / FPR from evaluation.rds for reference.
# The null covariate (x2) has no true edges: only density and false-positive
# summaries are meaningful for it.
# Outputs in <res>/summary/: repro_jaccard_pairs.csv, repro_summary.csv,
#   repro_edge_frequency.csv, Fig_repro_jaccard.{png,pdf},
#   Fig_repro_frequency.{png,pdf}, Fig_repro_density.{png,pdf}

suppressPackageStartupMessages({ library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args)) args[1] else "results_reproducibility"
out <- file.path(res_dir, "summary"); dir.create(out, showWarnings = FALSE)

rd <- function(f) {
  fl <- Sys.glob(file.path(res_dir, "n*_p*_q*", "rep*", f))
  do.call(rbind, lapply(fl, readRDS))
}
nets <- rd("networks.rds"); truth <- rd("truth.rds"); ev <- rd("evaluation.rds")
cfg <- unique(nets[, c("n", "p", "q", "rep")])
cat("Datasets per n:\n"); print(table(cfg$n))

rule_lab <- c(bic_relaxed = "BIC relaxed", ebic_relaxed = "EBIC relaxed", stability_bic = "Stability selection")
prior_lab <- c(none = "No prior", noisy = "Noisy prior", perfect = "Perfect prior")
comp_lab <- c(Baseline = "Baseline", x1 = "Continuous effect", x3_1 = "Binary effect", x2 = "Null covariate",
              x1_1 = "Binary effect")

cells <- unique(merge(unique(nets[, c("n", "p", "q", "rule", "prior", "symm")]),
                      unique(ev[, c("n", "p", "q", "component")])))
key_e <- function(d) paste(d$i, d$j)
pairs_all <- list(); summ <- list(); freq_all <- list()
for (k in seq_len(nrow(cells))) {
  cl <- cells[k, ]
  reps <- sort(unique(cfg$rep[cfg$n == cl$n & cfg$p == cl$p & cfg$q == cl$q]))
  R <- length(reps)
  sub <- nets[nets$n == cl$n & nets$p == cl$p & nets$q == cl$q & nets$rule == cl$rule &
                nets$prior == cl$prior & nets$symm == cl$symm & nets$component == cl$component, ]
  E <- lapply(reps, function(r) key_e(sub[sub$rep == r, ]))
  # true edges: identical across datasets (fixed graph), checked
  tr <- truth[truth$n == cl$n & truth$p == cl$p & truth$q == cl$q & truth$component == cl$component, ]
  T_by_rep <- tapply(key_e(tr), tr$rep, function(x) paste(sort(x), collapse = ";"))
  if (length(unique(T_by_rep)) > 1) stop("True support differs across datasets: ", paste(cl, collapse = " "))
  Tset <- unique(key_e(tr[tr$rep == reps[1], ]))
  # Jaccard over pairs
  jac <- c()
  if (R > 1) for (a in 1:(R - 1)) for (b in (a + 1):R) {
    u <- length(union(E[[a]], E[[b]]))
    if (u > 0) jac <- c(jac, length(intersect(E[[a]], E[[b]])) / u)
  }
  if (length(jac)) pairs_all[[k]] <- cbind(cl, jaccard = jac, row.names = NULL)
  # edge-selection frequency
  f <- table(unlist(E)) / R
  ft <- setNames(rep(0, length(Tset)), Tset); ft[intersect(names(f), Tset)] <- f[intersect(names(f), Tset)]
  ff <- f[!names(f) %in% Tset]
  p <- cl$p; n_non <- p * (p - 1) / 2 - length(Tset)
  freq_all[[k]] <- rbind(
    if (length(ft)) data.frame(cl, truth = "True edge", freq = as.numeric(ft), row.names = NULL),
    if (length(ff)) data.frame(cl, truth = "Non-edge", freq = as.numeric(ff), row.names = NULL))
  ne <- sapply(E, length)
  e_ev <- ev[ev$n == cl$n & ev$p == cl$p & ev$q == cl$q & ev$rule == cl$rule & ev$prior == cl$prior &
               ev$symm == cl$symm & ev$component == cl$component, ]
  summ[[k]] <- data.frame(cl, datasets = R, n_true_edges = length(Tset),
    jaccard_mean = if (length(jac)) mean(jac) else NA, jaccard_sd = if (length(jac)) sd(jac) else NA,
    jaccard_median = if (length(jac)) median(jac) else NA,
    true_freq_mean = if (length(ft)) mean(ft) else NA,
    true_share_freq80 = if (length(ft)) mean(ft >= 0.8) else NA,
    true_never = if (length(ft)) sum(ft == 0) else NA,
    false_freq_mean_x1000 = 1000 * sum(ff) / n_non,
    false_n_freq20 = sum(ff >= 0.2), false_n_freq50 = sum(ff >= 0.5),
    ever_selected = length(f), instability = if (length(f)) mean(4 * f * (1 - f)) else NA,
    edges_mean = mean(ne), edges_sd = sd(ne), edges_cv = if (mean(ne) > 0) sd(ne) / mean(ne) else NA,
    density_mean = mean(ne) / (p * (p - 1) / 2),
    F1_mean = mean(e_ev$F1), TPR_mean = mean(e_ev$TPR), FPR_mean_x1000 = 1000 * mean(e_ev$FPR),
    row.names = NULL)
}
pairs_all <- do.call(rbind, pairs_all); summ <- do.call(rbind, summ); freq_all <- do.call(rbind, freq_all)
summ <- summ[order(summ$component, summ$symm, summ$n, summ$rule, summ$prior), ]
write.csv(pairs_all, file.path(out, "repro_jaccard_pairs.csv"), row.names = FALSE)
write.csv(summ, file.path(out, "repro_summary.csv"), row.names = FALSE)
write.csv(freq_all, file.path(out, "repro_edge_frequency.csv"), row.names = FALSE)
cols <- c("n", "rule", "prior", "symm", "jaccard_mean", "true_share_freq80", "false_n_freq20",
          "instability", "edges_mean", "edges_cv", "F1_mean")
for (cn in unique(summ$component)) {
  cat("\n==", cn, "==\n"); print(summ[summ$component == cn, cols], row.names = FALSE, digits = 3)
}

# ---- Figures --------------------------------------------------------------------
lab <- function(d) {
  d$Rule <- factor(rule_lab[d$rule], rule_lab); d$Prior <- factor(prior_lab[d$prior], prior_lab)
  d$Component <- factor(comp_lab[d$component], unique(comp_lab)); d$N <- factor(paste0("n = ", d$n), paste0("n = ", sort(unique(d$n))))
  d
}
pal <- c(`No prior` = "#7f7f7f", `Noisy prior` = "#1f78b4", `Perfect prior` = "#33a02c")
save2 <- function(g, f, w, h) { ggsave(file.path(out, paste0(f, ".png")), g, width = w, height = h, dpi = 300, bg = "white")
                                ggsave(file.path(out, paste0(f, ".pdf")), g, width = w, height = h) }
pj <- lab(pairs_all[pairs_all$component != "x2", ])
g1 <- ggplot(pj, aes(N, jaccard, fill = Prior)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, position = position_dodge(0.8), width = 0.7) +
  facet_grid(Component + symm ~ Rule) + scale_fill_manual(values = pal) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Pairwise Jaccard similarity between datasets", fill = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
save2(g1, "Fig_repro_jaccard", 8, 9)

fq <- lab(freq_all[freq_all$component != "x2", ])
g2 <- ggplot(fq[fq$truth == "True edge", ], aes(freq, colour = Prior, linetype = symm)) +
  stat_ecdf(linewidth = 0.4) + facet_grid(Component + N ~ Rule) + scale_colour_manual(values = pal) +
  labs(x = "Selection frequency across datasets (true edges)", y = "Cumulative share of true edges",
       colour = NULL, linetype = NULL) + theme_bw(base_size = 8) + theme(legend.position = "bottom")
save2(g2, "Fig_repro_frequency", 8, 10)

sd_ <- lab(summ)
g3 <- ggplot(sd_, aes(N, edges_mean, colour = Prior, shape = symm, group = interaction(Prior, symm))) +
  geom_pointrange(aes(ymin = pmax(edges_mean - edges_sd, 0), ymax = edges_mean + edges_sd),
                  position = position_dodge(0.5), size = 0.2) +
  geom_hline(aes(yintercept = n_true_edges), linetype = 2, linewidth = 0.3) +
  facet_grid(Component ~ Rule, scales = "free_y") + scale_colour_manual(values = pal) +
  labs(x = NULL, y = "Number of estimated edges (mean +- SD across datasets); dashed: true", colour = NULL, shape = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
save2(g3, "Fig_repro_density", 8, 8)
cat("\nWritten to", out, "\n")
