#################################################
## Ablation of the selective prior: prior on the baseline only (prCGGM)
## vs prior on all terms (baseline and covariate-specific)
#################################################
# Usage (from the Simulation_analysis_v2 folder):
#   Rscript --vanilla analyse_ablation.R [results_ablation] [results]
# Arms (same datasets, p = 100, q = 3, replications 1-100, bic_relaxed):
#   No prior                      results/            prior "none"
#   Noisy / Perfect, baseline only results/            priors "noisy", "perfect" (prCGGM)
#   Noisy / Perfect, all terms     results_ablation/   same priors, --prior_on_covariates TRUE
# Only replications present in both folders are used (paired comparison).
# Outputs in <ablation>/summary/:
#   ablation_by_replication.csv  metrics per arm, replication, symmetrisation, component
#   ablation_summary.csv         mean and SD per arm, n, symmetrisation, component
#   ablation_paired.csv          all terms minus baseline only (per prior): mean
#                                paired difference and Wilcoxon signed-rank p
#   Fig_ablation.{png,pdf}       F1 (baseline, binary, continuous), TPR of the binary
#                                effect, FPR x 1000 (binary effect, null covariate), AND
#   Fig_ablation_OR.{png,pdf}    the same with OR

suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
abl_dir <- if (length(args) >= 1) args[1] else "results_ablation"
main_dir <- if (length(args) >= 2) args[2] else "results"
out <- file.path(abl_dir, "summary"); dir.create(out, showWarnings = FALSE, recursive = TRUE)

rd <- function(dir) {
  fl <- Sys.glob(file.path(dir, "n*_p100_q3", "rep*", "evaluation.rds"))
  do.call(rbind, lapply(fl, readRDS))
}
abl <- rd(abl_dir); main <- rd(main_dir)
abl$arm <- paste0(abl$prior, "_all")
main$arm <- ifelse(main$prior == "none", "none", paste0(main$prior, "_baseline"))
keep <- unique(abl[, c("n", "rep")])
cat("Ablation replications per n:\n"); print(table(keep$n))
main <- merge(main, keep)
d <- rbind(main[, names(abl)], abl)
arm_lab <- c(none = "No prior", noisy_baseline = "Noisy prior, baseline only (prCGGM)",
             noisy_all = "Noisy prior, all terms", perfect_baseline = "Perfect prior, baseline only (prCGGM)",
             perfect_all = "Perfect prior, all terms")
comp_lab <- c(Baseline = "Baseline", x3_1 = "Binary effect", x1 = "Continuous effect", x2 = "Null covariate")
d$Arm <- factor(arm_lab[d$arm], arm_lab)
d$Component <- factor(comp_lab[d$component], comp_lab)
d$FPR_x1000 <- 1000 * d$FPR
write.csv(d, file.path(out, "ablation_by_replication.csv"), row.names = FALSE)

metrics <- c("F1", "TPR", "Precision", "FPR_x1000", "n_est_edges", "Pearson_TP")
s <- do.call(rbind, lapply(split(d, list(d$arm, d$n, d$symm, d$component), drop = TRUE), function(x)
  data.frame(arm = x$arm[1], n = x$n[1], symm = x$symm[1], component = x$component[1], reps = nrow(x),
             t(unlist(lapply(metrics, function(m) setNames(c(mean(x[[m]], na.rm = TRUE), sd(x[[m]], na.rm = TRUE)),
                                                           paste0(m, c("_mean", "_sd")))))))))
s <- s[order(s$symm, s$component, s$n, match(s$arm, names(arm_lab))), ]
write.csv(s, file.path(out, "ablation_summary.csv"), row.names = FALSE)

# Paired differences: all terms minus baseline only, per prior
pr <- do.call(rbind, lapply(c("noisy", "perfect"), function(p0) {
  a <- d[d$arm == paste0(p0, "_all"), ]; b <- d[d$arm == paste0(p0, "_baseline"), ]
  m <- merge(a, b, by = c("n", "rep", "symm", "component"), suffixes = c("_all", "_base"))
  do.call(rbind, lapply(split(m, list(m$n, m$symm, m$component), drop = TRUE), function(x)
    do.call(rbind, lapply(c("F1", "TPR", "FPR_x1000", "n_est_edges"), function(v) {
      dv <- x[[paste0(v, "_all")]] - x[[paste0(v, "_base")]]
      ok <- !is.na(dv)
      data.frame(prior = p0, n = x$n[1], symm = x$symm[1], component = x$component[1], metric = v,
                 mean_baseline_only = mean(x[[paste0(v, "_base")]], na.rm = TRUE),
                 mean_all_terms = mean(x[[paste0(v, "_all")]], na.rm = TRUE),
                 mean_diff = mean(dv[ok]),
                 p_wilcoxon = if (sum(ok) > 1 && any(dv[ok] != 0)) suppressWarnings(wilcox.test(dv[ok])$p.value) else NA)
    }))))
}))
pr <- pr[order(pr$symm, pr$component, pr$metric, pr$prior, pr$n), ]
write.csv(pr, file.path(out, "ablation_paired.csv"), row.names = FALSE)
show <- pr[pr$symm == "AND" & ((pr$metric %in% c("F1", "TPR") & pr$component != "x2") |
                                 (pr$metric == "FPR_x1000")), ]
print(show, row.names = FALSE, digits = 3)

# Figures
pal <- c(`No prior` = "#7f7f7f", `Noisy prior, baseline only (prCGGM)` = "#1f78b4", `Noisy prior, all terms` = "#a6cee3",
         `Perfect prior, baseline only (prCGGM)` = "#33a02c", `Perfect prior, all terms` = "#b2df8a")
panels <- rbind(data.frame(component = c("Baseline", "x3_1", "x1"), metric = "F1"),
                data.frame(component = "x3_1", metric = "TPR"),
                data.frame(component = c("x3_1", "x1", "x2"), metric = "FPR_x1000"))
for (sy in c("AND", "OR")) {
  long <- do.call(rbind, lapply(seq_len(nrow(panels)), function(k) {
    x <- d[d$symm == sy & d$component == panels$component[k], ]
    data.frame(x[, c("n", "rep", "Arm")], value = x[[panels$metric[k]]],
               panel = paste0(sub("FPR_x1000", "FPR x 1000", panels$metric[k]), ": ", comp_lab[panels$component[k]]))
  }))
  long$panel <- factor(long$panel, unique(long$panel))
  long$N <- factor(paste0("n = ", long$n), paste0("n = ", sort(unique(long$n))))
  g <- ggplot(long, aes(N, value, fill = Arm)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.75, position = position_dodge(0.85)) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) + scale_fill_manual(values = pal) +
    labs(x = NULL, y = NULL, fill = NULL, title = sprintf("Selective prior vs prior on all terms (p = 100, q = 3, %s)", sy)) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
  f <- if (sy == "AND") "Fig_ablation" else "Fig_ablation_OR"
  ggsave(file.path(out, paste0(f, ".png")), g, width = 8, height = 9, dpi = 300, bg = "white")
  ggsave(file.path(out, paste0(f, ".pdf")), g, width = 8, height = 9)
}
cat("\nWritten to", out, "\n")
