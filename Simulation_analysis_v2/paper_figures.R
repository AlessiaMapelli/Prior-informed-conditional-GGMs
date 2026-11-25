#################################################
## Figures and tables of the simulation study for the paper
#################################################
# Usage (from the Simulation_analysis_v2 folder):
#   Rscript --vanilla paper_figures.R [out_dir] [ratio_max]
# out_dir:   output folder (default paper_figures)
# ratio_max: settings with P / n >= ratio_max are greyed out as outside the
#            operating range of the default rule, where
#            P = (p - 1)(q + 1) is the number of predictors per node
#            (default 0.5; see "Operating range" below)
#
# Inputs (summaries written by the analysis scripts):
#   results/summary/metrics_by_replication.csv       main study, BIC relaxed
#   results/summary/computing_time.csv               main study, time per node
#   results_ebic_relaxed/summary/metrics_by_replication.csv   n = 1000, EBIC relaxed
#   results_stability/summary/metrics_by_replication.csv      stability pilot (optional)
#   ../Add_simulation_analysis/results_v3/summary/overall_by_replication.csv  Part 1
#   ../Add_simulation_analysis/results_v3/summary/metrics_by_replication.csv  Part 1
#   selection_test/results_100/summary/overall_metrics.csv    Part 2
#
# Outputs in out_dir (PDF for the paper, PNG for inspection):
#   Main text
#     Fig2_F1_q3.pdf          F1, q = 3, AND: baseline, binary, continuous effect
#     Fig3_TPR_FPR_q3.pdf     TPR (baseline, binary effect) and FPR x 1000
#                             (baseline, null covariate), q = 3, AND
#   Supplementary (section numbers as in paper_revision_plan.md)
#     FigS_time.pdf           S2: computing time per node
#     FigS_q1.pdf             S4: F1, TPR, FPR x 1000, q = 1, AND
#     FigS_q3_other.pdf       S4: TPR of the continuous effect, FPR x 1000 of the
#                             binary and continuous effects, q = 3, AND
#     FigS_AND_OR.pdf         S4: F1 with AND vs OR, Noisy prior
#     FigS_Pearson.pdf        S4: Pearson correlation on the true positives
#     FigS_rules_overall.pdf  S5: overall F1 and FPR of the selection rules (Part 1)
#     FigS_rules_components.pdf S5: F1 per component of the selection rules (Part 1)
#     FigS_n1000.pdf          S7: BIC vs EBIC relaxed (and stability) at n = 1000
#     FigS_n1000_prior.pdf    S7: effect of the prior (None / Noisy / Perfect) at
#                             n = 1000, per rule and symmetrisation
#     TableS_main_full.csv    S4: all metrics, mean and SD, every setting
#     TableS_main_compact.csv S4: F1 and FPR x 1000, Noisy prior, AND
#     TableS_part2.csv        S6: overall F1 / FPR of BIC and EBIC relaxed across priors
#
# Operating range: the default rule (BIC relaxed) selects saturated supports
# when P approaches n (README, "EBIC relaxed at n = 1000: results"). In the
# study this happens at n = 1000 for P = 996-1996, while P = 498 (ratio 0.498)
# behaves well, so the cells with P / n >= 0.5 are greyed out: n = 1000 with
# p = 250, q = 3 (P = 996); p = 500, q = 1 (P = 998); p = 500, q = 3 (P = 1996).
# F1 and TPR are drawn faded there; FPR is not drawn (it would set the scale)
# and its mean over replications and priors is written in the cell instead.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "paper_figures"
ratio_max <- if (length(args) >= 2) as.numeric(args[2]) else 0.5
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(ggplot2))

prior_cols <- c(Perfect = "#1b9e77", Noisy = "#d95f02", None = "#7570b3")
base_theme <- theme_bw(base_size = 8) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        axis.text.x = element_text(angle = 30, hjust = 1))
save_fig <- function(g, name, width, height) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), g, width = width,
           height = height, units = "in", dpi = 300)
  }
  cat("Written", name, "\n")
}
n_levels <- c(1000, 5000, 10000, 20000)
n_lab <- function(n) factor(format(n, big.mark = ","),
                            levels = format(n_levels, big.mark = ","))

# Main study ---------------------------------------------------------------
main <- read.csv("results/summary/metrics_by_replication.csv")
main$P <- (main$p - 1) * (main$q + 1)
main$out_range <- main$P / main$n >= ratio_max
main$prior <- factor(main$prior, levels = names(prior_cols))
main$n_f <- n_lab(main$n)
main$p_f <- factor(paste("p =", main$p), levels = paste("p =", c(50, 100, 250, 500)))
cat("Out of range (P / n >=", ratio_max, "):\n")
print(unique(main[main$out_range, c("n", "p", "q", "P")]), row.names = FALSE)

# Long format of selected (metric, component) pairs; FPR per 1000
panel_data <- function(d, spec) {
  do.call(rbind, lapply(seq_len(nrow(spec)), function(i) {
    x <- d[d$ctype == spec$ctype[i], ]
    v <- x[[spec$metric[i]]]
    if (spec$metric[i] == "FPR") v <- 1000 * v
    data.frame(x[, c("n_f", "p_f", "prior", "out_range", "rep")], value = v,
               row = spec$label[i], metric = spec$metric[i])
  }))
}

# Grid of boxplots: rows = (metric, component), columns = p, x = n, fill = prior
grid_plot <- function(d, spec, ylab = NULL) {
  pd <- panel_data(d, spec)
  pd$row <- factor(pd$row, levels = spec$label)
  shade <- unique(pd[pd$out_range, c("n_f", "p_f", "row", "metric")])
  shade$x <- as.integer(shade$n_f)
  # FPR in the shaded cells: median written instead of the boxes
  fpr_out <- pd$metric == "FPR" & pd$out_range
  lab <- if (any(fpr_out)) aggregate(value ~ n_f + p_f + row, pd[fpr_out, ], mean) else
    data.frame()
  pd$value[fpr_out] <- NA
  g <- ggplot(pd, aes(x = n_f, y = value, fill = prior)) +
    geom_rect(data = shade, aes(xmin = x - 0.5, xmax = x + 0.5),
              ymin = -Inf, ymax = Inf, fill = "grey88", inherit.aes = FALSE) +
    geom_boxplot(aes(alpha = out_range), outlier.size = 0.3, linewidth = 0.25,
                 position = position_dodge(width = 0.8), na.rm = TRUE) +
    scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.3), guide = "none") +
    scale_fill_manual(values = prior_cols, name = "Prior") +
    facet_grid(row ~ p_f, scales = "free_y", switch = "y") +
    labs(x = "Sample size n", y = ylab) + base_theme +
    theme(strip.placement = "outside")
  if (nrow(lab) > 0) {
    g <- g + geom_text(data = lab, aes(x = n_f, y = Inf,
                                       label = paste0("P/n >= ", ratio_max,
                                                      "\nmean ", round(value))),
                       vjust = 1.3, size = 2, inherit.aes = FALSE)
  }
  g
}
spec <- function(metric, ctype, label) data.frame(metric, ctype, label)

and3 <- main[main$symm == "AND" & main$q == 3, ]
and1 <- main[main$symm == "AND" & main$q == 1, ]

# Figure 2: F1, q = 3
s2 <- spec("F1", c("Baseline", "Binary effect", "Continuous effect"),
           c("F1: baseline", "F1: binary effect", "F1: continuous effect"))
save_fig(grid_plot(and3, s2) + scale_y_continuous(limits = c(0, 1)),
         "Fig2_F1_q3", 7, 5.2)

# Figure 3: TPR and FPR x 1000, q = 3
s3 <- spec(c("TPR", "TPR", "FPR", "FPR"),
           c("Baseline", "Binary effect", "Baseline", "Null covariate"),
           c("TPR: baseline", "TPR: binary effect", "FPR x 1000: baseline",
             "FPR x 1000: null covariate"))
save_fig(grid_plot(and3, s3), "Fig3_TPR_FPR_q3", 7, 6.5)

# S4: q = 1
s4 <- spec(rep(c("F1", "TPR", "FPR"), each = 2), rep(c("Baseline", "Binary effect"), 3),
           c("F1: baseline", "F1: binary effect", "TPR: baseline",
             "TPR: binary effect", "FPR x 1000: baseline",
             "FPR x 1000: binary effect"))
save_fig(grid_plot(and1, s4), "FigS_q1", 7, 8.5)

# S4: remaining components of q = 3
s5 <- spec(c("TPR", "FPR", "FPR"),
           c("Continuous effect", "Binary effect", "Continuous effect"),
           c("TPR: continuous effect", "FPR x 1000: binary effect",
             "FPR x 1000: continuous effect"))
save_fig(grid_plot(and3, s5), "FigS_q3_other", 7, 5.2)

# S4: AND vs OR, Noisy prior
ao <- main[main$prior == "Noisy" & main$ctype != "Null covariate", ]
ao$row <- factor(paste0("q = ", ao$q, ": ", ao$ctype),
                 levels = c("q = 1: Baseline", "q = 1: Binary effect", "q = 3: Baseline",
                            "q = 3: Binary effect", "q = 3: Continuous effect"))
g <- ggplot(ao, aes(x = n_f, y = F1, fill = symm)) +
  geom_boxplot(aes(alpha = out_range), outlier.size = 0.3, linewidth = 0.25) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.3), guide = "none") +
  scale_fill_manual(values = c(AND = "#4d4d4d", OR = "#bdbdbd"),
                    name = "Symmetrisation") +
  facet_grid(row ~ p_f) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Sample size n", y = "F1 (Noisy prior)") + base_theme
save_fig(g, "FigS_AND_OR", 7, 7.5)

# S4: Pearson correlation on the true positives, AND
pe <- main[main$symm == "AND" & main$ctype != "Null covariate" & !main$out_range, ]
pe$row <- factor(paste0("q = ", pe$q, ": ", pe$ctype), levels = levels(ao$row))
g <- ggplot(pe, aes(x = n_f, y = Pearson_TP, fill = prior)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.25, na.rm = TRUE) +
  scale_fill_manual(values = prior_cols, name = "Prior") +
  facet_grid(row ~ p_f, scales = "free_y") +
  labs(x = "Sample size n",
       y = "Pearson correlation with the true weights (true positives)") +
  base_theme
save_fig(g, "FigS_Pearson", 7, 7.5)

# S2: computing time per node
tm <- read.csv("results/summary/computing_time.csv")
tm$q_f <- paste("q =", tm$q)
g <- ggplot(tm, aes(x = n, y = time_s.mean, colour = factor(p), group = factor(p))) +
  geom_line() + geom_point(size = 1) +
  geom_errorbar(aes(ymin = time_s.mean, ymax = time_s.max), width = 0.03, alpha = 0.5) +
  scale_x_log10(breaks = n_levels, labels = format(n_levels, big.mark = ",")) +
  scale_y_log10() + facet_wrap(~ q_f) +
  scale_colour_viridis_d(name = "p", end = 0.85) +
  labs(x = "Sample size n", y = "Seconds per node (mean; bar up to max)") + base_theme
save_fig(g, "FigS_time", 7, 3)

# Tables S4
agg_metrics <- c("density", "TPR", "FPR", "Precision", "F1", "Spearman_TP", "Pearson_TP")
keys <- c("n", "p", "q", "prior", "symm", "ctype")
mean_d <- aggregate(main[, agg_metrics], main[, keys], mean, na.rm = TRUE)
sd_d <- aggregate(main[, agg_metrics], main[, keys], sd, na.rm = TRUE)
names(sd_d)[-(1:6)] <- paste0(agg_metrics, "_sd")
reps <- aggregate(list(replications = main$rep), main[, keys], length)
full <- merge(merge(mean_d, sd_d, by = keys), reps, by = keys)
full$P <- (full$p - 1) * (full$q + 1)
full <- full[order(full$q, full$p, full$n, full$prior, full$symm, full$ctype), ]
write.csv(full, file.path(out_dir, "TableS_main_full.csv"), row.names = FALSE)
cmp <- full[full$prior == "Noisy" & full$symm == "AND", ]
cmp <- data.frame(q = cmp$q, p = cmp$p, n = cmp$n, P = cmp$P, component = cmp$ctype,
                  F1 = sprintf("%.2f (%.2f)", cmp$F1, cmp$F1_sd),
                  FPR_x1000 = sprintf("%.2f", 1000 * cmp$FPR))
write.csv(cmp, file.path(out_dir, "TableS_main_compact.csv"), row.names = FALSE)

# S5: selection rules (Part 1: p = 100, q = 3, Noisy prior) -----------------
p1_dir <- "../Add_simulation_analysis/results_v3/summary"
rule_order <- c("CV lambda.min", "CV lambda.1se", "lambda.min + pcor threshold",
                "lambda.min + backbone", "BIC (penalised RSS)", "EBIC (penalised RSS)",
                "BIC relaxed", "EBIC relaxed",
                "Stability selection (CV hyperparameters)",
                "Stability selection (BIC hyperparameters)",
                "Stability selection (BIC relaxed hyperparameters)")
mean_sd <- function(d, value, by) {
  m <- aggregate(d[[value]], d[by], mean, na.rm = TRUE)
  s <- aggregate(d[[value]], d[by], sd, na.rm = TRUE)
  m$sd <- s$x
  names(m)[names(m) == "x"] <- "mean"
  m
}
rule_plot <- function(d, facet, xlab) {
  d$method <- factor(d$method, levels = rev(rule_order))
  d$n_f <- n_lab(d$n)
  ggplot(d, aes(y = method, x = mean, colour = symm)) +
    geom_pointrange(aes(xmin = mean - sd, xmax = mean + sd), size = 0.15,
                    position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = c(AND = "black", OR = "grey60"),
                        name = "Symmetrisation") +
    facet + labs(y = NULL, x = xlab) + base_theme
}
if (file.exists(file.path(p1_dir, "overall_by_replication.csv"))) {
  ov <- read.csv(file.path(p1_dir, "overall_by_replication.csv"))
  ov$FPR <- 1000 * ov$FPR
  ov_l <- rbind(cbind(mean_sd(ov, "F1", c("n", "method", "symm")), metric = "Overall F1"),
                cbind(mean_sd(ov, "FPR", c("n", "method", "symm")),
                      metric = "Overall FPR x 1000"))
  save_fig(rule_plot(ov_l, facet_grid(. ~ metric + n_f, scales = "free_x"),
                     "Mean +/- SD over 100 replications"),
           "FigS_rules_overall", 7.5, 3.6)
  cm <- read.csv(file.path(p1_dir, "metrics_by_replication.csv"))
  cm <- cm[cm$component_type != "Null covariate", ]
  cm_l <- mean_sd(cm, "F1", c("n", "method", "symm", "component_type"))
  save_fig(rule_plot(cm_l, facet_grid(component_type ~ n_f), "F1, mean +/- SD"),
           "FigS_rules_components", 7.5, 6)
} else cat("Part 1 summaries not found: S5 figures skipped\n")

# S6: Part 2 table ---------------------------------------------------------
p2_file <- "selection_test/results_100/summary/overall_metrics.csv"
if (file.exists(p2_file)) {
  p2 <- read.csv(p2_file)
  p2 <- p2[p2$rule %in% c("BIC relaxed", "EBIC relaxed", "BIC (current)") & !is.na(p2$n), ]
  p2$rule[p2$rule == "BIC (current)"] <- "BIC (penalised RSS)"
  w <- reshape(p2[, c("symm", "prior", "rule", "n", "F1", "FPR")], idvar = c("symm", "prior", "rule"),
               timevar = "n", direction = "wide")
  write.csv(w[order(w$symm, w$prior, w$rule), ], file.path(out_dir, "TableS_part2.csv"),
            row.names = FALSE)
  cat("Written TableS_part2\n")
}

# S7: n = 1000, BIC relaxed vs EBIC relaxed (and stability if available) ---
small <- list(BIC = main[main$n == 1000, ])
if (file.exists("results_ebic_relaxed/summary/metrics_by_replication.csv")) {
  small$EBIC <- read.csv("results_ebic_relaxed/summary/metrics_by_replication.csv")
}
if (file.exists("results_stability/summary/metrics_by_replication.csv")) {
  small$Stability <- read.csv("results_stability/summary/metrics_by_replication.csv")
}
rule_names <- c(BIC = "BIC relaxed", EBIC = "EBIC relaxed", Stability = "Stability selection")
sm <- do.call(rbind, lapply(names(small), function(r) {
  x <- small[[r]][small[[r]]$n == 1000 & small[[r]]$prior == "Noisy",
                  c("n", "p", "q", "rep", "symm", "ctype", "F1")]
  x$rule <- paste(rule_names[r], "+", x$symm)
  x
}))
sm <- sm[sm$ctype != "Null covariate", ]
sm$P <- (sm$p - 1) * (sm$q + 1)
sm_m <- mean_sd(sm, "F1", c("P", "p", "q", "ctype", "rule"))
sm_m$setting <- factor(sprintf("%d\n(p=%d, q=%d)", sm_m$P, sm_m$p, sm_m$q))
sm_m$setting <- factor(sm_m$setting, levels = unique(sm_m$setting[order(sm_m$P)]))
g <- ggplot(sm_m, aes(x = setting, y = mean, colour = rule, group = rule)) +
  geom_vline(xintercept = sum(sort(unique(sm_m$P)) / 1000 < ratio_max) + 0.5,
             linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)), size = 0.2,
                  position = position_dodge(width = 0.6)) +
  facet_wrap(~ ctype, ncol = 1) +
  scale_colour_brewer(palette = "Paired", name = NULL) +
  labs(x = "Predictors per node P (n = 1000; dashed line: P / n = 0.5)",
       y = "F1, Noisy prior (mean +/- SD)") + base_theme
save_fig(g, "FigS_n1000", 7, 6)

# S7b: n = 1000, effect of the prior for each rule and symmetrisation ------
sp <- do.call(rbind, lapply(names(small), function(r) {
  x <- small[[r]][small[[r]]$n == 1000, c("n", "p", "q", "rep", "prior", "symm", "ctype", "F1")]
  x$rule <- rule_names[r]
  x
}))
sp <- sp[sp$ctype != "Null covariate", ]
sp$P <- (sp$p - 1) * (sp$q + 1)
sp_m <- mean_sd(sp, "F1", c("P", "p", "q", "ctype", "rule", "symm", "prior"))
sp_m$prior <- factor(sp_m$prior, levels = names(prior_cols))
sp_m$row <- factor(paste(sp_m$rule, "+", sp_m$symm),
                   levels = as.vector(outer(rule_names, c("AND", "OR"), paste, sep = " + ")))
sp_m$row <- droplevels(sp_m$row)
sp_m$setting <- factor(sprintf("%d\n(p=%d, q=%d)", sp_m$P, sp_m$p, sp_m$q))
sp_m$setting <- factor(sp_m$setting, levels = unique(sp_m$setting[order(sp_m$P)]))
g <- ggplot(sp_m, aes(x = setting, y = mean, colour = prior)) +
  geom_vline(xintercept = sum(sort(unique(sp_m$P)) / 1000 < ratio_max) + 0.5,
             linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)), size = 0.15,
                  position = position_dodge(width = 0.6)) +
  facet_grid(row ~ ctype) +
  scale_colour_manual(values = prior_cols, name = "Prior") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "Predictors per node P (n = 1000; dashed line: P / n = 0.5)",
       y = "F1 (mean +/- SD over replications)") + base_theme +
  theme(strip.text.y = element_text(size = 6.5))
save_fig(g, "FigS_n1000_prior", 9, 9)

cat("Done:", normalizePath(out_dir), "\n")
