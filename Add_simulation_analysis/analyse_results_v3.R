#################################################
## Summary of the 100-replication run (results_v3)
#################################################
# Usage (from the Add_simulation_analysis folder):
#   Rscript analyse_results_v3.R [results_v3]
# Reads results_v3/<setting>/rep*/evaluation.rds and times.rds and writes to
# results_v3/summary/:
#   metrics_by_replication.csv  one row per replication, rule, symmetrisation
#                               and component type (covariate components of
#                               the same type averaged within a replication)
#   metrics_summary.csv         mean and SD over the replications
#   table_<metric>_<symm>.csv   "mean (SD)" tables, rules x component type x n
#   Fig_<metric>.pdf/.png       boxplots over the replications
#   covariate_effects_by_type.csv  effect metrics, binary vs continuous covariates
#   computing_time.csv          seconds per node and minutes per setting
#   overall_<metric>_<symm>.csv overall F1, FPR, Pearson: all components of a
#                               fit pooled (see below)
# All networks use the OLS refit on the selected predictors. Component types:
#   baseline  network at the mean of the continuous covariates and the
#             reference level of the binary one
#   binary / continuous effect  covariate with an effect on the network
#   null      covariate without an effect (only the density / FPR matter)
# Passes: evaluation.rds (CV rules, BIC/EBIC, stability at the CV
# hyperparameters, post-hoc thresholds), evaluation_stab_bic.rds (stability at
# the BIC hyperparameters), evaluation_relaxed.rds (BIC/EBIC on the relaxed
# RSS, stability at the relaxed-BIC hyperparameters). The bic/ebic rows of
# the later passes duplicate the first pass (checked) and are not used.

suppressPackageStartupMessages({
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) > 0) args[1] else "results_v3"
out_dir <- file.path(res_dir, "summary")
dir.create(out_dir, showWarnings = FALSE)

eval_files <- Sys.glob(file.path(res_dir, "n*_p*", "rep*", "evaluation.rds"))
time_files <- Sys.glob(file.path(res_dir, "n*_p*", "rep*", "times.rds"))
ev <- do.call(rbind, lapply(eval_files, readRDS))
# Second pass (--mode stab_bic): stability selection at the BIC
# hyperparameters, no cross-validation. Its bic/ebic rows duplicate the first
# pass (checked to be identical) and are not used.
sb_files <- Sys.glob(file.path(res_dir, "n*_p*", "rep*", "evaluation_stab_bic.rds"))
if (length(sb_files) > 0) {
  sb <- do.call(rbind, lapply(sb_files, readRDS))
  chk <- merge(ev[ev$rule == "bic", ], sb[sb$rule == "bic", ],
               by = c("n", "rep", "post", "symm", "component"))
  cat("bic rows of the second pass identical to the first:",
      isTRUE(all.equal(chk$F1.x, chk$F1.y)), "(", nrow(chk), "rows )\n")
  ev <- rbind(ev, sb[sb$rule == "stab_bic", ])
}
rx_files <- Sys.glob(file.path(res_dir, "n*_p*", "rep*", "evaluation_relaxed.rds"))
if (length(rx_files) > 0) {
  rx <- do.call(rbind, lapply(rx_files, readRDS))
  chk <- merge(ev[ev$rule %in% c("bic", "ebic"), ], rx[rx$rule %in% c("bic", "ebic"), ],
               by = c("n", "rep", "rule", "post", "symm", "component"))
  cat("bic/ebic rows of the relaxed pass identical to the first:",
      isTRUE(all.equal(chk$F1.x, chk$F1.y)), "(", nrow(chk), "rows )\n")
  ev <- rbind(ev, rx[rx$rule %in% c("bic_relaxed", "ebic_relaxed", "stab_bic_relaxed"), ])
}
cat("Evaluated replications per n:\n")
print(table(unique(ev[, c("n", "rep")])$n))

ev$method <- ifelse(ev$post == "none", ev$rule, paste0(ev$rule, "+", ev$post))
method_levels <- c("min", "1se", "bic", "bic_relaxed", "ebic", "ebic_relaxed", "stab",
                   "stab_bic", "stab_bic_relaxed", "min+pcor", "min+backbone")
method_labels <- c("CV lambda.min", "CV lambda.1se", "BIC (penalised RSS)", "BIC relaxed",
                   "EBIC (penalised RSS)", "EBIC relaxed",
                   "Stability selection (CV hyperparameters)",
                   "Stability selection (BIC hyperparameters)",
                   "Stability selection (BIC relaxed hyperparameters)",
                   "lambda.min + pcor threshold", "lambda.min + backbone")
keep <- method_levels %in% ev$method
ev$method <- factor(ev$method, levels = method_levels[keep], labels = method_labels[keep])
ev$symm <- factor(ev$symm, levels = c("OR", "AND"))
# Effect components split by covariate type (x3 is the binary covariate)
binary_comp <- ev$component %in% c("x3", "x3_1")
ev$component_type <- ifelse(ev$component_type == "effect",
                            ifelse(binary_comp, "binary_effect", "continuous_effect"),
                            ev$component_type)
ev$component_type <- factor(ev$component_type,
                            levels = c("baseline", "binary_effect", "continuous_effect", "null"),
                            labels = c("Baseline", "Binary effect", "Continuous effect",
                                       "Null covariate"))
metrics <- c("density", "TPR", "FPR", "Precision", "F1", "Spearman_TP", "Pearson_TP", "Frobenius")

# Average the covariate components of the same type within a replication
by_rep <- aggregate(ev[, metrics],
                    by = ev[, c("n", "rep", "method", "symm", "component_type")],
                    FUN = mean, na.rm = TRUE)
for (m in metrics) by_rep[[m]][is.nan(by_rep[[m]])] <- NA
write.csv(by_rep, file.path(out_dir, "metrics_by_replication.csv"), row.names = FALSE)

keys <- by_rep[, c("n", "method", "symm", "component_type")]
summ <- aggregate(by_rep$F1, keys, length)
names(summ)[names(summ) == "x"] <- "n_rep"
for (m in metrics) {
  mu <- aggregate(by_rep[[m]], keys, mean, na.rm = TRUE)
  s <- aggregate(by_rep[[m]], keys, sd, na.rm = TRUE)
  stopifnot(identical(mu[, names(keys)], summ[, names(keys)]))
  summ[[paste0(m, "_mean")]] <- mu$x
  summ[[paste0(m, "_sd")]] <- s$x
}
summ <- summ[order(summ$component_type, summ$symm, summ$method, summ$n), ]
write.csv(summ, file.path(out_dir, "metrics_summary.csv"), row.names = FALSE)

# "mean (SD)" tables: one row per component type and method, one column per n
fmt <- function(mu, s) ifelse(is.na(mu), "-", sprintf("%.2f (%.2f)", mu, s))
for (m in metrics) for (sy in c("OR", "AND")) {
  d <- summ[summ$symm == sy, ]
  d$value <- fmt(d[[paste0(m, "_mean")]], d[[paste0(m, "_sd")]])
  w <- reshape(d[, c("component_type", "method", "n", "value")],
               idvar = c("component_type", "method"), timevar = "n", direction = "wide")
  names(w) <- sub("^value\\.", "n = ", names(w))
  w <- w[order(w$component_type, w$method), ]
  write.csv(w, file.path(out_dir, sprintf("table_%s_%s.csv", m, sy)), row.names = FALSE)
}

# Covariate effects split by covariate type: the binary covariate (SD 0.5)
# carries a stronger signal per unit of standardised effect than the uniform
# ones (SD 0.29), which explains most of the spread of the effect metrics
eff <- ev[ev$component_type %in% c("Binary effect", "Continuous effect"), ]
eff$covariate_type <- ifelse(eff$component %in% c("x3", "x3_1"), "binary", "continuous")
groups <- split(eff, eff[, c("n", "method", "symm", "covariate_type")], drop = TRUE)
by_type <- do.call(rbind, lapply(groups, function(d) {
  out <- d[1, c("n", "method", "symm", "covariate_type")]
  out$n_components <- nrow(d)
  for (m in c("TPR", "Precision", "F1", "Pearson_TP", "Spearman_TP")) {
    out[[m]] <- fmt(mean(d[[m]], na.rm = TRUE), sd(d[[m]], na.rm = TRUE))
  }
  out
}))
by_type <- by_type[order(by_type$symm, by_type$covariate_type, by_type$method, by_type$n), ]
write.csv(by_type, file.path(out_dir, "covariate_effects_by_type.csv"), row.names = FALSE)

# Boxplots over the replications
plot_metric <- function(m, label, types) {
  d <- by_rep[by_rep$component_type %in% types & !is.na(by_rep[[m]]), ]
  d$n <- factor(d$n, levels = sort(unique(d$n)),
                labels = paste0("n = ", format(sort(unique(d$n)), big.mark = ",", trim = TRUE)))
  g <- ggplot(d, aes(x = method, y = .data[[m]], fill = n)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
    facet_grid(component_type ~ symm, scales = "free_y") +
    scale_fill_brewer(palette = "Blues") +
    labs(x = NULL, y = label, fill = NULL) +
    theme_bw(base_size = 10) +
    scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")
  for (ext in c("pdf", "png")) {
    ggsave(file.path(out_dir, sprintf("Fig_%s.%s", m, ext)), g,
           width = 10, height = 2.6 * length(types) + 1.2, dpi = 300)
  }
}
eff_types <- c("Baseline", "Binary effect", "Continuous effect")
all_types <- c(eff_types, "Null covariate")
plot_metric("F1", "F1 score", eff_types)
plot_metric("TPR", "True positive rate", eff_types)
plot_metric("FPR", "False positive rate", all_types)
plot_metric("density", "Density", all_types)
plot_metric("Pearson_TP", "Pearson correlation (true positives)", eff_types)
plot_metric("Spearman_TP", "Spearman correlation (true positives)", eff_types)

# Overall metrics per fit: all components (baseline, covariate effects, null)
# pooled. F1 and FPR from the summed TP / FP / FN / TN over the components;
# Pearson on the true positives as the TP-weighted mean of the per-component
# correlations (the edge values are not stored).
pairs <- ev$p * (ev$p - 1) / 2
ev$TP <- round(ifelse(is.na(ev$TPR), 0, ev$TPR) * ev$n_true_edges)
ev$FP <- ev$n_est_edges - ev$TP
ev$FN <- ev$n_true_edges - ev$TP
ev$TN <- pairs - ev$n_true_edges - ev$FP
ev$wP <- ifelse(is.na(ev$Pearson_TP), 0, ev$TP)
ev$P_TP <- ifelse(is.na(ev$Pearson_TP), 0, ev$Pearson_TP * ev$TP)
ok <- c("n", "rep", "method", "symm")
ov <- aggregate(ev[, c("TP", "FP", "FN", "TN", "wP", "P_TP")], ev[, ok], sum)
ov$F1 <- 2 * ov$TP / (2 * ov$TP + ov$FP + ov$FN)
ov$FPR <- ov$FP / (ov$FP + ov$TN)
ov$Pearson <- ifelse(ov$wP > 0, ov$P_TP / ov$wP, NA)
write.csv(ov, file.path(out_dir, "overall_by_replication.csv"), row.names = FALSE)
for (m in c("F1", "FPR", "Pearson")) for (sy in c("OR", "AND")) {
  d <- ov[ov$symm == sy, ]
  g <- split(d, d[, c("method", "n")], drop = TRUE)
  s <- do.call(rbind, lapply(g, function(x) {
    v <- if (m == "FPR") 1000 * x[[m]] else x[[m]]
    data.frame(method = x$method[1], n = x$n[1],
               value = sprintf(if (m == "Pearson") "%.3f (%.3f)" else if (m == "FPR") "%.1f (%.1f)" else "%.2f (%.2f)",
                               mean(v, na.rm = TRUE), sd(v, na.rm = TRUE)))
  }))
  w <- reshape(s, idvar = "method", timevar = "n", direction = "wide")
  names(w) <- sub("^value\\.", "n = ", names(w))
  w <- w[order(w$method), c("method", paste0("n = ", sort(unique(s$n))))]
  write.csv(w, file.path(out_dir, sprintf("overall_%s_%s.csv", m, sy)), row.names = FALSE)
  cat("\nOverall", m, if (m == "FPR") "(x 1000)" else "", sy, "\n"); print(w, row.names = FALSE)
}

# Computing time. Per rule, the fits each one needs: BIC/EBIC the lambda
# paths only; CV rules the paths and the cross-validation; stability
# selection additionally the subsamples (at the lambda.min hyperparameters).
tm <- do.call(rbind, lapply(time_files, function(f) {
  t <- readRDS(f)
  data.frame(n = t$n[1], rep = t$rep[1], path = mean(t$path), cv = mean(t$cv),
             stab = mean(t$stab), refit = mean(t$refit), wall_min = t$wall_total[1] / 60)
}))
tm$bic_ebic <- tm$path
tm$cv_rules <- tm$path + tm$cv
tm$stability <- tm$path + tm$cv + tm$stab
sb_time_files <- Sys.glob(file.path(res_dir, "n*_p*", "rep*", "times_stab_bic.rds"))
cols <- c("bic_ebic", "cv_rules", "stability", "refit", "wall_min")
if (length(sb_time_files) > 0) {
  tsb <- do.call(rbind, lapply(sb_time_files, function(f) {
    t <- readRDS(f)
    data.frame(n = t$n[1], rep = t$rep[1], stability_bic = mean(t$path + t$stab))
  }))
  tm <- merge(tm, tsb, by = c("n", "rep"), all.x = TRUE)
  cols <- c("bic_ebic", "cv_rules", "stability", "stability_bic", "refit", "wall_min")
}
rx_time_files <- Sys.glob(file.path(res_dir, "n*_p*", "rep*", "times_relaxed.rds"))
if (length(rx_time_files) > 0) {
  trx <- do.call(rbind, lapply(rx_time_files, function(f) {
    t <- readRDS(f)
    data.frame(n = t$n[1], rep = t$rep[1], relaxed = mean(t$path + t$relaxed),
               stability_bic_relaxed = mean(t$path + t$relaxed + t$stab))
  }))
  tm <- merge(tm, trx, by = c("n", "rep"), all.x = TRUE)
  cols <- append(cols, c("relaxed", "stability_bic_relaxed"), after = which(cols == "bic_ebic"))
}
time_summary <- do.call(rbind, lapply(split(tm, tm$n), function(d) {
  data.frame(n = d$n[1], n_rep = nrow(d),
             t(sapply(cols, function(cl) sprintf("%.1f (%.1f)", mean(d[[cl]]), sd(d[[cl]])))))
}))
rownames(time_summary) <- NULL
time_labels <- c(bic_ebic = "BIC/EBIC (s/node)", relaxed = "BIC/EBIC relaxed (s/node)",
                 stability_bic_relaxed = "Stability, BIC relaxed hyperparameters (s/node)",
                 cv_rules = "CV rules (s/node)",
                 stability = "Stability, CV hyperparameters (s/node)",
                 stability_bic = "Stability, BIC hyperparameters (s/node)",
                 refit = "OLS refit (s/node)", wall_min = "Setting wall time (min, 16 cores)")
names(time_summary) <- c("n", "n_rep", time_labels[cols])
write.csv(time_summary, file.path(out_dir, "computing_time.csv"), row.names = FALSE)
print(time_summary)

cat("Summary written to", out_dir, "\n")
