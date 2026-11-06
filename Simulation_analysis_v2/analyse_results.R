#################################################
## Summary of the full simulation study
#################################################
# Usage (from the Simulation_analysis_v2 folder):
#   Rscript analyse_results.R [results] [p values separated by "," or ":", e.g. 50:100]
# Reads results/n*_p*_q*/rep*/{evaluation,nodes,config}.rds and writes to
# results/summary/:
#   completeness.csv            replications found per setting
#   metrics_by_replication.csv  one row per replication, prior, symmetrisation
#                               and component type (components of the same type
#                               averaged within a replication)
#   metrics_summary.csv         mean and SD over the replications
#   table_<metric>_<symm>.csv   "mean (SD)" tables: rows p, q, component type,
#                               prior; columns n
#   prior_gain_<symm>.csv       paired differences Perfect - None and
#                               Noisy - None (mean (SD) over replications)
#   Fig_<metric>_<symm>.pdf/png boxplots over the replications
#   computing_time.csv          seconds per node, wall minutes per task
#   tuning.csv                  selected hyperparameters and support size
# Component types: baseline (T0), binary / continuous covariate with an
# effect, null covariate (x2 at q = 3).

has_gg <- suppressPackageStartupMessages(requireNamespace("ggplot2", quietly = TRUE))
if (has_gg) library(ggplot2) else cat("ggplot2 not available: figures skipped\n")
args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) > 0) args[1] else "results"
p_keep <- if (length(args) > 1) as.numeric(strsplit(args[2], "[,:]")[[1]]) else NULL
out_dir <- file.path(res_dir, "summary")
dir.create(out_dir, showWarnings = FALSE)

set_dirs <- Sys.glob(file.path(res_dir, "n*_p*_q*"))
par_of <- function(d, k) as.numeric(sub(paste0(".*", k, "([0-9]+).*"), "\\1", basename(d)))
if (!is.null(p_keep)) set_dirs <- set_dirs[par_of(set_dirs, "_p") %in% p_keep]

read_all <- function(file) {
  fs <- unlist(lapply(set_dirs, function(d) Sys.glob(file.path(d, "rep*", file))))
  list(files = fs, data = lapply(fs, readRDS))
}
ev_l <- read_all("evaluation.rds")
nd_l <- read_all("nodes.rds")
cf_l <- read_all("config.rds")
ev <- do.call(rbind, ev_l$data)
nd <- do.call(rbind, nd_l$data)

comp_tab <- aggregate(rep ~ n + p + q, unique(ev[, c("n", "p", "q", "rep")]), length)
names(comp_tab)[4] <- "n_rep"
comp_tab <- comp_tab[order(comp_tab$p, comp_tab$q, comp_tab$n), ]
write.csv(comp_tab, file.path(out_dir, "completeness.csv"), row.names = FALSE)
cat("Replications per setting:\n"); print(comp_tab, row.names = FALSE)
chk <- aggregate(F1 ~ n + p + q + rep, ev, length)
cat("Rows per replication (should be constant within q):\n")
print(table(q = chk$q, rows = chk$F1))

ev$prior <- factor(ev$prior, levels = c("perfect", "noisy", "none"),
                   labels = c("Perfect", "Noisy", "None"))
ev$symm <- factor(ev$symm, levels = c("OR", "AND"))
ev$ctype <- ifelse(ev$type == "baseline", "Baseline",
                   ifelse(!ev$effect, "Null covariate",
                          ifelse(ev$type == "binary", "Binary effect", "Continuous effect")))
ev$ctype <- factor(ev$ctype, levels = c("Baseline", "Binary effect", "Continuous effect",
                                        "Null covariate"))
metrics <- c("density", "TPR", "FPR", "Precision", "F1", "Spearman_TP", "Pearson_TP",
             "Frobenius", "n_est_edges", "n_true_edges")
keys_rep <- c("n", "p", "q", "rep", "prior", "symm", "ctype")
groups <- split(ev, ev[, keys_rep], drop = TRUE)
by_rep <- do.call(rbind, lapply(groups, function(d) {
  out <- d[1, keys_rep]
  for (m in metrics) out[[m]] <- mean(d[[m]], na.rm = TRUE)
  out
}))
for (m in metrics) by_rep[[m]][is.nan(by_rep[[m]])] <- NA
rownames(by_rep) <- NULL
write.csv(by_rep, file.path(out_dir, "metrics_by_replication.csv"), row.names = FALSE)

keys <- c("n", "p", "q", "prior", "symm", "ctype")
summ <- do.call(rbind, lapply(split(by_rep, by_rep[, keys], drop = TRUE), function(d) {
  out <- d[1, keys]; out$n_rep <- nrow(d)
  for (m in metrics) {
    out[[paste0(m, "_mean")]] <- mean(d[[m]], na.rm = TRUE)
    out[[paste0(m, "_sd")]] <- sd(d[[m]], na.rm = TRUE)
  }
  out
}))
summ <- summ[order(summ$symm, summ$p, summ$q, summ$ctype, summ$prior, summ$n), ]
rownames(summ) <- NULL
write.csv(summ, file.path(out_dir, "metrics_summary.csv"), row.names = FALSE)

fmt <- function(mu, s) ifelse(is.na(mu), "-", sprintf("%.2f (%.2f)", mu, s))
wide <- function(d, value) {
  d$value <- value
  w <- reshape(d[, c("p", "q", "ctype", "prior", "n", "value")],
               idvar = c("p", "q", "ctype", "prior"), timevar = "n", direction = "wide")
  names(w) <- sub("^value\\.", "n = ", names(w))
  w[, c("p", "q", "ctype", "prior", paste0("n = ", sort(unique(d$n))))]
}
for (m in c("F1", "TPR", "FPR", "Precision", "Pearson_TP", "Spearman_TP", "density"))
  for (sy in c("OR", "AND")) {
    d <- summ[summ$symm == sy, ]
    w <- wide(d, fmt(d[[paste0(m, "_mean")]], d[[paste0(m, "_sd")]]))
    w <- w[order(w$p, w$q, w$ctype, w$prior), ]
    write.csv(w, file.path(out_dir, sprintf("table_%s_%s.csv", m, sy)), row.names = FALSE)
  }

# Paired gain of the prior: same dataset fitted under the three priors
gain_metrics <- c("F1", "TPR", "FPR", "Precision")
base <- by_rep[by_rep$prior == "None", c("n", "p", "q", "rep", "symm", "ctype", gain_metrics)]
gain <- do.call(rbind, lapply(c("Perfect", "Noisy"), function(pr) {
  a <- by_rep[by_rep$prior == pr, c("n", "p", "q", "rep", "symm", "ctype", gain_metrics)]
  m <- merge(a, base, by = c("n", "p", "q", "rep", "symm", "ctype"), suffixes = c("", ".none"))
  out <- m[, c("n", "p", "q", "rep", "symm", "ctype")]
  out$comparison <- paste(pr, "- None")
  for (g in gain_metrics) out[[g]] <- m[[g]] - m[[paste0(g, ".none")]]
  out
}))
gkeys <- c("n", "p", "q", "symm", "ctype", "comparison")
gsum <- do.call(rbind, lapply(split(gain, gain[, gkeys], drop = TRUE), function(d) {
  out <- d[1, gkeys]
  for (g in gain_metrics) out[[g]] <- fmt(mean(d[[g]], na.rm = TRUE), sd(d[[g]], na.rm = TRUE))
  out$F1_improved <- sprintf("%.0f%%", 100 * mean(d$F1 > 0, na.rm = TRUE))
  out
}))
for (sy in c("OR", "AND")) {
  g <- gsum[gsum$symm == sy, ]
  g <- g[order(g$comparison, g$p, g$q, g$ctype, g$n), ]
  write.csv(g, file.path(out_dir, sprintf("prior_gain_%s.csv", sy)), row.names = FALSE)
}

# Figures: boxplots over replications, x = n, fill = prior, facets ctype x (p, q)
plot_metric <- function(m, label, types, sy) {
  d <- by_rep[by_rep$symm == sy & by_rep$ctype %in% types & !is.na(by_rep[[m]]), ]
  if (nrow(d) == 0) return(invisible())
  d$n <- factor(format(d$n, big.mark = ",", trim = TRUE),
                levels = format(sort(unique(d$n)), big.mark = ",", trim = TRUE))
  pq <- unique(d[order(d$p, d$q), c("p", "q")])
  d$pq <- factor(sprintf("p = %d, q = %d", d$p, d$q),
                 levels = sprintf("p = %d, q = %d", pq$p, pq$q))
  g <- ggplot(d, aes(x = n, y = .data[[m]], fill = prior)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.25) +
    facet_grid(ctype ~ pq, scales = "free_y") +
    scale_fill_manual(values = c(Perfect = "#1b9e77", Noisy = "#d95f02", None = "#7570b3")) +
    labs(x = "Sample size n", y = label, fill = "Prior", title = paste("Symmetrisation:", sy)) +
    theme_bw(base_size = 9) + theme(legend.position = "top",
                                    axis.text.x = element_text(angle = 30, hjust = 1))
  for (ext in c("pdf", "png"))
    ggsave(file.path(out_dir, sprintf("Fig_%s_%s.%s", m, sy, ext)), g,
           width = 2.2 * nlevels(d$pq) + 1.5, height = 2.2 * length(unique(d$ctype)) + 1.2,
           dpi = 300)
}
eff_types <- c("Baseline", "Binary effect", "Continuous effect")
if (has_gg) for (sy in c("OR", "AND")) {
  plot_metric("F1", "F1 score", eff_types, sy)
  plot_metric("TPR", "True positive rate", eff_types, sy)
  plot_metric("FPR", "False positive rate", c(eff_types, "Null covariate"), sy)
  plot_metric("Pearson_TP", "Pearson correlation (true positives)", eff_types, sy)
}

# Computing time and tuning
cf <- do.call(rbind, lapply(seq_along(cf_l$data), function(k) {
  x <- cf_l$data[[k]]; f <- cf_l$files[[k]]
  data.frame(n = par_of(dirname(dirname(f)), "^n"), p = par_of(dirname(dirname(f)), "_p"),
             q = par_of(dirname(dirname(f)), "_q"),
             data_s = x$time_data_s, wall_min = sum(x$wall_s) / 60)
}))
tkeys <- c("n", "p", "q")
node_time <- aggregate(time_s ~ n + p + q, nd, function(v) c(mean = mean(v), max = max(v)))
node_time <- do.call(data.frame, node_time)
wall <- aggregate(cbind(data_s, wall_min) ~ n + p + q, cf, function(v) c(mean = mean(v), max = max(v)))
wall <- do.call(data.frame, wall)
tt <- merge(node_time, wall, by = tkeys)
tt <- tt[order(tt$p, tt$q, tt$n), ]
write.csv(tt, file.path(out_dir, "computing_time.csv"), row.names = FALSE)
cat("\nComputing time (node seconds; task wall minutes over the three priors):\n")
print(format(tt, digits = 3), row.names = FALSE)

tun <- do.call(rbind, lapply(split(nd, nd[, c("n", "p", "q", "prior")], drop = TRUE), function(d) {
  data.frame(n = d$n[1], p = d$p[1], q = d$q[1], prior = d$prior[1],
             lambda_median = median(d$lambda), n_selected_mean = mean(d$n_selected),
             asparse = paste(names(table(d$asparse)), round(100 * prop.table(table(d$asparse))),
                             sep = ":", collapse = " "),
             weight = paste(names(table(d$weight)), round(100 * prop.table(table(d$weight))),
                            sep = ":", collapse = " "))
}))
tun <- tun[order(tun$p, tun$q, tun$prior, tun$n), ]
write.csv(tun, file.path(out_dir, "tuning.csv"), row.names = FALSE)
cat("Summary written to", out_dir, "\n")
