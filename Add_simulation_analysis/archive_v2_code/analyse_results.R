#################################################
## Summaries and figures of the exploratory simulation
#################################################
# Usage (from the Add_simulation_analysis folder, after the array job):
#   Rscript analyse_results.R [results_dir]
# Reads evaluation.rds, cv_diag.rds, times.rds and cv_curves.rds of every
# completed setting and writes tables and figures to <results_dir>/summary/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

res_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(res_dir)) res_dir <- "results"
out_dir <- file.path(res_dir, "summary")
dir.create(out_dir, showWarnings = FALSE)

read_all <- function(name) {
  f <- list.files(res_dir, pattern = paste0("^", name, "\\.rds$"), recursive = TRUE,
                  full.names = TRUE)
  f <- f[!grepl("/summary/", f)]
  lapply(f, function(x) { d <- readRDS(x); attr(d, "file") <- x; d })
}
setting_of <- function(file) {
  m <- regmatches(file, regexec("n(\\d+)_p(\\d+)_([a-z]+)/rep(\\d+)", file))[[1]]
  data.frame(n = as.integer(m[2]), p = as.integer(m[3]), rep = as.integer(m[5]))
}

evaluation <- bind_rows(read_all("evaluation"))
cv_diag <- bind_rows(read_all("cv_diag"))
times <- bind_rows(lapply(read_all("times"), function(d) cbind(d, setting_of(attr(d, "file")))))
cat("Completed settings:", nrow(distinct(evaluation, n, rep)), "\n")

rule_levels <- c("min", "1se", "bic", "ebic", "stab", "min+backbone", "min+pcor")
evaluation <- evaluation %>%
  mutate(method = ifelse(post == "none", rule, paste0(rule, "+", post)),
         method = factor(method, levels = rule_levels),
         n_lab = factor(n, levels = sort(unique(n))))
# The covariate without an effect (all-zero truth) is labelled "null"
evaluation <- evaluation %>%
  mutate(component_type = case_when(
    component == "Baseline" ~ "Baseline",
    component == "Baseline_centred" ~ "Baseline (centred target)",
    n_true_edges == 0 ~ "Covariate (null)",
    TRUE ~ "Covariate (effect)"))

#################################################
## 1. CV diagnostics
#################################################
cv_sum <- cv_diag %>%
  mutate(at_end_min = idx_min == L, at_end_1se = idx_1se == L,
         pos_min = log(lambda_min / lambda_max) / log(lambda_end / lambda_max),
         pos_1se = log(lambda_1se / lambda_max) / log(lambda_end / lambda_max)) %>%
  group_by(n) %>%
  summarise(share_min_at_path_end = mean(at_end_min),
            share_1se_at_path_end = mean(at_end_1se),
            median_pos_min = median(pos_min), median_pos_1se = median(pos_1se),
            median_df_min = median(df_min), median_df_1se = median(df_1se),
            median_cv_gain_1se_to_min = median(cv_gain_1se_to_min), .groups = "drop")
write.csv(cv_sum, file.path(out_dir, "cv_diagnostics.csv"), row.names = FALSE)
print(cv_sum, width = Inf)

p_cv <- cv_diag %>%
  mutate(pos_min = log(lambda_min / lambda_max) / log(lambda_end / lambda_max),
         pos_1se = log(lambda_1se / lambda_max) / log(lambda_end / lambda_max)) %>%
  pivot_longer(c(pos_min, pos_1se), names_to = "rule", values_to = "position") %>%
  mutate(rule = sub("pos_", "lambda.", rule)) %>%
  ggplot(aes(factor(n), position, fill = rule)) +
  geom_boxplot(outlier.size = 0.4) +
  labs(x = "n", y = "Position on the lambda path (0 = lambda_max, 1 = end)",
       title = "Where the CV-selected lambda falls") +
  theme_minimal()
ggsave(file.path(out_dir, "Fig_cv_lambda_position.png"), p_cv, width = 7, height = 4.5, dpi = 300)

curves <- read_all("cv_curves")
if (length(curves) > 0) {
  cdf <- bind_rows(lapply(curves, function(cl) {
    s <- setting_of(attr(cl, "file"))
    bind_rows(lapply(cl, function(nd) bind_rows(lapply(seq_along(nd$curves), function(k) {
      cc <- nd$curves[[k]]
      data.frame(s, node = nd$node, combo = k, lambda = cc$lambda,
                 cvm_rel = cc$cvm / min(cc$cvm), df = cc$df)
    }))))
  })) %>% filter(rep == 1, combo == 1)
  p_curve <- ggplot(cdf, aes(log10(lambda), cvm_rel, colour = factor(node))) +
    geom_line() + facet_wrap(~ n, scales = "free_x", labeller = label_both) +
    labs(x = "log10(lambda)", y = "CV error / minimum", colour = "Node",
         title = "Cross-validation curves (replication 1, first hyperparameter combination)") +
    coord_cartesian(ylim = c(1, 1.1)) + theme_minimal()
  ggsave(file.path(out_dir, "Fig_cv_curves.png"), p_curve, width = 9, height = 6, dpi = 300)
}

#################################################
## 2. Selection rules: density and reconstruction
#################################################
metric_sum <- evaluation %>%
  group_by(n, method, coef_type, symm, component_type) %>%
  summarise(across(c(density, TPR, FPR, F1, MagCor, Spearman_TP, Pearson_TP, Frobenius),
                   list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
            .groups = "drop")
write.csv(metric_sum, file.path(out_dir, "selection_rules_metrics.csv"), row.names = FALSE)

plot_metric <- function(metric, components, coef = "pen", file, ylab = metric) {
  d <- evaluation %>% filter(coef_type == coef, component_type %in% components, !is.na(method))
  g <- ggplot(d, aes(n_lab, .data[[metric]], fill = method)) +
    geom_boxplot(outlier.size = 0.3) +
    facet_grid(component_type ~ symm, scales = "free_y") +
    labs(x = "n", y = ylab, fill = "Rule") + theme_minimal() +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, file), g, width = 10, height = 2.8 * length(components) + 1, dpi = 300)
}
comps <- c("Baseline", "Baseline (centred target)", "Covariate (effect)")
plot_metric("density", c(comps, "Covariate (null)"), file = "Fig_density.png", ylab = "Network density")
plot_metric("F1", comps, file = "Fig_F1.png")
plot_metric("FPR", c(comps, "Covariate (null)"), file = "Fig_FPR.png")
plot_metric("TPR", comps, file = "Fig_TPR.png")
plot_metric("MagCor", comps, file = "Fig_MagCor.png", ylab = "Spearman (all pairs, as in the paper)")

#################################################
## 3. Refitting: penalised vs OLS-refitted coefficients
#################################################
refit <- evaluation %>%
  filter(post == "none", component_type %in% comps) %>%
  select(n, rep, method, symm, component, component_type, coef_type, Spearman_TP, Pearson_TP, Frobenius) %>%
  pivot_longer(c(Spearman_TP, Pearson_TP, Frobenius), names_to = "metric") %>%
  pivot_wider(names_from = coef_type, values_from = value) %>%
  mutate(diff_refit_minus_pen = refit - pen)
refit_sum <- refit %>%
  group_by(n, method, symm, component_type, metric) %>%
  summarise(pen = mean(pen, na.rm = TRUE), refit = mean(refit, na.rm = TRUE),
            diff = mean(diff_refit_minus_pen, na.rm = TRUE), .groups = "drop")
write.csv(refit_sum, file.path(out_dir, "refit_vs_penalised.csv"), row.names = FALSE)
g <- refit %>% filter(symm == "OR", metric != "Frobenius") %>%
  ggplot(aes(factor(n), diff_refit_minus_pen, fill = method)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.size = 0.3) +
  facet_grid(component_type ~ metric, scales = "free_y") +
  labs(x = "n", y = "Refit - penalised (true-positive edges)", fill = "Rule",
       title = "Effect of the OLS refit on the edge weights (OR)") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "Fig_refit.png"), g, width = 9, height = 8, dpi = 300)

#################################################
## 4. Computing time per node
#################################################
# min / 1se need the path and the CV fits; bic / ebic only the path;
# stab adds the subsample fits to "min"; refit is reported separately
time_sum <- times %>%
  mutate(min_1se = path + cv, bic_ebic = path, stab = path + cv + stab) %>%
  group_by(n) %>%
  summarise(across(c(min_1se, bic_ebic, stab, refit), list(mean = mean, sd = sd)),
            nodes = n(), .groups = "drop")
write.csv(time_sum, file.path(out_dir, "computing_time_per_node_seconds.csv"), row.names = FALSE)
print(time_sum, width = Inf)
cat("Outputs written to", out_dir, "\n")
