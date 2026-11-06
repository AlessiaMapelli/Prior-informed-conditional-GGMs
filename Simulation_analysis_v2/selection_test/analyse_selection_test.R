#################################################
## Summary of the selection-rule test
#################################################
# Usage (from the Simulation_analysis_v2 folder):
#   Rscript selection_test/analyse_selection_test.R [selection_test/results]
# Writes to <results>/summary/:
#   metrics.csv     mean (SD) of F1, TPR, Precision, FPR per variant, n, prior,
#                   symmetrisation and component type
#   time.csv        seconds per node of each variant (paths, relaxed RSS,
#                   stability), mean over nodes, priors and replications
#   choices.csv     share of nodes at each alpha / theta, per variant and n

args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) > 0) args[1] else "selection_test/results"
out_dir <- file.path(res_dir, "summary"); dir.create(out_dir, showWarnings = FALSE)
rd <- function(f) do.call(rbind, lapply(Sys.glob(file.path(res_dir, "n*", "rep*", f)), readRDS))
ev <- rd("evaluation.rds"); ch <- rd("choices.rds")
tm <- lapply(Sys.glob(file.path(res_dir, "n*", "rep*", "timing.rds")), readRDS)
cat("Replications per n:\n"); print(table(unique(ev[, c("n", "rep")])$n))

variants <- c(bic_current = "BIC (penalised RSS), current grid",
              bic_extended = "BIC (penalised RSS), extended grid",
              bic_relaxed_current = "BIC relaxed, current grid",
              bic_relaxed_extended = "BIC relaxed, extended grid",
              ebic_relaxed_current = "EBIC relaxed, current grid",
              ebic_relaxed_extended = "EBIC relaxed, extended grid",
              stab_bic = "Stability at BIC (penalised) hyperparameters",
              stab_bic_relaxed = "Stability at BIC relaxed hyperparameters")
ev$ctype <- ifelse(ev$type == "baseline", "Baseline", ifelse(!ev$effect, "Null covariate",
                   ifelse(ev$type == "binary", "Binary effect", "Continuous effect")))
fmt <- function(v) sprintf("%.2f (%.2f)", mean(v, na.rm = TRUE), sd(v, na.rm = TRUE))
keys <- c("variant", "n", "prior", "symm", "ctype")
g <- split(ev, ev[, keys], drop = TRUE)
met <- do.call(rbind, lapply(g, function(d) {
  out <- d[1, keys]
  for (m in c("F1", "TPR", "Precision", "FPR")) out[[m]] <- fmt(d[[m]])
  out$F1_mean <- mean(d$F1, na.rm = TRUE)
  out
}))
met$label <- variants[met$variant]
met <- met[order(met$symm, met$ctype, met$prior, match(met$variant, names(variants)), met$n), ]
write.csv(met, file.path(out_dir, "metrics.csv"), row.names = FALSE)

# Compact F1 table: rows variant, columns n; per symmetrisation, component, prior
for (sy in c("OR", "AND")) {
  d <- met[met$symm == sy & met$ctype != "Null covariate", ]
  w <- reshape(d[, c("ctype", "prior", "label", "n", "F1")], idvar = c("ctype", "prior", "label"),
               timevar = "n", direction = "wide")
  names(w) <- sub("^F1\\.", "n = ", names(w))
  write.csv(w, file.path(out_dir, sprintf("F1_%s.csv", sy)), row.names = FALSE)
}

# Time per node of each variant
paths <- do.call(rbind, lapply(tm, `[[`, "paths"))
stab <- do.call(rbind, lapply(tm, `[[`, "stability"))
per_node <- function(sel, relaxed) {
  d <- paths[sel, ]
  agg <- aggregate(cbind(t_path, t_relaxed) ~ n + rep + prior + node, d, sum)
  gram <- aggregate(t_gram ~ n + rep + prior + node, d, max)  # once per node
  agg <- merge(agg, gram, by = c("n", "rep", "prior", "node"))
  agg$t <- agg$t_path + if (relaxed) agg$t_relaxed + agg$t_gram else 0
  agg
}
tt <- list(
  bic_current = per_node(paths$in_current, FALSE),
  bic_relaxed_current = per_node(paths$in_current, TRUE),
  bic_extended = per_node(paths$in_extended, FALSE),
  bic_relaxed_extended = per_node(paths$in_extended, TRUE))
for (s in intersect(c("stab_bic", "stab_bic_relaxed"), unique(stab$variant))) {
  base <- if (s == "stab_bic") tt$bic_current else tt$bic_relaxed_current
  st <- stab[stab$variant == s, c("n", "rep", "prior", "node", "t")]
  names(st)[5] <- "t_stab"
  m <- merge(base, st, by = c("n", "rep", "prior", "node"))
  m$t <- m$t + m$t_stab
  tt[[s]] <- m
}
time_tab <- do.call(rbind, lapply(names(tt), function(v) {
  a <- aggregate(t ~ n, tt[[v]], function(x) c(mean = mean(x), sd = sd(x)))
  data.frame(variant = v, label = variants[v], n = a$n,
             seconds_per_node = sprintf("%.1f (%.1f)", a$t[, "mean"], a$t[, "sd"]))
}))
time_w <- reshape(time_tab, idvar = c("variant", "label"), timevar = "n", direction = "wide")
names(time_w) <- sub("^seconds_per_node\\.", "n = ", names(time_w))
write.csv(time_w, file.path(out_dir, "time.csv"), row.names = FALSE)
cat("\nSeconds per node (mean (SD) over nodes, priors, replications):\n")
print(time_w[, -1], row.names = FALSE)
cat("\nShare of nodes where both BICs choose the same combination (current grid):\n")
same <- merge(ch[ch$variant == "bic_current", ], ch[ch$variant == "bic_relaxed_current", ],
              by = c("n", "rep", "prior", "node"))
same$same <- same$asparse.x == same$asparse.y & same$weight.x == same$weight.y
print(aggregate(same ~ n, same, mean))

# Hyperparameter choices
cc <- do.call(rbind, lapply(split(ch, ch[, c("variant", "n")], drop = TRUE), function(d) {
  data.frame(variant = d$variant[1], n = d$n[1],
             alpha = paste(names(table(d$asparse)), sprintf("%.0f%%", 100 * prop.table(table(d$asparse))), sep = ":", collapse = " "),
             theta = paste(names(table(d$weight)), sprintf("%.0f%%", 100 * prop.table(table(d$weight))), sep = ":", collapse = " "),
             df_mean = mean(d$df))
}))
cc <- cc[order(match(cc$variant, names(variants)), cc$n), ]
write.csv(cc, file.path(out_dir, "choices.csv"), row.names = FALSE)

cat("\nF1 (OR), mean over replications:\n")
d <- met[met$symm == "OR" & met$ctype != "Null covariate", ]
print(xtabs(F1_mean ~ label + n + ctype + prior, d)[, , , "none"], digits = 2)
cat("Summary written to", out_dir, "\n")
