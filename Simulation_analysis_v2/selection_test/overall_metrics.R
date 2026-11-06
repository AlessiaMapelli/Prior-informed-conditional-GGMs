# Overall metrics per fit: all components (baseline, binary, continuous, null)
# pooled. F1 and FPR from the summed TP / FP / FN / TN; Pearson on true
# positives as the TP-weighted mean of the per-component correlations.
# Usage: Rscript selection_test/overall_metrics.R [selection_test/results]
args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) > 0) args[1] else "selection_test/results"
ev <- do.call(rbind, lapply(Sys.glob(file.path(res_dir, "n*", "rep*", "evaluation.rds")), readRDS))
pairs <- ev$p * (ev$p - 1) / 2
ev$TP <- round(ifelse(is.na(ev$TPR), 0, ev$TPR) * ev$n_true_edges)
ev$FP <- ev$n_est_edges - ev$TP
ev$FN <- ev$n_true_edges - ev$TP
ev$TN <- pairs - ev$n_true_edges - ev$FP
ev$wP <- ifelse(is.na(ev$Pearson_TP), 0, ev$TP)
ev$P_TP <- ifelse(is.na(ev$Pearson_TP), 0, ev$Pearson_TP * ev$TP)
k <- c("n", "rep", "prior", "variant", "symm")
a <- aggregate(ev[, c("TP", "FP", "FN", "TN", "wP", "P_TP")], ev[, k], sum)
a$F1 <- 2 * a$TP / (2 * a$TP + a$FP + a$FN)
a$FPR <- a$FP / (a$FP + a$TN)
a$Pearson <- a$P_TP / a$wP
labels <- c(bic_current = "BIC (current)", bic_extended = "BIC, extended grid",
            bic_relaxed_current = "BIC relaxed", bic_relaxed_extended = "BIC relaxed, extended grid",
            ebic_relaxed_current = "EBIC relaxed", ebic_relaxed_extended = "EBIC relaxed, extended grid",
            stab_bic = "Stability (BIC)", stab_bic_relaxed = "Stability (BIC relaxed)")
fmt <- function(v, d) sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(v), sd(v))
out <- do.call(rbind, lapply(split(a, a[, c("variant", "n", "prior", "symm")], drop = TRUE), function(d)
  data.frame(symm = d$symm[1], prior = d$prior[1], rule = labels[d$variant[1]], n = d$n[1],
             F1 = fmt(d$F1, 2), FPR = fmt(1000 * d$FPR, 1), Pearson = fmt(d$Pearson, 3))))
out$symm <- factor(out$symm, c("OR", "AND")); out$prior <- factor(out$prior, c("perfect", "noisy", "none"))
out$rule <- factor(out$rule, labels)
out <- out[order(out$symm, out$prior, out$rule, out$n), ]
dir.create(file.path(res_dir, "summary"), showWarnings = FALSE)
write.csv(out, file.path(res_dir, "summary", "overall_metrics.csv"), row.names = FALSE)
for (m in c("F1", "FPR", "Pearson")) {
  w <- reshape(out[, c("symm", "prior", "rule", "n", m)], idvar = c("symm", "prior", "rule"),
               timevar = "n", direction = "wide")
  names(w) <- sub(paste0("^", m, "\\."), "n=", names(w))
  write.csv(w, file.path(res_dir, "summary", sprintf("overall_%s.csv", m)), row.names = FALSE)
  cat("\n=====", m, if (m == "FPR") "(x 1000)" else "", "=====\n")
  print(w, row.names = FALSE)
}
