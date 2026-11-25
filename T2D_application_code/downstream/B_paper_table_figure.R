#################################################
## Part B for the paper: prediction table and figure (two drafts of the same
## content, to choose the clearer one).
## Rows = protein sets; columns = base model (age + sex; QDiabetes-2018
## variables; QDiabetes + HbA1c). Each cell: test AUC (95% CI, DeLong) on the
## held-out 30% and the likelihood-ratio test against the base model, with
## both models refitted on the held-out data (valid: the protein sets were
## fixed on the selection set).
## Needs B_prediction.R outputs (B_model_performance.csv, B_model_vs_base.csv).
## Run from the private project folder: Rscript --vanilla ../T2D_application_code/downstream/B_paper_table_figure.R
## Outputs in results_bic_relaxed_10y/downstream/ (tables) and figures/:
##   B_table_paper.csv, B_table_paper.tex (main-text draft)
##   B_table_full.csv (all models, supplementary)
##   figures/Fig4b_prediction_by_base.{png,pdf}
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))
suppressPackageStartupMessages({ library(ggplot2) })

perf <- read.csv(file.path(OUT_DIR, "B_model_performance.csv"))
vb <- read.csv(file.path(OUT_DIR, "B_model_vs_base.csv"))
d <- merge(perf, vb[, c("model", "reference", "delta_AUC_test", "p_DeLong", "LRT_df", "LRT_chisq_test", "p_LRT_test")],
           by = "model", all.x = TRUE)

rows <- data.frame(set = c("", "M2", "M3b", "M3e", "M3c", "M4b_k12", "M4b_k15", "M4b_k21", "M8b_k12", "M8b_k15", "M8b_k21", "M7"),
  label = c("Base model alone", "AND proteins (7) + 5 products", "Candidate hubs (12)",
            "Top 10% central (15)", "Rewired proteins (21)", "Top-12 DE", "Top-15 DE", "Top-21 DE",
            "Lasso path, first 12", "Lasso path, first 15", "Lasso path, first 21", "+ polygenic risk score"),
  type = c("Reference", rep("Network-based", 4), rep("Differential expression", 3), rep("Lasso path", 3), "Reference"))
bases <- c(M0 = "Age + sex", M5 = "QDiabetes", M6 = "QDiabetes + HbA1c")
model_id <- function(set, b) {
  if (set == "") return(b)
  if (set == "M7") return(if (b == "M6") "M7" else NA)
  if (b == "M0") set else paste0(b, "+", set)
}
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, sub("e-0?", "e-", sprintf("%.0e", p)), sprintf("%.3f", p)))
tab <- do.call(rbind, lapply(seq_len(nrow(rows)), function(k) do.call(rbind, lapply(names(bases), function(b) {
  m <- model_id(rows$set[k], b)
  x <- if (is.na(m)) NULL else d[d$model == m, ]
  if (is.null(x) || nrow(x) == 0) return(data.frame(row = k, label = rows$label[k], type = rows$type[k], base = b, model = NA,
                                                     n_proteins = NA, AUC = NA, lo = NA, hi = NA, dAUC = NA, p_DeLong = NA,
                                                     LRT_chisq = NA, LRT_df = NA, p_LRT = NA, reference = NA))
  data.frame(row = k, label = rows$label[k], type = rows$type[k], base = b, model = m, n_proteins = x$n_proteins,
             AUC = x$AUC_test, lo = x$AUC_test_lo, hi = x$AUC_test_hi, dAUC = x$delta_AUC_test, p_DeLong = x$p_DeLong,
             LRT_chisq = x$LRT_chisq_test, LRT_df = x$LRT_df, p_LRT = x$p_LRT_test, reference = x$reference)
}))))
# base-alone rows are the references; their tests against the smaller base
# (QDiabetes vs age + sex, + HbA1c vs QDiabetes) stay in B_table_full.csv
tab[tab$row == 1, c("dAUC", "p_DeLong", "LRT_chisq", "LRT_df", "p_LRT", "reference")] <- NA
write.csv(tab, file.path(OUT_DIR, "B_table_paper.csv"), row.names = FALSE)

# Full table (supplementary): every model, with base, AUC, CI, dAUC, DeLong, LRT
full <- d[, c("model", "base", "n_proteins", "n_products", "AUC_train", "AUC_test", "AUC_test_lo", "AUC_test_hi",
              "sens_test", "spec_test", "reference", "delta_AUC_test", "p_DeLong", "LRT_df", "LRT_chisq_test", "p_LRT_test")]
write.csv(full[order(match(full$base, names(bases)), full$model), ], file.path(OUT_DIR, "B_table_full.csv"), row.names = FALSE)

# ---- LaTeX draft (main text) ---------------------------------------------------------
cell <- function(x) if (is.na(x$AUC)) c("--", "") else
  c(sprintf("%.3f (%.3f--%.3f)", x$AUC, x$lo, x$hi), if (is.na(x$p_LRT)) "ref." else fmt_p(x$p_LRT))
lines <- sapply(seq_len(nrow(rows)), function(k) {
  cells <- unlist(lapply(names(bases), function(b) cell(tab[tab$row == k & tab$base == b, ])))
  lab <- if (rows$type[k] == "Reference") paste0("\\textit{", rows$label[k], "}") else rows$label[k]
  paste(c(lab, cells), collapse = " & ")
})
grp <- function(t) sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", t)
body <- c(paste(lines[1], "\\\\ \\midrule"), grp("Network-based sets"), paste(lines[2:5], "\\\\"),
          "\\midrule", grp("Size-matched benchmarks"), paste(lines[6:11], "\\\\"), "\\midrule", paste(lines[12], "\\\\"))
base_auc <- sapply(names(bases), function(b) tab$AUC[tab$row == 1 & tab$base == b])
tex <- c("\\begin{table*}[htbp]",
  "\\caption{\\textbf{Discrimination of 10-year incident T2D in held-out data.} Logistic models fitted on the selection set (30,240 participants, 1,028 cases) and evaluated on the held-out 30\\% (12,984 participants, 432 cases). Each protein set is added to three base models: age and sex; the QDiabetes-2018 variables; and the QDiabetes-2018 variables plus HbA1c. AUC: area under the ROC curve with DeLong 95\\% confidence interval. LRT $p$: likelihood-ratio test of the protein terms against the base model, with both models refitted on the held-out data. DE: differential expression (proteins ranked by $|\\log \\mathrm{OR}|$ on the selection set). Network-based sets were defined on the T2D-specific network; benchmarks of the same size were selected on the same data.}",
  "\\label{tab:prediction}",
  "\\centering\\small",
  "\\begin{tabular}{lcccccc}", "\\toprule",
  " & \\multicolumn{2}{c}{Age + sex} & \\multicolumn{2}{c}{QDiabetes} & \\multicolumn{2}{c}{QDiabetes + HbA1c} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  "Protein set & AUC (95\\% CI) & LRT $p$ & AUC (95\\% CI) & LRT $p$ & AUC (95\\% CI) & LRT $p$ \\\\", "\\midrule",
  body, "\\bottomrule", "\\end{tabular}", "\\end{table*}")
tex <- gsub("ref\\.", "--", tex)
writeLines(tex, file.path(OUT_DIR, "B_table_paper.tex"))

# ---- Figure draft: same content -------------------------------------------------------
fig <- tab[!is.na(tab$AUC), ]
fig$Base <- factor(bases[fig$base], bases)
fig$label <- factor(fig$label, rev(rows$label))
fig$type <- factor(fig$type, c("Reference", "Network-based", "Differential expression", "Lasso path"))
fig$txt <- ifelse(is.na(fig$p_LRT), "", paste0("p = ", fmt_p(fig$p_LRT)))
ref <- data.frame(Base = factor(bases, bases), AUC = base_auc)
g <- ggplot(fig, aes(AUC, label, colour = type)) +
  geom_vline(data = ref, aes(xintercept = AUC), linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.25, linewidth = 0.5) +
  geom_text(aes(x = Inf, label = txt), hjust = 1.05, size = 2.3, colour = "grey20", show.legend = FALSE) +
  facet_wrap(~ Base, nrow = 1) +
  scale_colour_manual(values = c(Reference = "grey35", `Network-based` = "#d95f02", `Differential expression` = "#7570b3",
                                 `Lasso path` = "#1b9e77"), name = NULL) +
  scale_x_continuous(limits = c(0.62, 0.95), breaks = seq(0.65, 0.9, 0.05)) +
  labs(x = "Test AUC, held-out 30% (95% CI, DeLong); dashed: base model; text: LRT p vs base model", y = NULL,
       title = "10-year incident T2D: added value of protein sets") +
  theme_bw(base_size = 9) + theme(legend.position = "top", panel.grid.minor = element_blank())
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE)
ggsave(file.path(OUT_DIR, "figures", "Fig4b_prediction_by_base.png"), g, width = 10, height = 4.6, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "figures", "Fig4b_prediction_by_base.pdf"), g, width = 10, height = 4.6)
print(tab[, c("label", "base", "AUC", "lo", "hi", "dAUC", "p_DeLong", "p_LRT")], digits = 3, row.names = FALSE)
