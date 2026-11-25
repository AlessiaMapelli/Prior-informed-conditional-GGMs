#################################################
## Part B. Prediction of 10-year incident T2D (README.md §1, Part B).
## Design (decided 2026-11-23): protein sets are chosen AND the models are fitted on
## the selection set (70%, unmatched); evaluation on the whole held-out 30%
## (former training + test sets, a random unmatched sample of the cohort). Both are
## restricted to the common complete cases (all clinical variables incl. HbA1c and PRS).
## The first run (fit on the age/sex-matched training set) is kept as *_v1_matched.
## Run from the private project folder: Rscript --vanilla ../T2D_application_code/downstream/B_prediction.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))
suppressPackageStartupMessages({ library(pROC); library(glmnet) })
set.seed(123)

e <- new.env(); load(INPUT_RDATA, envir = e)
clin <- e$clinical_model_df
clin_vars <- e$clinical_var                         # QDiabetes-2018 variables (incl. sex, age)
all_clin  <- e$clinical_var_biomarkers_PRS          # + hba1c + PRS_T2D

prep <- function(d) {
  d <- merge(d[, c("id", "Diab", prot)], clin[, c("id", all_clin)], by = "id")
  d$y <- as.integer(d$Diab == "Incident")
  d[complete.cases(d[, c(prot, all_clin)]), ]
}
tr <- prep(selection_data); te <- prep(rbind(train_data, test_data))
cat("Training:", nrow(tr), "participants,", sum(tr$y), "cases | Test:", nrow(te), "participants,", sum(te$y), "cases\n")

# Standardise proteins with the training mean and SD
mu <- colMeans(tr[, prot]); s <- apply(tr[, prot], 2, sd)
for (j in prot) { tr[[j]] <- (tr[[j]] - mu[j]) / s[j]; te[[j]] <- (te[[j]] - mu[j]) / s[j] }
add_products <- function(d, pairs) { for (k in seq_len(nrow(pairs))) d[[paste0(pairs[k, 1], "_x_", pairs[k, 2])]] <- d[[pairs[k, 1]]] * d[[pairs[k, 2]]]; d }

# ---- Protein sets ------------------------------------------------------------
nodes_tab <- read.csv(file.path(OUT_DIR, "A2_node_table.csv"))
de <- read.csv(file.path(OUT_DIR, "A2_DE_selection_set.csv"))
and_pairs <- as.matrix(and_edges)
timp1_pairs <- cbind("TIMP1", c("ACY1", "ADH4", "DPP4", "SELE"))
stopifnot(all(timp1_pairs[, 2] %in% candidate_hubs))
central <- nodes_tab$protein[nodes_tab$top10_central]
de_all <- de$protein[de$q_BH < 0.05]
de_top <- function(k) de$protein[order(de$rank_abs_effect)][1:k]

# Lasso on the selection set: proteins standardised, age and sex unpenalised
sel <- selection_data
Xs <- cbind(scale(as.matrix(sel[, prot])), age = sel$Age, male = as.integer(sel$Sex == "Male"))
ys <- as.integer(sel$Diab == "Incident")
pf <- c(rep(1, length(prot)), 0, 0)
cvl <- cv.glmnet(Xs, ys, family = "binomial", penalty.factor = pf, nfolds = 5)
cf <- coef(cvl, s = "lambda.min")[-1, 1]
lasso_all <- intersect(names(cf)[cf != 0], prot)
path <- cvl$glmnet.fit$beta[prot, ]
entry <- apply(path != 0, 1, function(z) if (any(z)) which(z)[1] else Inf)
size_at_entry <- abs(sapply(prot, function(j) if (is.finite(entry[j])) path[j, entry[j]] else 0))
lasso_order <- prot[order(entry, -size_at_entry)]
lasso_top <- function(k) lasso_order[1:k]
cat("Lasso (selection set): lambda.min selects", length(lasso_all), "proteins | DE q < 0.05:", length(de_all), "\n")
cat("First 12 to enter the lasso path:", lasso_top(12), "\n")
cat("Top-12 DE by |log OR|:", de_top(12), "\n")
cat("First 21 in the lasso path:", lasso_top(21), "\nTop-21 DE:", de_top(21), "\n")

sets <- list(
  M0 = list(prot = character(0)),
  M1 = list(prot = and_prot),
  M2 = list(prot = and_prot, pairs = and_pairs),
  M3a = list(prot = "TIMP1"),
  M3b = list(prot = candidate_hubs),
  M3c = list(prot = rewired),
  M3d = list(prot = candidate_hubs, pairs = timp1_pairs),
  M3e = list(prot = central),
  M4 = list(prot = de_all),
  M4b_k7 = list(prot = de_top(7)), M4b_k12 = list(prot = de_top(12)),
  M4b_k15 = list(prot = de_top(15)), M4b_k21 = list(prot = de_top(21)),
  M8 = list(prot = lasso_all),
  M8b_k7 = list(prot = lasso_top(7)), M8b_k12 = list(prot = lasso_top(12)),
  M8b_k15 = list(prot = lasso_top(15)), M8b_k21 = list(prot = lasso_top(21)))
writeLines(unlist(lapply(names(sets), function(m) paste0(m, ": ", paste(c(sets[[m]]$prot,
  if (!is.null(sets[[m]]$pairs)) paste0(sets[[m]]$pairs[, 1], "x", sets[[m]]$pairs[, 2])), collapse = ", ")))),
  file.path(OUT_DIR, "B_protein_sets.txt"))

# ---- Model fitting and evaluation --------------------------------------------
base_terms <- list(M0 = c("age", "sex"), M5 = clin_vars, M6 = c(clin_vars, "hba1c"), M7 = c(clin_vars, "hba1c", "PRS_T2D"))
fits <- list(); fits_te <- list(); preds <- list()
fit_model <- function(name, base, set) {
  pairs <- set$pairs; d_tr <- tr; d_te <- te
  prod <- character(0)
  if (!is.null(pairs)) { d_tr <- add_products(d_tr, pairs); d_te <- add_products(d_te, pairs); prod <- paste0(pairs[, 1], "_x_", pairs[, 2]) }
  terms <- c(base_terms[[base]], set$prot, prod)
  f <- reformulate(terms, "y")
  m <- suppressWarnings(glm(f, family = binomial, data = d_tr, control = glm.control(maxit = 100)))
  p_tr <- predict(m, d_tr, type = "response"); p_te <- predict(m, d_te, type = "response")
  r_tr <- roc(d_tr$y, p_tr, quiet = TRUE, direction = "<"); r_te <- roc(d_te$y, p_te, quiet = TRUE, direction = "<")
  thr <- coords(r_tr, "best", best.method = "youden", ret = "threshold", transpose = FALSE)$threshold[1]
  ci <- ci.auc(r_te, method = "delong")
  # Same model refitted on the held-out data (for the test-set LRT, valid because the sets are fixed)
  m_te <- suppressWarnings(glm(f, family = binomial, data = d_te, control = glm.control(maxit = 100)))
  fits[[name]] <<- m; fits_te[[name]] <<- m_te; preds[[name]] <<- r_te
  data.frame(model = name, base = base, n_proteins = length(set$prot), n_products = length(prod),
             n_param = length(coef(m)) - 1, converged = m$converged,
             separation = any(p_tr < 1e-8 | p_tr > 1 - 1e-8),
             AUC_train = as.numeric(auc(r_tr)), AUC_test = ci[2], AUC_test_lo = ci[1], AUC_test_hi = ci[3],
             threshold_train = thr, sens_test = mean(p_te[d_te$y == 1] >= thr), spec_test = mean(p_te[d_te$y == 0] < thr))
}
res <- list()
for (m in names(sets)) res[[m]] <- fit_model(m, "M0", sets[[m]])
for (b in c("M5", "M6", "M7")) res[[b]] <- fit_model(b, b, list(prot = character(0)))
for (b in c("M5", "M6")) for (m in c("M2", "M3a", "M3b", "M3c", "M3e", "M4b_k7", "M4b_k12", "M4b_k15", "M4b_k21", "M8b_k7", "M8b_k12", "M8b_k15", "M8b_k21"))
  res[[paste0(b, "+", m)]] <- fit_model(paste0(b, "+", m), b, sets[[m]])
res <- do.call(rbind, res); rownames(res) <- NULL
print(res[, c("model", "n_proteins", "n_param", "converged", "separation", "AUC_train", "AUC_test", "AUC_test_lo", "AUC_test_hi", "sens_test", "spec_test")],
      row.names = FALSE, digits = 3)
write.csv(res, file.path(OUT_DIR, "B_model_performance.csv"), row.names = FALSE)

# ---- Comparisons ---------------------------------------------------------------
cmp <- function(small, big, nested = TRUE) {
  dl <- roc.test(preds[[big]], preds[[small]], method = "delong", paired = TRUE)
  lr <- function(fl) { a <- anova(fl[[small]], fl[[big]], test = "LRT"); c(a$Deviance[2], a$Df[2], a$`Pr(>Chi)`[2]) }
  lrt <- if (nested) lr(fits) else c(NA, NA, NA)
  lrt_te <- if (nested) lr(fits_te) else c(NA, NA, NA)
  data.frame(model = big, reference = small, delta_AUC_test = as.numeric(auc(preds[[big]]) - auc(preds[[small]])),
             p_DeLong = dl$p.value, LRT_df = lrt[2], LRT_chisq_train = lrt[1], p_LRT_train = lrt[3],
             LRT_chisq_test = lrt_te[1], p_LRT_test = lrt_te[3])
}
comparisons <- rbind(
  # products and network sets vs age + sex
  cmp("M0", "M1"), cmp("M1", "M2"), cmp("M0", "M3a"), cmp("M0", "M3b"), cmp("M0", "M3c"),
  cmp("M3b", "M3d"), cmp("M0", "M3e"),
  # network sets vs benchmark selection at matched size (not nested)
  cmp("M4b_k7", "M1", FALSE), cmp("M8b_k7", "M1", FALSE), cmp("M4b_k7", "M2", FALSE), cmp("M8b_k7", "M2", FALSE),
  cmp("M4b_k12", "M3b", FALSE), cmp("M8b_k12", "M3b", FALSE), cmp("M4b_k12", "M3d", FALSE), cmp("M8b_k12", "M3d", FALSE),
  cmp("M4b_k15", "M3e", FALSE), cmp("M8b_k15", "M3e", FALSE), cmp("M4b_k21", "M3c", FALSE), cmp("M8b_k21", "M3c", FALSE),
  # same, on top of the clinical models
  do.call(rbind, lapply(c("M5", "M6"), function(b) rbind(
    cmp(paste0(b, "+M4b_k12"), paste0(b, "+M3b"), FALSE), cmp(paste0(b, "+M8b_k12"), paste0(b, "+M3b"), FALSE),
    cmp(paste0(b, "+M4b_k15"), paste0(b, "+M3e"), FALSE), cmp(paste0(b, "+M8b_k15"), paste0(b, "+M3e"), FALSE),
    cmp(paste0(b, "+M4b_k21"), paste0(b, "+M3c"), FALSE), cmp(paste0(b, "+M8b_k21"), paste0(b, "+M3c"), FALSE)))),
  cmp("M4", "M3b", FALSE), cmp("M8", "M3b", FALSE), cmp("M4", "M3c", FALSE), cmp("M8", "M3c", FALSE),
  # incremental value over the clinical benchmark
  cmp("M5", "M6"), cmp("M6", "M7"),
  do.call(rbind, lapply(c("M5", "M6"), function(b) do.call(rbind, lapply(c("M2", "M3a", "M3b", "M3c", "M3e", "M4b_k7", "M4b_k12", "M4b_k15", "M4b_k21", "M8b_k7", "M8b_k12", "M8b_k15", "M8b_k21"),
    function(m) cmp(b, paste0(b, "+", m)))))))
print(comparisons, row.names = FALSE, digits = 3)
write.csv(comparisons, file.path(OUT_DIR, "B_model_comparisons.csv"), row.names = FALSE)

# Every model against its base model (for the paper table and figure): protein sets
# vs age + sex (M0); M5 vs M0; M6 vs M5; M7 vs M6; M5/M6 + set vs M5/M6. All nested.
stopifnot(all(c("age", "sex") %in% clin_vars))
base_of <- function(m, b) switch(m, M0 = NA, M5 = "M0", M6 = "M5", M7 = "M6", b)
vs_base <- do.call(rbind, lapply(seq_len(nrow(res)), function(k) {
  r <- base_of(res$model[k], res$base[k]); if (is.na(r)) NULL else cmp(r, res$model[k]) }))
write.csv(vs_base, file.path(OUT_DIR, "B_model_vs_base.csv"), row.names = FALSE)

# Products: coefficients and Wald tests in M2 and M3d
for (m in c("M2", "M3d")) { s <- summary(fits[[m]])$coefficients; cat("\n", m, "product terms:\n"); print(s[grepl("_x_", rownames(s)), , drop = FALSE], digits = 3) }

# ---- Figure: test AUC with 95% CI ----------------------------------------------
res$group <- ifelse(res$base == "M0", ifelse(grepl("^M[123]", res$model), "Network-based", ifelse(res$model == "M0", "Age + sex", "Benchmark selection")),
                    "Clinical (QDiabetes) based")
res$model <- factor(res$model, levels = rev(res$model))
g <- ggplot(res, aes(AUC_test, model, colour = group)) +
  geom_vline(xintercept = res$AUC_test[res$model == "M0"], linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = AUC_test_lo, xmax = AUC_test_hi), size = 0.25) +
  geom_text(aes(label = paste0(n_proteins, " prot.")), x = Inf, hjust = 1.1, size = 2.3, colour = "grey30") +
  scale_colour_manual(values = c(`Network-based` = "#d95f02", `Benchmark selection` = "#7570b3", `Age + sex` = "grey40", `Clinical (QDiabetes) based` = "#1b9e77"), name = NULL) +
  labs(x = "Test AUC (95% CI, DeLong)", y = NULL, title = "10-year incident T2D: test-set discrimination") +
  theme(legend.position = "top")
save_fig(g, "Fig_B_AUC", 7, 7.5)
