#################################################
## A3. TIMP1 neighbourhood
## A4.1 Replication of the T2D OR edges in the held-out data (training + test)
## A4.2 STRING support of the T2D edges vs the permutation null
## Run from the private project folder (see README.md):
##   Rscript --vanilla downstream/A3_A4_TIMP1_replication_STRING.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))
suppressPackageStartupMessages({ use_prcggm()
                                 library(igraph); library(ggraph); library(ggrepel) })
covariate_design <- prCGGM:::covariate_design

edges <- edge_list(A_or)
maxT <- read.csv(file.path(PERM_DIR, "T2D_OR_edges_maxT.csv"))
key <- function(a, b) paste(pmin(a, b), pmax(a, b))
edges$p_maxT <- maxT$p_maxT[match(key(edges$protein_1, edges$protein_2), key(maxT$p1, maxT$p2))]
stopifnot(!anyNA(edges$p_maxT))
# Direction of change: slope in controls = baseline, in cases = baseline + T2D
edges$in_cases <- ifelse(abs(edges$baseline + edges$weight) > abs(edges$baseline), "stronger", "weaker")
edges$in_cases[edges$baseline == 0] <- "case-specific"

# ---- STRING prior ------------------------------------------------------------
S <- readMM("STRING_W"); dimnames(S) <- list(read.table("STRING_W_rownames.txt")$V1, read.table("STRING_W_colnames.txt")$V1)
S <- as.matrix(S[prot, prot])
edges$STRING_w <- S[cbind(match(edges$protein_1, prot), match(edges$protein_2, prot))]

# ---- A4.1 Replication in the held-out data ------------------------------------
# Rebuild prCGGM's node regressions on new data: residuals from the mean step,
# design [r_-j, U_c * r_-j], restricted to the support selected in the selection set.
U_of <- function(d) covariate_design(code_covariates(d), fit$covariates)
resid_of <- function(d, U) as.matrix(d[, prot]) - cbind(1, U) %*% t(fit$mean_effects[prot, ])
U_sel <- U_of(selection_data); R_sel <- resid_of(selection_data, U_sel)
cat("Check residuals on the selection set (max |diff| vs fit$residuals):",
    signif(max(abs(R_sel - fit$residuals)), 3), "\n")

ho <- rbind(train_data, test_data)
U <- U_of(ho); R <- resid_of(ho, U)
cat("Held-out data:", nrow(ho), "participants,", sum(U[, "Diab_Incident"]), "cases\n")
design_node <- function(j, cols) {        # columns named as in fit$coefficients
  X <- sapply(cols, function(cn) {
    parts <- strsplit(cn, ":", fixed = TRUE)[[1]]
    if (length(parts) == 1) R[, parts] else U[, parts[1]] * R[, parts[2]]
  })
  matrix(X, nrow = nrow(R), dimnames = list(NULL, cols))
}
# Directed T2D terms behind the OR edges
dir_terms <- which(B_t2d != 0, arr.ind = TRUE)
rep_dir <- do.call(rbind, lapply(sort(unique(dir_terms[, 1])), function(j) {
  cf <- fit$coefficients[j, ]; cols <- setdiff(names(cf)[cf != 0], "(Intercept)")
  X <- design_node(j, cols)
  m <- lm(R[, j] ~ X + U)            # covariate main effects kept as adjustment
  s <- summary(m)$coefficients
  ks <- dir_terms[dir_terms[, 1] == j, 2]
  do.call(rbind, lapply(ks, function(k) {
    term <- paste0("X", "Diab_Incident:", prot[k])
    data.frame(response = prot[j], predictor = prot[k], beta_sel = B_t2d[j, k],
               beta_ho = s[term, 1], se_ho = s[term, 2], p_ho = s[term, 4], support_size = length(cols))
  }))
}))
rep_dir$concordant <- sign(rep_dir$beta_ho) == sign(rep_dir$beta_sel)
cat("\nDirected T2D terms:", nrow(rep_dir), "| sign-concordant:", sum(rep_dir$concordant),
    sprintf("(%.0f%%), binomial p = %.2g", 100 * mean(rep_dir$concordant),
            binom.test(sum(rep_dir$concordant), nrow(rep_dir), alternative = "greater")$p.value),
    "| p < 0.05:", sum(rep_dir$p_ho < 0.05), "(expected", round(0.05 * nrow(rep_dir), 1), ")",
    "| p < 0.05 and concordant:", sum(rep_dir$p_ho < 0.05 & rep_dir$concordant), "\n")
write.csv(rep_dir, file.path(OUT_DIR, "A4_replication_directed_terms.csv"), row.names = FALSE)

# Edge level: the direction that defines the OR weight (larger |beta|)
pick <- function(a, b) {
  ab <- rep_dir[rep_dir$response == a & rep_dir$predictor == b, ]
  ba <- rep_dir[rep_dir$response == b & rep_dir$predictor == a, ]
  r <- rbind(ab, ba); r[which.max(abs(r$beta_sel)), ]
}
rep_edge <- do.call(rbind, Map(pick, edges$protein_1, edges$protein_2))
edges$beta_heldout <- rep_edge$beta_ho; edges$p_heldout <- rep_edge$p_ho; edges$concordant <- rep_edge$concordant
bt <- binom.test(sum(edges$concordant), nrow(edges), alternative = "greater")
cat("Edges:", nrow(edges), "| sign-concordant:", sum(edges$concordant), sprintf("(%.0f%%), binomial p = %.2g", 100 * mean(edges$concordant), bt$p.value),
    "| p < 0.05:", sum(edges$p_heldout < 0.05), "(expected 8.5) | p < 0.05 and concordant:",
    sum(edges$p_heldout < 0.05 & edges$concordant), "\n")
for (cl in list(c("AND edges", "in_AND"), c("TIMP1 edges", "TIMP1"))) {
  s <- if (cl[2] == "TIMP1") edges$protein_1 == "TIMP1" | edges$protein_2 == "TIMP1" else edges$in_AND
  cat(sprintf("  %s: %d edges, concordant %d, p < 0.05 %d\n", cl[1], sum(s), sum(edges$concordant[s]), sum(edges$p_heldout[s] < 0.05)))
}
# Calibration of the replication under the null: same test on the null edges of 5 permutations is
# not possible without their supports; instead, compare with a sign test on random pairs (see A4 notes).

# ---- A4.2 STRING support vs the permutation null ------------------------------
o <- readRDS(file.path(PERM_DIR, "original.rds"))
string_stats <- function(Bt) {
  diag(Bt) <- 0; Aor <- (Bt != 0) | (t(Bt) != 0); Aand <- (Bt != 0) & (t(Bt) != 0); up <- upper.tri(Aor)
  c(OR_n = sum(Aor[up]), OR_frac = mean(S[Aor & up] > 0), OR_meanW = mean(S[Aor & up]),
    AND_n = sum(Aand[up]), AND_frac = if (any(Aand[up])) mean(S[Aand & up] > 0) else NA)
}
fl <- sort(list.files(PERM_DIR, "^perm_[0-9]+_chunk"))
st <- do.call(rbind, lapply(fl, function(f) {
  x <- readRDS(file.path(PERM_DIR, f)); C <- matrix(0, p, ncol(o$coef)); C[x$nodes, ] <- as.matrix(x$coef)
  c(b = x$b, string_stats(C[, 1 + 3 * p + seq_len(p)]))
}))
st <- as.data.frame(st); obs <- st[st$b == 0, ]; nul <- st[st$b > 0, ]
panel_frac <- mean(S[upper.tri(S)] > 0)
cat("\nSTRING support (W > 0): panel pairs", sprintf("%.3f", panel_frac),
    "| observed OR edges", sprintf("%.3f (%d/%d)", obs$OR_frac, round(obs$OR_frac * obs$OR_n), obs$OR_n),
    "| null OR edges mean", sprintf("%.3f [95%%: %.3f-%.3f]", mean(nul$OR_frac), quantile(nul$OR_frac, .025), quantile(nul$OR_frac, .975)),
    "| p =", round((1 + sum(nul$OR_frac >= obs$OR_frac)) / (nrow(nul) + 1), 3), "\n")
cat("Mean STRING weight: observed", round(obs$OR_meanW, 4), "| null mean", round(mean(nul$OR_meanW), 4),
    "| p =", round((1 + sum(nul$OR_meanW >= obs$OR_meanW)) / (nrow(nul) + 1), 3), "\n")
cat("Pooled null STRING fraction (all null edges):", round(sum(nul$OR_frac * nul$OR_n) / sum(nul$OR_n), 3), "\n")
write.csv(st, file.path(OUT_DIR, "A4_STRING_vs_null.csv"), row.names = FALSE)

g <- ggplot(nul, aes(OR_frac)) + geom_histogram(bins = 25, fill = "grey70", colour = "white") +
  geom_vline(xintercept = obs$OR_frac, colour = "#d95f02", linewidth = 0.8) +
  geom_vline(xintercept = panel_frac, linetype = 2) +
  labs(x = "Fraction of T2D OR edges with STRING support", y = "Permutations",
       title = "STRING support of the T2D edges: observed (orange), null, panel (dashed)")
save_fig(g, "FigS_A4_STRING_vs_null", 6, 4)

# ---- A3. TIMP1 neighbourhood -----------------------------------------------------
tn <- edges[edges$protein_1 == "TIMP1" | edges$protein_2 == "TIMP1", ]
tn$partner <- ifelse(tn$protein_1 == "TIMP1", tn$protein_2, tn$protein_1)
tn$partner_class <- all_class <- ifelse(tn$partner %in% candidate_hubs, "Candidate hub",
                                 ifelse(tn$partner %in% rewired, "Rewired", "Other"))
tn <- tn[order(tn$p_maxT, -abs(tn$weight)), ]
cat("\nTIMP1: mean effect (prCGGM)", round(fit$mean_effects["TIMP1", "Diab_Incident"], 3),
    "| partners", nrow(tn), "| stronger in cases", sum(tn$in_cases == "stronger"),
    "| weaker", sum(tn$in_cases == "weaker"), "| case-specific", sum(tn$in_cases == "case-specific"),
    "| STRING-supported", sum(tn$STRING_w > 0), "| held-out concordant", sum(tn$concordant), "\n")
print(tn[, c("partner", "weight", "baseline", "in_cases", "p_maxT", "STRING_w", "beta_heldout", "p_heldout", "partner_class")],
      row.names = FALSE, digits = 3)
write.csv(tn, file.path(OUT_DIR, "A3_TIMP1_neighbourhood.csv"), row.names = FALSE)

gt <- graph_from_data_frame(data.frame(from = "TIMP1", to = tn$partner, weight = abs(tn$weight),
                                       change = tn$in_cases, string = tn$STRING_w > 0),
                            vertices = data.frame(name = c("TIMP1", tn$partner),
                                                  cls = c("Confirmed hub", tn$partner_class)), directed = FALSE)
set.seed(1)
g <- ggraph(gt, layout = "star", center = 1) +
  geom_edge_link(aes(edge_width = weight, edge_colour = change, edge_linetype = string)) +
  scale_edge_width(range = c(0.3, 2.5), name = "|T2D coefficient|") +
  scale_edge_colour_manual(values = c(stronger = "#b2182b", weaker = "#2166ac", `case-specific` = "grey35"),
                           name = "Dependence in cases") +
  scale_edge_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), name = "STRING support") +
  geom_node_point(aes(colour = cls), size = 4) +
  scale_colour_manual(values = c(`Confirmed hub` = "#d95f02", `Candidate hub` = "#7570b3", Rewired = "#1b9e77", Other = "grey60"), name = NULL) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA)) +
  labs(title = "TIMP1: T2D-specific neighbourhood (OR graph)",
       subtitle = "Edge colour: change of the partial regression slope in cases vs controls; dashed: no STRING support")
save_fig(g, "Fig_A3_TIMP1_neighbourhood", 8, 7)

write.csv(edges, file.path(OUT_DIR, "A4_edges_with_replication_STRING.csv"), row.names = FALSE)
