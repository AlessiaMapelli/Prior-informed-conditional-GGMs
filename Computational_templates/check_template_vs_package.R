#################################################
## Check that this template and the prCGGM package give the same networks
#################################################
# Runs GGReg_full_estimation() (template) and prcggm() (package) on the
# example data with the same seed, for every selection rule shared by the two
# implementations and both symmetrisations, and compares the adjacency
# matrices of every component.
#
# Usage (from the Computational_templates folder; needs prCGGM >= 0.2.0):
#   Rscript --vanilla check_template_vs_package.R
# Optional: PRCGGM_LIB=<library with prCGGM 0.2> to use a non-default library.
#
# Expected result (2026-11-25, prCGGM v0.2 commit e408e8d):
#   - all deterministic rules (bic_relaxed, ebic_relaxed, bic, ebic, cv_min,
#     cv_1se) give the same supports, with coefficients equal to ~1e-15;
#   - stability_bic: identical with AND on the example data, may differ by an
#     edge with OR, because the random subsamples are drawn differently (the
#     package draws one set shared by all nodes, the template one set per
#     node); the estimator is the same.

lib <- Sys.getenv("PRCGGM_LIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(prCGGM))
if (packageVersion("prCGGM") < "0.1.0.9000") {
  stop("prCGGM ", packageVersion("prCGGM"), " found: the selection rules need version >= 0.2.0")
}
cat("prCGGM", as.character(packageVersion("prCGGM")), "from", find.package("prCGGM"), "\n")
suppressPackageStartupMessages(source("ggReg_main_functions.R"))
load("Example_data.RData")  # expression_data, covariate_data, ppi_weight_matrix

rules <- c("bic_relaxed", "ebic_relaxed", "bic", "ebic", "cv_min", "cv_1se", "stability_bic")
edges <- function(m) { m <- abs(as.matrix(m)) > 1e-10; diag(m) <- FALSE; unname(m) }
res <- list()
for (sel in rules) for (sy in c("AND", "OR")) {
  set.seed(2024)
  tp <- suppressWarnings(GGReg_full_estimation(
    x = expression_data, known_ppi = ppi_weight_matrix, covariates = covariate_data,
    selection = sel, symm_method = sy, output_path = paste0(tempdir(), "/")))
  set.seed(2024)
  pk <- suppressWarnings(prcggm(expression_data, covariate_data, ppi_weight_matrix,
                                selection = sel, symm_method = sy))
  A <- tp$results$Dic_adj_matrics; B <- pk$networks
  stopifnot(setequal(names(A), names(B)))
  res[[length(res) + 1]] <- data.frame(
    selection = sel, symm = sy,
    same_support = all(sapply(names(A), function(k) identical(edges(A[[k]]), edges(B[[k]])))),
    max_abs_diff = max(sapply(names(A), function(k) max(abs(as.matrix(A[[k]]) - as.matrix(B[[k]]))))),
    edges_template = sum(sapply(A, function(m) sum(edges(m)) / 2)),
    edges_package = sum(sapply(B, function(m) sum(edges(m)) / 2)))
}
res <- do.call(rbind, res)
print(res, row.names = FALSE)
det <- res$selection != "stability_bic"
if (all(res$same_support[det]) && max(res$max_abs_diff[det]) < 1e-8) {
  cat("\nOK: template and package agree for all deterministic rules.\n")
} else {
  stop("Template and package differ for a deterministic rule.")
}
