#################################################
## Shared set-up for the downstream analysis of the 10-year T2D fit
## (see README.md). Sourced by every downstream script, which is run from
## the private project folder, e.g.
##   Rscript --vanilla ../T2D_application_code/downstream/A1_global_null.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "config.R"))
suppressPackageStartupMessages({ library(Matrix); library(ggplot2) })

RES_DIR  <- Sys.getenv("T2D_RES", "results_bic_relaxed_10y")
PERM_DIR <- file.path(RES_DIR, "permutation_null")
OUT_DIR  <- file.path(RES_DIR, "downstream")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

load(file.path(RES_DIR, "PPI_T2D_app.RData"))   # fit, selection_data, train_data, test_data, ...
prot <- rownames(fit$coefficients)
p <- length(prot)

# Node-wise coefficient blocks (row = response node, column = predictor node)
coef_block <- function(prefix) {
  C <- fit$coefficients
  B <- if (prefix == "") C[, prot] else C[, paste0(prefix, ":", prot)]
  B <- as.matrix(B); dimnames(B) <- list(prot, prot); diag(B) <- 0; B
}
B_base <- coef_block("")
B_t2d  <- coef_block("Diab_Incident")

# T2D component: OR / AND adjacency; OR weight = coefficient with the larger |beta|
A_or  <- (B_t2d != 0) | (t(B_t2d) != 0)
A_and <- (B_t2d != 0) & (t(B_t2d) != 0)
W_or  <- ifelse(abs(B_t2d) >= abs(t(B_t2d)), B_t2d, t(B_t2d)); W_or[!A_or] <- 0
# Baseline coefficient in the same direction as the OR weight (for "stronger / weaker in cases")
Wb_or <- ifelse(abs(B_t2d) >= abs(t(B_t2d)), B_base, t(B_base)); Wb_or[!A_or] <- 0

edge_list <- function(A) {
  idx <- which(A & upper.tri(A), arr.ind = TRUE)
  data.frame(protein_1 = prot[idx[, 1]], protein_2 = prot[idx[, 2]],
             weight = W_or[idx], baseline = Wb_or[idx],
             in_AND = A_and[idx], stringsAsFactors = FALSE)
}

# Lists fixed in README.md §4 (recomputed here as a check)
hub <- read.csv(file.path(PERM_DIR, "T2D_hub_degree_test.csv"))
confirmed_hubs <- hub$protein[hub$p_maxT <= 0.05]
candidate_hubs <- hub$protein[hub$p_marg <= 0.01 & hub$obs >= 5]
rewired        <- hub$protein[hub$p_marg <= 0.01 & hub$obs >= 3]
stopifnot(identical(confirmed_hubs, "TIMP1"), length(candidate_hubs) == 12, length(rewired) == 21)
and_edges <- edge_list(A_and)[, c("protein_1", "protein_2")]
and_prot  <- sort(unique(unlist(and_edges)))

save_fig <- function(g, name, width, height) {
  for (ext in c("pdf", "png"))
    ggsave(file.path(OUT_DIR, paste0(name, ".", ext)), g, width = width, height = height,
           units = "in", dpi = 300)
  cat("Written", name, "\n")
}
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
                                          strip.background = element_rect(fill = "grey95")))
