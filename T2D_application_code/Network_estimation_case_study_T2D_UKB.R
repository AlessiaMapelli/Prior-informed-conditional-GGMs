###############################################################################
# T2D case study (UK Biobank proteomics) - Step 1: network estimation
#
# Prior-informed conditional GGM fitted with prCGGM (>= 0.2.0) using the
# selection guideline of the simulation study (P/n < 0.5, n >= 5000 ->
# selection = "bic_relaxed", symm = "AND"); p = 366 proteins, q = 3 covariates
# (Sex, Age, incident T2D) -> P = (p - 1)(q + 1) = 1,460 predictors per node.
#
# Endpoint (setting ENDPOINT below):
#   "10y" (current): 10-year incident T2D (surv_status, surv_time capped at
#         3,650 days); participants censored before 10 years (death / lost to
#         follow-up without T2D) are dropped; 70% selection set
#         (n = 32,439, 1,112 cases); outputs in results_bic_relaxed_10y/.
#         Primary analysis of the paper.
#   "5y"  5-year incident T2D (surv_status_5y); 50% selection set
#         (n = 24,446, 270 cases); outputs in results_bic_relaxed/.
#         Sensitivity analysis (Supplementary): T2D component not above the
#         permutation null.
#
# Run from the private project folder (data and results), with the code in
# T2D_CODE_DIR (see config.R and README.md):
#   Rscript --vanilla ../T2D_application_code/Network_estimation_case_study_T2D_UKB.R
# The script runs in two passes on the cluster login node (or an interactive
# session on a compute node):
#   1st run: data preparation, cohort split, mean step and SLURM submission
#            (one array task per protein).
#   2nd run: (once all 366 tasks have finished) collection of the node results
#            and saving of <OUT_DIR>/PPI_T2D_app.RData.
# The pass is chosen automatically from the files present in OUT_DIR.
###############################################################################

rm(list = ls(all = TRUE))
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "config.R"))

# ---------------------------------------------------------------------------
# 0. Settings
# ---------------------------------------------------------------------------
# prCGGM >= 0.2.0 is needed (selection = "bic_relaxed"); set PRCGGM_LIB to
# use a personal library (config.R).
ENDPOINT   <- Sys.getenv("T2D_ENDPOINT", "10y")      # "10y" (primary) or "5y"
EP <- switch(ENDPOINT,
  "10y" = list(status = "surv_status", time = "surv_time", horizon = 3650,
               drop_censored = TRUE, p_sel = 0.7,
               out_dir = "results_bic_relaxed_10y",
               # ID file to reproduce (written by the first run of this setting)
               ref_ids = file.path("results_bic_relaxed_10y", "Selection_data_id.txt")),
  "5y"  = list(status = "surv_status_5y", time = NULL, horizon = NULL,
               drop_censored = FALSE, p_sel = 0.5,
               out_dir = "results_bic_relaxed",
               ref_ids = "Selection_data_id.txt"),
  stop("ENDPOINT must be \"10y\" or \"5y\""))
OUT_DIR    <- EP$out_dir                             # all outputs of this version
NODE_DIR   <- file.path(OUT_DIR, "prcggm_nodes")     # SLURM inputs/node results
SBATCH     <- file.path(CODE_DIR, "slurm_prcggm_node.sbatch")
RUN_CONFOUNDER_ANALYSIS <- FALSE   # exploratory choice of the covariates (section 2), run once
SUBMIT     <- FALSE                 # FALSE: only print the sbatch command
SEED_SPLIT <- 123
SEED_FIT   <- 2026

use_prcggm()
library(Matrix)   # readMM()
library(dplyr)
cat("prCGGM version:", as.character(packageVersion("prCGGM")),
    "from", find.package("prCGGM"), "\n")

source(file.path(CODE_DIR, "helper_functions.R"))
dir.create(OUT_DIR, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------------
load(INPUT_RDATA)   # built by preprocessing/Case_specific_data_preprocessing_T2D_app_case.R
stopifnot(length(full_omic_var) == 366, !anyDuplicated(model_data$id))

# ---------------------------------------------------------------------------
# 2. Choice of the covariates (exploratory, run once)
#    a) factors explaining systematic proteomic variation (PC R^2), on the 366
#       proteins that enter the network;
#    b) factors associated with 5-year incident T2D.
#    Selected: Sex, Age, incident T2D.
# ---------------------------------------------------------------------------
potential_conf_factors <- c(
  "sex", "Age_at_recruitment", "BMI", "Smoking_status",
  "Cholesterol_lowering_medication", "Fasting_time",
  "Month_blood_sample_collected", "Time_of_day_blood_sample_collected",
  "C_reactive_protein", "genPC1", "genPC2", "genPC3", "genPC4", "genPC5")

if (RUN_CONFOUNDER_ANALYSIS) {
  results_conf_analysis_df <- perform_confounder_analysis(
    prot_data[, c("id", full_omic_var)], pheno_data, col_id = "id",
    conf_factor = potential_conf_factors)
  results_t2d_conf <- perform_t2d_confounder_analysis(
    model_data = model_data, pheno_data = pheno_data, col_id = "id",
    conf_factor = potential_conf_factors, outcome = EP$status,
    adjust_vars = c("sex", "Age_at_recruitment"),
    output_dir = "T2D_confounder_analysis")
}

# ---------------------------------------------------------------------------
# 3. Study cohort and split (EP$p_sel selection / matched training / test)
#    Complete cases are defined on the variables of the network model only;
#    the clinical variables are carried along so that split_dataset() can move
#    test participants with categorical levels unseen in training.
#    10y: participants censored before the horizon without T2D are dropped
#    (as the 5-year definition did, where they were NA).
# ---------------------------------------------------------------------------
cohort_data <- model_data
if (EP$drop_censored) {
  censored <- cohort_data[[EP$status]] == 0 & cohort_data[[EP$time]] < EP$horizon
  cat("Dropping", sum(censored), "participant(s) censored before",
      EP$horizon, "days without T2D\n")
  cohort_data <- cohort_data[!censored, ]
}
clin_extra <- setdiff(clinical_var, c("id", "sex", "Age_at_recruitment",
                                      "surv_status_5y", "surv_status",
                                      "surv_time", full_omic_var))
study_cohort <- cohort_data[, c("id", "sex", "Age_at_recruitment",
                                EP$status, full_omic_var, clin_extra)]
names(study_cohort)[2:4] <- c("Sex", "Age", "Diab")
study_cohort <- study_cohort[stats::complete.cases(
  study_cohort[, c("id", "Sex", "Age", "Diab", full_omic_var)]), ]
stopifnot(is.logical(study_cohort$Diab) || all(study_cohort$Diab %in% c(0, 1)))
study_cohort$Diab <- factor(ifelse(as.logical(study_cohort$Diab),
                                   "Incident", "Non diabetic"),
                            levels = c("Non diabetic", "Incident"))
study_cohort$Sex <- as.factor(study_cohort$Sex)
cat("Endpoint ", ENDPOINT, " (", EP$status, "): ", nrow(study_cohort),
    " complete cases, ", sum(study_cohort$Diab == "Incident"), " incident\n",
    sep = "")

# Reference selection set (read before split_dataset() rewrites the ID files
# in OUT_DIR): reruns must reproduce the split of the first run.
old_sel <- if (file.exists(EP$ref_ids)) read.table(EP$ref_ids)$V1

list_df <- split_dataset(
  data = study_cohort, sex_col = "Sex", age_col = "Age",
  endpoint_col = "Diab", incident_name = "Incident",
  non_incident_name = "Non diabetic", p_sel = EP$p_sel, p_case_train = 0.5,
  clinical_var = clin_extra, controls_per_case = 4, seed = SEED_SPLIT,
  output_dir = OUT_DIR)

if (!is.null(old_sel) && !setequal(old_sel, list_df$selection_data$id)) {
  stop("Selection set differs from ", EP$ref_ids, ": check the data and the ",
       "seed before fitting.")
}

summary_table <- summary_cohort(list_df, endpoint_col = "Diab",
                                incident_name = "Incident",
                                sex_col = "Sex", age_col = "Age")
print(summary_table)
write.csv(summary_table, file.path(OUT_DIR, "Cohort_summary.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Covariate coding (selection set)
#    - factors: reference level = most prevalent category, so that the
#      baseline network refers to the largest group (non-diabetic, majority
#      sex);
#    - Age: centred at the selection-set mean. prCGGM additionally divides
#      numeric covariates by their SD, so the Age component is the change in
#      the network per 1 SD of age and the baseline is the network at the mean
#      age.
# ---------------------------------------------------------------------------
selection_data <- list_df$selection_data
train_data     <- list_df$train_data
test_data      <- list_df$test_data

ref_most_prevalent <- function(f) {
  f <- droplevels(as.factor(f))
  stats::relevel(f, ref = names(which.max(table(f))))
}
covariate_levels <- list(
  Sex  = levels(ref_most_prevalent(selection_data$Sex)),
  Diab = levels(ref_most_prevalent(selection_data$Diab)))
stopifnot(covariate_levels$Diab[1] == "Non diabetic")
age_mean <- mean(selection_data$Age)

code_covariates <- function(d) {
  data.frame(Sex   = factor(d$Sex,  levels = covariate_levels$Sex),
             Age_c = d$Age - age_mean,
             Diab  = factor(d$Diab, levels = covariate_levels$Diab))
}
covariates <- code_covariates(selection_data)
stopifnot(!anyNA(covariates))
cat("Reference levels: Sex =", covariate_levels$Sex[1],
    "| Diab =", covariate_levels$Diab[1],
    "| Age centred at", round(age_mean, 2), "years\n")
print(table(Sex = covariates$Sex, Diab = covariates$Diab))

expression_data <- as.matrix(selection_data[, full_omic_var])

# ---------------------------------------------------------------------------
# 5. Prior knowledge: STRING confidence scores in [0, 1]
# ---------------------------------------------------------------------------
STRING_prior <- readMM("STRING_W")
dimnames(STRING_prior) <- list(read.table("STRING_W_rownames.txt")$V1,
                               read.table("STRING_W_colnames.txt")$V1)
stopifnot(all(full_omic_var %in% rownames(STRING_prior)),
          all(full_omic_var %in% colnames(STRING_prior)))
STRING_prior <- as.matrix(STRING_prior[full_omic_var, full_omic_var])
stopifnot(identical(colnames(expression_data), colnames(STRING_prior)),
          isSymmetric(unname(STRING_prior)),
          min(STRING_prior) >= 0, max(STRING_prior) <= 1)
cat("Prior-supported pairs:",
    sum(STRING_prior[upper.tri(STRING_prior)] > 0), "of",
    choose(length(full_omic_var), 2), "\n")

# ---------------------------------------------------------------------------
# 6. Estimation: mean step + one SLURM array task per protein
#    Defaults of prCGGM v.0.2: selection = "bic_relaxed" (lambda, alpha, theta
#    jointly by the BIC of the least-squares refit), refit = TRUE,
#    alpha_grid = c(0.75, 0.5, 0.9, 0.99), theta_grid = c(1, 0.8, 1.1, 1.3)
#    (full grid, no random search), symm_method = "AND"; mean step by lasso
#    with 5-fold CV (lambda.min). All given explicitly for the record.
# ---------------------------------------------------------------------------
input_file <- file.path(NODE_DIR, "prcggm_input.rds")
node_files <- file.path(NODE_DIR, sprintf("prcggm_node_%d.rds",
                                          seq_along(full_omic_var)))

if (!file.exists(input_file)) {
  # ---- 1st run: prepare and submit -----------------------------------------
  set.seed(SEED_FIT)
  job <- prcggm_slurm_prepare(
    x = expression_data,
    covariates = covariates,
    prior = STRING_prior,
    dir = NODE_DIR,
    screening = FALSE,
    mean_estimation = TRUE,
    lambda_mean_type = "min",
    selection = "bic_relaxed",
    refit = TRUE,
    tune = TRUE,
    alpha_grid = c(0.75, 0.5, 0.9, 0.99),
    theta_grid = c(1, 0.8, 1.1, 1.3),
    random_search = FALSE,
    nfolds = 5,
    symm_method = "AND",
    copy_templates = TRUE,       # copies prcggm_node_job.R (and a template sbatch)
    verbose = TRUE)
  dir.create("logs", showWarnings = FALSE)
  # Our sbatch (cluster R setup + optional personal library), run from the
  # project folder; logs go to ./logs.
  # NOTE: ignore the command printed above by prcggm_slurm_prepare() (it uses
  # the package's generic template, without the singularity R and the
  # personal library); use this one, which runs from the project folder.
  cmd <- sprintf("cd %s && sbatch --array=1-%d %s %s %s", getwd(), job$p,
                 SBATCH, normalizePath(job$input), normalizePath(NODE_DIR))
  cat("\nSubmit with (NOT the command above):\n  ", cmd, "\n")
  if (SUBMIT) system(cmd)
  cat("\nRe-run this script when all", job$p, "tasks have completed.\n")

} else if (!all(file.exists(node_files))) {
  # ---- tasks still running (or failed) ------------------------------------
  missing_nodes <- which(!file.exists(node_files))
  cat(length(missing_nodes), "node(s) not finished:",
      paste(head(missing_nodes, 30), collapse = ", "), "\n",
      "Resubmit failed ones with e.g.\n  cd ", getwd(), " && sbatch --array=",
      paste(missing_nodes, collapse = ","), " ", SBATCH, " ",
      normalizePath(input_file), " ", normalizePath(NODE_DIR), "\n", sep = "")

} else {
  # ---- 2nd run: collect and save -------------------------------------------
  fit <- prcggm_slurm_collect(NODE_DIR)
  print(fit)
  print(summary(fit))

  # Diagnostics for the report / paper
  tun <- fit$tuning
  cat("\nSelected predictors per node (of", (fit$p - 1) * (fit$q + 1), "):\n")
  print(summary(tun$df))
  cat("Selected alpha:\n");  print(table(tun$alpha))
  cat("Selected theta:\n");  print(table(tun$theta))
  cat("Node time (min):\n"); print(summary(tun$time))
  n_edges <- sapply(fit$networks, function(m) sum(m[upper.tri(m)] != 0))
  cat("Edges per component:\n"); print(n_edges)
  base_edges <- fit$networks$Baseline[upper.tri(fit$networks$Baseline)] != 0
  prior_up   <- STRING_prior[upper.tri(STRING_prior)] > 0
  cat("Baseline edges with STRING support:", sum(base_edges & prior_up), "of",
      sum(base_edges), "\n")
  write.csv(tun, file.path(OUT_DIR, "Node_tuning.csv"), row.names = FALSE)
  write.csv(data.frame(component = names(n_edges), edges = n_edges,
                       density = n_edges / choose(fit$p, 2)),
            file.path(OUT_DIR, "Edges_per_component.csv"), row.names = FALSE)

  # Objects used by the downstream reports
  estimation_info <- list(
    prcggm_version = as.character(packageVersion("prCGGM")),
    endpoint = ENDPOINT, endpoint_def = EP,
    selection = "bic_relaxed", symm_method = "AND",
    covariate_levels = covariate_levels, age_mean = age_mean,
    seed_split = SEED_SPLIT, seed_fit = SEED_FIT, date = Sys.Date(),
    session = sessionInfo())
  save(fit, selection_data, train_data, test_data, covariates,
       code_covariates, covariate_levels, age_mean, estimation_info,
       file = file.path(OUT_DIR, "PPI_T2D_app.RData"))
  cat("Saved", file.path(OUT_DIR, "PPI_T2D_app.RData"), "\n")
}
