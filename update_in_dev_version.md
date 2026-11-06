# Changes in the development version of the paper repository

This file tracks every change made to the repository since the last commit on `main`:

- **Reference commit:** `fbc43aa` ("Addition to gitignore").
- **Matching package release:** the code corresponds to the CRAN submission of `prCGGM` 0.1.0.
- **Status:** committed on branch `dev` (2026-11-06; not pushed). The matching package changes are on the package branch `dev-selection-rules` (see `Rpackages/prCGGM/update_in_dev_version.md`).
- **Detail:** the READMEs referenced below hold the full findings, tables and run logs.

Last updated: 2026-11-06.

---

## Summary

| # | Area | What changed | Effect on published results |
|---|---|---|---|
| 1 | Template copies in `Simulation_analysis/` and `T2D_application_code/` | Synced to the committed template, which carries the fixes of commit `bb80ca1` | Simulation: TPR/FPR/F1 identical, MagCor ±0.001 in 7 of 72 cells. T2D: unchanged |
| 2 | `Simulation_analysis/` | Results regenerated with the δ fix; `scr_matrix` added to the saved inputs | As in 1 |
| 3 | `T2D_application_code/` | Covariate network name `DiabIncident` → `Diab_Incident` | None (naming only) |
| 4 | `Computational_templates/ggReg_main_functions.R` | Selection rules, OLS refit, unified node fitting | New default rule; previous behaviour reproducible |
| 5 | `Computational_templates/ggReg_main_functions.R` | Relaxed BIC/EBIC, new default `ebic_relaxed` | New default rule |
| 5b | `Computational_templates/ggReg_main_functions.R` | Default grids aligned with the simulations | Default fits change (α 0.95 → 0.99, ϑ 1.5 dropped) |
| 6 | `Computational_templates/README.md` | Documents 4, 5 and 5b | — |
| 7 | `Add_simulation_analysis/` (git-ignored) | Exploratory and sensitivity simulations (Supplementary Part 1) | New supplementary material |
| 8 | `Simulation_analysis_v2/` (new) | Full simulation study and selection test (main text and Supplementary Part 2) | New main-text simulation |

---

## 1. Template copies synced to the committed template (2026-11-02)

`Simulation_analysis/ggReg_main_functions.R` and `T2D_application_code/ggReg_main_functions.R` ran on template commit `4cfcecd`; `T2D_application_prot/` ran on `330a36e`. They now equal the template of commit `bb80ca1` (the committed `Computational_templates` version), which includes these fixes:

- **Precision scale:** δ_jk = −β_jk / σ²_j. The old code used −σ_j·β_jk (`No_sim_Delta_hat_row`).
- **Symmetrisation:** the OR/AND choice is made on the regression coefficients, not on δ (`build_delta_dictionary(..., ref = Prec_reg_matrix)`, Eq. 18 of the paper).
- **Predictor selection by position** (`node_columns()`) instead of by column-name patterns.
  - Name patterns could select the wrong columns when protein names are prefixes of each other.
  - For T2D the selected predictors are identical.
- **Screening:**
  - optional user-supplied `scr_matrix`;
  - the γ threshold is the 20th percentile of the off-diagonal correlations only;
  - prior-supported pairs are never screened out.
- **`known_ppi` checks:** negative values are set to 0 and values > 1 are rescaled (with a warning).
- **Mean step:** with a fixed `lambda_mean`, the coefficients come from `coef(mod, s = lambda_mean)` (intercept included). The old code used `relax = TRUE` and `mod$beta`, which dropped the intercept.
- **Plots:** symmetric colour range, and no `corrplot` error on all-zero matrices.

Checks (`Add_simulation_analysis/README.md`, findings §1):
- T2D networks are unchanged by these fixes. Node 279 refitted with the current template matches the saved fit to about 1e-16.
- The downstream T2D analyses use only `Dic_adj_matrics` and `Cov_effect`.

## 2. `Simulation_analysis/` (paper simulation)

- **`Sim_functions.R`:** `scr_matrix = NULL` is now saved in `input_data_nodes.rda`; the current node code expects it.
- **New `Sim_regenerate_results.R` and `slurm_regenerate_results.sbatch`** (job 45334610, 2026-11-03).
  - They re-collect all 240 replications from the saved node files with δ = −β/σ² and the OR/AND choice on β. Nothing is refitted.
  - TPR, FPR, F1 and accuracy are identical to before. 7 of 72 MagCor cells change by at most ±0.001.
  - Previous outputs: `Simulation_results/archive_before_delta_fix/` (git-ignored).
- **Known problems of the paper's generator** (not changed here, documented for the revision; `Add_simulation_analysis/README.md` §3):
  - `sample_pa(power = 2.5)` gives a star network;
  - no subject's precision matrix is positive definite without a subject-specific diagonal correction;
  - the covariate matrices carry a diagonal.

  The new simulation (8) replaces this generator.

## 3. `T2D_application_code/`

- The current code names covariate networks `<variable>_<level>`, so the T2D incidence network is `Diab_Incident`, not `DiabIncident`.
- Updated in `Network_description_case_study_T2D_UKB.Rmd` (parameter `outcome_of_interest`), `enrichment_communities.py` and `README.md`.
- `T2D_application_prot/` (private, git-ignored) keeps `DiabIncident`, because its saved results use that name.

## 4. Template: selection rules and OLS refit (step 4, 2026-11-04)

File: `Computational_templates/ggReg_main_functions.R`. A copy is in `Simulation_analysis_v2/`.

- **New arguments** of `GGReg_full_estimation()`, `GGReg_cov_estimation_SLURM()`, `GGReg_cov_estimation_sequential()` and `GGReg_cov_estimation_parallel()`:
  - `selection` (rule for λ and the hyperparameters);
  - `refit` (default `TRUE`: OLS refit of the selected predictors);
  - `selection_control` (`ebic_gamma` 0.5, `stab_B` 50, `stab_fraction` 0.5, `stab_pi` 0.75, `stab_pfer` 1).

  They are saved in `input_data_nodes.rda`.
- **`lambda_prec_type` deprecated:** a warning is given and `"min"`/`"1se"` are mapped to `selection = "cv_min"`/`"cv_1se"`.
- **Legacy input files:** a node file written by an earlier version (no `selection` in the `.rda`) runs with `cv_min` (or its `lambda_prec_type`) and `refit = FALSE`, i.e. the previous behaviour.
- **New function `GGReg_node_model_selection()`.** One code path for all cases:
  - fits the λ path of each (α, ϑ) combination, with `cv.sparsegl` (seed 2024, shared folds) only when a CV rule needs it;
  - computes BIC/EBIC along each path;
  - selects (k, λ);
  - runs stability selection (seed `5000 + i`) when requested;
  - refits by OLS.

  It returns the coefficients, σ² = RSS/(n − df), the penalised coefficients and the selected parameters.
- **Rules:**
  - `cv_min` / `cv_1se`: λ by CV within each combination, combination by BIC;
  - `bic` / `ebic`: joint minimum over (α, ϑ, λ);
  - `stability_cv` / `stability_bic`: Meinshausen & Bühlmann; 50 half-subsamples, π = 0.75, path region with average selection ≤ √(PFER(2π − 1)P).
- **Without covariates:** the lasso is now fitted with `sparsegl` (one group per predictor, no intercept) instead of `cv.glmnet` with intercept.
- **`tune_hyperparams = FALSE`:** the covariate groups get group penalty factor 1, as in the tuned case. Before, they got √(p − 1).
- **Node output:** gains `Penalised_reg_matrix_row`. `collect_node_results()`'s `optimal_params` gains `n_selected`, `selection` and `refit`.

## 5. Template: relaxed BIC / EBIC, new default (2026-11-06)

- **Problem.** BIC and EBIC computed on the RSS of the penalised fit select ever smaller λ as n grows, because the lasso shrinkage inflates the RSS more and more. The selected λ scales like 1/n, so false positives increase with n, and the OLS refit keeps them.
- **New rules** in `GGReg_node_model_selection()`:
  - `bic_relaxed` and `ebic_relaxed`: BIC/EBIC computed on the RSS of the OLS refit of each distinct support along the path;
  - RSS = y'y − b'X_s'y with b = (X_s'X_s)⁻¹X_s'y, from the Gram matrix X'X computed once per node, with a QR fallback if X_s'X_s is singular;
  - identical to `lm.fit` to 1e-14;
  - extra cost about 1 s per node at p = 100.
- **Default `selection = "ebic_relaxed"`** in all entry points (was `"bic"`).
- **`stability_bic`** now takes its hyperparameters from `bic_relaxed`, which adds under 2% to the cost.
- **`bic` / `ebic`** on the penalised RSS are kept for comparison.
- **Checks:**
  - **Against the selection test:** `run_task.R` with `ebic_relaxed` / `bic_relaxed` reproduces its `ebic_relaxed_current` / `bic_relaxed_current` variants exactly (see "Verification" below).
  - **Against the package:** same λ, α, ϑ and network supports for `ebic_relaxed`, `bic_relaxed`, `bic` and `stability_bic` (`prcggm_example`, `mean_estimation = FALSE`).

## 5b. Template: default hyperparameter grids (2026-11-06)

- `asparse_grid = c(0.75, 0.5, 0.9, 0.99)` and `weight_grid = c(1.0, 0.8, 1.1, 1.3)` in all entry points.
  - They were `c(0.75, 0.5, 0.9, 0.95)` and `c(1.0, 0.8, 1.1, 1.3, 1.5)`.
  - These are the values evaluated in the simulations, as in the package.
- The first values are unchanged, so `tune_hyperparams = FALSE` is unaffected.
- `random_hyper_search` stays `FALSE`: the full 16 combinations. The paper's simulations use 50% random search (8 combinations per node), which the paper should state.
- The extended grid of the selection test, α {0.5, 0.75, 0.9, 0.95, 0.99, 1} × ϑ {0.4, 0.6, 0.8, 1, 1.3}, is used only in `Simulation_analysis_v2/selection_test/`.

## 6. `Computational_templates/README.md`

- Method summary: the table of selection rules and the refit.
- Parameter table: `selection` (default `"ebic_relaxed"`), `refit`, `selection_control`, `K`, `lambda_prec`.
- `optimal_params` description updated.

## 7. `Add_simulation_analysis/` (exploratory and sensitivity simulations; Supplementary Part 1)

**Git (2026-11-06):** the code, `README.md`, `code/` (copy of the original paper scripts loaded by `explore_functions.R`) and `archive_v2_code/` are tracked. `results*/`, `logs/` and the Zhang & Li PDF stay ignored (`.gitignore`).

- **Generator** (`R/explore_functions.R`, `simulate_setting()`):
  - Zhang & Li design: power-law baseline; ER covariate networks with probability 0.01; common diagonal (scheme B);
  - x1, x2 ~ U[0, 1], x3 binary; two of the three covariates carry an effect;
  - noisy prior from the support of T₀ = Δ₀ + Σ_continuous x̄ₕΔₕ;
  - truth on the scale of the fit.
- **Runs:**
  - `results/` (first, superseded);
  - `results_v2/` (20 settings);
  - `results_v3/` (100 replications × n ∈ {1k, 5k, 10k, 20k}, p = 100, q = 3) in three passes:
    1. all CV-based rules (job 45365482);
    2. stability at the BIC hyperparameters (job 45378868);
    3. relaxed rules: `--mode relaxed`, `slurm/run_settings_v3_relaxed.sbatch` — ready, to be submitted.
- **`analyse_results_v3.R`:** per-component tables (baseline, binary effect, continuous effect, null) and overall pooled metrics (F1, FPR, Pearson), plus the cost per rule.
- **`README.md`:**
  - hyperparameter tuning and selection rules (reference for the Methods and Supplementary);
  - findings log §1–8.

## 8. `Simulation_analysis_v2/` (new; main text and Supplementary Part 2)

- **Design:**
  - n ∈ {1k, 5k, 10k, 20k} × p ∈ {50, 100, 250, 500} × q ∈ {1, 3} × prior ∈ {perfect, noisy, none}, 100 replications;
  - each dataset is fitted under all three priors;
  - covariate edge probability 1/(p − 1);
  - q = 1: x1 binary; q = 3: x1 continuous, x2 continuous null, x3 binary.
- **Files:**
  - `R/sim_design.R` (generator, prior, evaluation);
  - `run_task.R` (one dataset × three priors via the template; default `selection = "ebic_relaxed"` since 2026-11-06);
  - `analyse_results.R`;
  - `slurm/` (settings, timing test, task and analysis sbatch);
  - copy of the template.
- **`selection_test/`:** rule × grid comparison.
  - 5 replications, with stability selection: `run_selection_test.sbatch`.
  - Part 2 at 100 replications: `run_part2.sbatch`.
  - Analysis scripts: `analyse_selection_test.R`, `overall_metrics.R`.
- **Status:**
  - p = 50 / 100 were run with `selection = "bic"` (penalised RSS); those results were moved to `results_bic_penalised/` (2026-11-06), and the blocks are to be rerun with `ebic_relaxed` into `results/`.
  - p = 250 / 500 are not yet run.
- **Git (2026-11-06):** code and `README.md` tracked; `results*/` (incl. `results_bic_penalised/`, `results_timing/`), `logs/` and `selection_test/results*/` ignored (`.gitignore`).

---

## Verification

| Check | Result |
|---|---|
| Template (step 4) vs package, BIC/EBIC with and without refit | identical (2026-11-04) |
| Template relaxed rules vs package (`ebic_relaxed`, `bic_relaxed`, `bic`, `stability_bic`) | same λ, α, ϑ and network supports (2026-11-06) |
| Template relaxed rules vs selection test (n = 1k, p = 100, q = 3, rep 1, three priors, OR/AND) | identical F1 and Pearson for `ebic_relaxed` and `bic_relaxed` (24/24 rows each) |
| v3 relaxed pass vs first pass (BIC/EBIC rows) | identical (local, n = 1k, rep 1) |
| Selection test `bic_current` vs main study | identical in all 20 cluster tasks |
| Relaxed RSS (Gram) vs `lm.fit` | identical to 1e-14 |

## Decisions (2026-11-06)

1. Add `bic_relaxed` and `ebic_relaxed` to the template and the package.
2. Stability selection with BIC hyperparameters uses the relaxed BIC.
3. Keep the hyperparameter grid α {0.5, 0.75, 0.9, 0.99} × ϑ {0.8, 1, 1.1, 1.3}, with 50% random search, in the simulations.
4. `ebic_relaxed` is the default and the main rule of the paper; OR and AND are both reported.

## Open points

- Commit plan: branch name. Settled: `Add_simulation_analysis/` and `Simulation_analysis_v2/` code tracked, outputs ignored; `.DS_Store` not committed.
- Do not push README/installation changes before CRAN accepts 0.1.0.
- Default grids: aligned with the simulations in the template and the package (5b, decided 2026-11-06).
- Rerun T2D with the new default (plan step 5): to be decided.
- `.DS_Store` is modified; it is not to be committed (decided 2026-11-06).
