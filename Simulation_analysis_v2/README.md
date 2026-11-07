# Simulation study for the paper (v2, full factorial design)

Replaces the simulation in `Simulation_analysis/` (published version), whose generator had the problems documented in `Add_simulation_analysis/README.md` §3. The design choices (selection rule, refit, baseline target) follow the sensitivity analysis in `Add_simulation_analysis/` (§8 of its README, to be reported in the Supplementary material).

Status: **design agreed and timing test done (2026-11-04); ready to submit.** Decisions are marked **[DECIDED]** with the date.

## Objectives
1. **Main text:** show the gain from including prior knowledge (perfect, noisy, none), across sample sizes, numbers of nodes and numbers of covariates.
2. **Supplementary material:** show that the method, as implemented in the code and the package, recovers the baseline and covariate networks and their edge weights, and how this scales with n, p and q (accuracy and computing time).

## Factorial design

| Factor | Levels |
|---|---|
| n (subjects) | 1,000 / 5,000 / 10,000 / 20,000 |
| p (nodes) | 50 / 100 / 250 / 500 |
| q (covariates) | 1 / 3 |
| Prior | perfect / noisy / none |
| Replications | 100 per setting |

4 × 4 × 2 × 3 = 96 settings, 9,600 fitted datasets.

**[DECIDED 2026-11-04] Paired design across priors:** for each (n, p, q, replication), one dataset is generated and fitted under the three priors. The prior comparison is therefore within-dataset, with the same seed, the same data and the same truth, which removes the between-dataset variance from the main comparison. Seed: function of (n, p, q, replication), not of the prior.

## Data generation
As in `Add_simulation_analysis/R/explore_functions.R` (`simulate_setting()`, scheme B):
- **Baseline network:** power-law degree distribution with exponent 2.5 (`igraph::sample_fitness_pl`), p edges, weights ±U[0.35, 0.5].
- **Covariate networks:** Erdős–Rényi, weights ±U[0.175, 0.25], zero diagonal.
- **Diagonal:** common constant d, the smallest that keeps Ω(u) positive definite (minimum eigenvalue 0.1) at all 2^q corners of the covariate box. This guarantees positive definiteness for every subject, and the conditional variances do not depend on u.
- **Baseline target:** network at the mean of the continuous covariates and at the reference level of the binary ones, T₀ = Δ₀ + Σ_{h continuous} x̄ₕΔₕ. The prior is generated from the support of T₀.
- **Noisy prior:** sensitivity 0.8, specificity 0.9, confidence of true edges U[0.7, 1], of false edges U[0.1, 0.6].
- **Truth on the scale of the fit:** standardised proteins and continuous covariates.

**[DECIDED 2026-11-04] Density of the covariate networks: edge probability 1/(p − 1).** With a fixed edge probability of 0.01, the expected number of edges per covariate is 0.01·p(p−1)/2:

| p | Baseline edges | Covariate edges (prob 0.01) | Covariate edges (prob 1/(p−1)) |
|---|---|---|---|
| 50 | 50 | 12 | 25 |
| 100 | 100 | 50 | 50 |
| 250 | 250 | 311 | 125 |
| 500 | 500 | 1,248 | 250 |

With a fixed probability, the covariate networks grow quadratically with p. At p = 500 they would have 2.5× more edges than the baseline network. The common diagonal d would also grow (more negative corner eigenvalues), shrinking all partial correlations, so that p would be confounded with effect size and density. Edge probability 1/(p − 1), i.e. p/2 expected edges per effective covariate. This equals 0.01 at p = 100 (the Zhang & Li setting) and keeps the ratio of covariate to baseline edges fixed at 1/2 for every p.

**[DECIDED 2026-11-04] Covariates.**
- *q = 1:* x1 ~ Bernoulli(0.5) with an effect (case/control, as in the T2D application).
- *q = 3:* x1 ~ U[0,1] with an effect, x2 ~ U[0,1] without effect (null covariate, measures false positives), x3 ~ Bernoulli(0.5) with an effect. The effective covariates are fixed rather than drawn at random: in the sensitivity analysis the covariate type explained most of the variability of the covariate-effect metrics, so each type is reported separately.
- Binary covariates are passed to the template as factors (levels 0/1, dummy `x_1`); continuous ones are standardised by the template.

## Model fitting
- **Code:** the template functions (`ggReg_main_functions.R`, a copy of `Computational_templates/ggReg_main_functions.R` after the step-4 update, see below): `GGReg_cov_single_node_processing()` for each node, run in parallel on the cores of the task, and `collect_node_results()` for the symmetrisation. Data are generated with mean zero and standardised, so the mean-estimation step is skipped (`mean_estimation = FALSE`).
- **[DECIDED 2026-11-04, REVISED 2026-11-06 and 2026-11-07] Selection rule:** `selection = "bic_relaxed"`: BIC on the RSS of the OLS refit of each support along the λ path, jointly with α and ϑ; AND is the main symmetrisation, OR also reported (the 2026-11-04 choice, BIC on the penalised RSS, over-selects as n grows; `ebic_relaxed` was the default on 2026-11-06–07, before the final decision; see "Decisions on the selection rule" below and `Add_simulation_analysis/README.md` §8). `run_task.R` defaults to it since 2026-11-07.
- **Hyperparameter grid:** α ∈ {0.5, 0.75, 0.9, 0.99}, ϑ ∈ {0.8, 1, 1.1, 1.3}, with a 50% random search (8 of the 16 pairs per node, node-specific seed). The same grid is used with and without the prior (without a prior the baseline penalty factor is 1 and ϑ still scales the covariate terms), so the only difference between the three fits is the prior.
- **[DECIDED 2026-11-04] OLS refit** of the selected predictors (`refit = TRUE`); edge weights on the precision scale, δ_jk = −β̂_jk/σ̂²_j.
- **[DECIDED 2026-11-04] Symmetrisation:** both OR and AND are evaluated, from the same node fits (`collect_node_results()` is called twice). Both will be reported; they also inform the guidelines.
- **[DECIDED 2026-11-04] Code used for the fits:** the updated template (step 4, see "Changes to the template"). The package was updated with the same options and gives identical networks for BIC (see "Changes to the package"), which is the evidence for the Supplementary material.

## Evaluation
For each prior, symmetrisation and component (baseline, each covariate):
- density, TPR, FPR, precision, F1;
- Pearson and Spearman correlation with the truth on the true-positive edges;
- Frobenius error (off-diagonal);
- computing time per node.

Covariate components are reported by type (binary, continuous, null).

## Computing cost (timing test, job 45381082, 2026-11-04)
Nodes 1–16 (the hubs of the fitness model, so slightly conservative) of five tasks, q = 3, all three priors, 16 cores (cnode13). No errors.

| Task | Data generation | Wall time per prior (16 nodes) | Mean node time | Peak memory (16 workers) | Elapsed |
|---|---|---|---|---|---|
| p = 500, n = 20k | 13.3 s | 205–209 s | 190–196 s | 32.3 GB | 10:51 |
| p = 500, n = 1k | 1.5 s | 72–76 s | 61–64 s | 5.7 GB | 3:49 |
| p = 250, n = 20k | 3.9 s | 90–95 s | 82–88 s | 20.4 GB | 4:49 |
| p = 250, n = 1k | 0.4 s | 84–89 s | 72–76 s | 4.4 GB | 4:23 |
| p = 100, n = 20k | 1.4 s | 34–35 s | 30–31 s | 9.8 GB | 1:49 |

- The cost is the same for the three priors.
- Each node also has 10–14 s of overhead: the template reloads the input file and rebuilds the design matrix for every node.
- Large n scales as expected from the sensitivity analysis (cost roughly linear in the number of predictors P). Small n with many predictors (n = 1k, P ≈ 1,000–2,000) is 2–4× slower than that extrapolation, because the λ path is harder to fit when P approaches or exceeds n.
- Largest task (p = 500, n = 20k, q = 3): 500 nodes × 3 priors × ~200 s ≈ 5 h on 16 cores.
- Total for 100 replications: about **43,000 core-hours**, of which about 70% is p = 500. With 50 replications at p = 500: about 28,000 core-hours.

Resources per block (with margin):

| Block | Tasks | OFFSET | `--time` | `--mem-per-cpu` |
|---|---|---|---|---|
| p = 50 | 2401–3200 | 2400 | 00:30:00 | 2G |
| p = 100 | 1601–2400 | 1600 | 01:00:00 | 3G |
| p = 250 | 801–1600 | 800 | 03:00:00 | 3G |
| p = 500 | 1–800 | 0 | 10:00:00 | 4G |

**[DECIDED 2026-11-04] Replications at p = 500:** 100, as for the other p. All four blocks submitted together with at most 50 concurrent tasks per block (`%50`).

## Changes to the template (`Computational_templates/ggReg_main_functions.R`, step 4, 2026-11-04)
- New arguments of `GGReg_full_estimation()` and of the estimation functions:
  - `selection`: one of `"bic"` (default), `"ebic"`, `"cv_min"`, `"cv_1se"`, `"stability_bic"`, `"stability_cv"`;
  - `refit = TRUE`: OLS refit of the selected predictors;
  - `selection_control`: `ebic_gamma` (0.5), `stab_B` (50), `stab_fraction` (0.5), `stab_pi` (0.75), `stab_pfer` (1).
- New function `GGReg_node_model_selection()`, used by `GGReg_cov_single_node_processing()` for all cases (with/without covariates, with/without prior, tuned or fixed hyperparameters, fixed `lambda_prec`).
- `lambda_prec_type` is deprecated: `"min"`/`"1se"` are mapped to `selection = "cv_min"`/`"cv_1se"` with a warning. Input files written by the previous version (no `selection`) are read as `"cv_min"` without refit.
- `optimal_params` gains `n_selected`, `selection` and `refit`. The node results also store the penalised coefficients (`Penalised_reg_matrix_row`).
- Behaviour changes not related to the selection rule:
  - Without covariates and without a prior, the lasso is now fitted with `sparsegl` (no intercept, as in all other cases) instead of `cv.glmnet` (with intercept).
  - With `tune_hyperparams = FALSE`, the covariate groups have penalty factor 1 as in the tuned case; before, they had sqrt(p − 1) in this case only.
- Checks (local, p = 12, n = 600, two covariates, prior):
  - `selection = "cv_min", refit = FALSE` reproduces the previous template exactly (networks, weights, residual variances);
  - all six rules run;
  - the refit keeps the support and changes only the weights;
  - the no-prior and no-covariate cases run;
  - the deprecated argument still works, with a warning.
- The template README documents the new options.

## Changes to the package (`Rpackages/prCGGM`, branch `dev-selection-rules`, version 0.1.0.9000, 2026-11-04)
Developed on a separate branch, so that `main` stays at the version submitted to CRAN (0.1.0) until it is accepted; to be released as 0.2.0. Not committed yet.
- `prcggm()` and `prcggm_slurm_prepare()`: new arguments `selection` (same six rules and default `"bic"` as the template), `refit = TRUE` and `selection_control` (same settings); `lambda_prec_type` deprecated (mapped to `cv_min`/`cv_1se` with a warning). New element `coefficients_penalised`; `tuning` gains `df`; `print()`/`summary()` report the rule.
- Stability selection: the subsamples are drawn in `prepare_fit()` before any fit and shared by all nodes, so fits stay reproducible with `set.seed()`, also in parallel. (The template uses node-specific seeds instead, so the stability and CV rules give different random draws in the two implementations; BIC and EBIC are deterministic.)
- Without covariates and prior: `sparsegl` without intercept instead of `cv.glmnet` (as in the template).
- `simulate_prcggm_data()` rewritten to follow this design: fitness-model baseline with power-law exponent (`exponent`, default 2.5; `power` deprecated) and p edges; covariate networks with edge probability `1/(p − 1)` and weights `effect_range` (±U[0.175, 0.25]); common diagonal from the corners (positive definite for every subject, no subject-specific correction); continuous covariates U[0, 1]; true baseline T₀ and prior from T₀; continuous effects per SD; returns `diagonal`. Every simulated covariate has an effect.
- Tests added (default rule and refit, refit keeps the support and is a least-squares fit, all rules, EBIC sparser than BIC, reproducibility of stability selection, deprecated argument, invalid settings, simulation: prior from T₀, positive definiteness at the corners). `devtools::test()`: all pass. `R CMD check --as-cran`: 0 errors, 0 warnings, 3 notes (development version number, no internet time check, no pandoc for NEWS).
- Vignette, NEWS (0.2.0 entry with the breaking changes), DESCRIPTION and WORDLIST updated.
- **Template–package equivalence** (n = 2000, p = 30, one binary and one continuous covariate, noisy prior, full 4 × 4 grid): identical coefficients and selected (α, ϑ), precision matrices equal to 2 × 10⁻¹⁶, for `"bic"` and `"ebic"`, with and without refit. This is the check of the Supplementary material; it can be repeated on a few settings of the full design.

## Files

| File | Content |
|---|---|
| `ggReg_main_functions.R` | Copy of the updated template |
| `R/sim_design.R` | `dataset_seed()`, `covariate_design()`, `simulate_dataset()`, `make_prior()`, `evaluate_networks()` |
| `run_task.R` | One task: one dataset (n, p, q, rep), fitted under the three priors, evaluated with OR and AND |
| `slurm/settings.csv` | The 3,200 tasks (4 n × 4 p × 2 q × 100 replications), ordered by p (500, 250, 100, 50), then n (largest first) |
| `slurm/settings_timing.csv` | The 5 timing tasks |
| `slurm/timing_test.sbatch` | Timing test: nodes 1–16 of each timing task, 16 cores |
| `slurm/run_tasks.sbatch` | Array job of the full study (index + `OFFSET`); resources set per block of p on the command line |

Outputs per task in `results/n<n>_p<p>_q<q>/rep<rep>/`: `evaluation.rds` (metrics per prior, symmetrisation and component), `nodes.rds` (per prior and node: computing time, selected λ, α, ϑ, number of selected predictors) and `config.rds` (seed, diagonal, edge counts, data-generation and fitting times). Data and node files are not stored: each dataset is regenerated from its seed (`dataset_seed()`); the noisy prior uses seed + 7. Input and node files are written to the local temporary directory of the compute node and deleted at the end of the task.

Local test (2026-11-04; n = 1000, p = 50, q = 3 and 1, 6 cores): runs end to end in under 10 s per prior. Baseline F1 (OR): perfect 0.98, noisy 0.88, none 0.72 (q = 3); Pearson on the true-positive edges ≈ 0.99.

## Running on the cluster

From this folder:

```bash
mkdir -p logs
# 1. Timing test: nodes 1-16 of p = 500 and 250 (n = 1000, 20000) and p = 100 (n = 20000), q = 3
sbatch --array=1-5 slurm/timing_test.sbatch
sacct -j <jobid> --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,State
grep -h "Task\|Data generated\|Prior" logs/timing_<jobid>_*.out
# 2. Full study, one block per p (resources from the timing test, see "Computing cost")
sbatch --array=1-800%50 --export=ALL,OFFSET=2400 --time=00:30:00 --mem-per-cpu=2G slurm/run_tasks.sbatch   # p = 50
sbatch --array=1-800%50 --export=ALL,OFFSET=1600 --time=01:00:00 --mem-per-cpu=3G slurm/run_tasks.sbatch   # p = 100
sbatch --array=1-800%50 --export=ALL,OFFSET=800  --time=03:00:00 --mem-per-cpu=3G slurm/run_tasks.sbatch   # p = 250
sbatch --array=1-800%50 --export=ALL,OFFSET=0    --time=10:00:00 --mem-per-cpu=4G slurm/run_tasks.sbatch   # p = 500
```

Note on the timing test: in `sample_fitness_pl` the low node indices have the largest fitness, so nodes 1–16 are the hubs. Their regressions are not more expensive per λ, but they select more predictors, so the timing is slightly conservative.

## Log
- 2026-11-04: design draft; decisions agreed (paired priors, edge probability 1/(p − 1), fixed effective covariates, BIC + refit, OR and AND, template first); template updated (step 4) and checked; simulation code, timing test and submission scripts written and tested locally.
- 2026-11-04: package updated on branch `dev-selection-rules` (selection rules, refit, simulator); tests and R CMD check pass; template–package equivalence verified for BIC/EBIC.
- 2026-11-04: timing test (job 45381082) completed without errors; resources per block set; total cost about 43,000 core-hours (open: 50 or 100 replications at p = 500).
- 2026-11-04: full study submitted, four blocks (p = 50, 100, 250, 500), 100 replications each, `%50` per block.
- 2026-11-05: p = 50 and p = 100 blocks submitted with `%20` (not `%50`) and completed: 1,600 tasks, no errors in the logs, 100 replications in all 16 settings (`results/summary/completeness.csv`, analysis job 45426276). p = 250 and p = 500 not yet submitted.
- 2026-11-05: **issue found: BIC over-selects as n grows.** F1 falls and FPR rises with n (e.g. p = 50, q = 1, Perfect prior, baseline F1 0.95 → 0.87 from n = 1k to 20k) while TPR is already 1. Cause: BIC is computed on the RSS of the penalised fit, which the lasso shrinkage inflates by an amount growing with n; BIC then prefers ever smaller lambda (selected lambda ∝ 1/n, `results/summary/tuning.csv`) and admits noise edges, which the OLS refit keeps. Local check (one dataset per n, p = 50, q = 1, no prior): BIC on the refitted ("relaxed") RSS gives baseline F1 0.92 / 0.97 / 0.99 at n = 1k / 5k / 20k, against 0.85 / 0.79 / 0.65. Same pattern in the v3 sensitivity runs (Add_simulation_analysis README §8). Plan: add `bic_relaxed` and `ebic_relaxed` to the template and the package after the selection test below; p = 50 / 100 to be rerun.
- 2026-11-05: selection test written (`selection_test/`, see below); local check at n = 1000, p = 50, q = 3: the current BIC variant reproduces the main study exactly.
- 2026-11-07: Part 1 (Add_simulation_analysis third pass, job 45436972) and Part 2 (job 45436977) completed without errors. Both confirm the relaxed criteria, the guideline (EBIC relaxed + OR best at small n, BIC relaxed + AND best for n ≥ 5000) and the current grid.
- 2026-11-07: **final decision: default `bic_relaxed` + AND** (best for n ≥ 5000 and for binary effects, the case-study setting). Template, package (`v.0.2`) and `run_task.R` defaults changed; template copy in this folder refreshed. Main study p = 50 / 100 to be rerun into `results/`, then p = 250 / 500.

## Selection test (2026-11-05)
Question: (i) how much does BIC on the relaxed RSS improve over BIC on the penalised RSS, and at what cost; (ii) does an extended hyperparameter grid help, given that the selections pile up at alpha = 0.99 and theta = 0.8; (iii) does stability selection gain from taking its hyperparameters from the relaxed BIC instead of the penalised one, or is the extra cost not worth it.

Design: p = 100, q = 3, n ∈ {1k, 5k, 10k, 20k}, replications 1–5 (the same datasets as the main study), three priors. For each node the lambda path of every combination is fitted once and all variants are applied to the same paths; all report the OLS refit on the selected support:

| Variant | Lambda and combination | Grid |
|---|---|---|
| `bic_current` | BIC on penalised RSS (rule of the main study) | current |
| `bic_extended` | BIC on penalised RSS | extended |
| `bic_relaxed_current` / `_extended` | BIC on relaxed RSS | current / extended |
| `ebic_relaxed_current` / `_extended` | EBIC (gamma 0.5) on relaxed RSS | current / extended |
| `stab_bic` | stability selection at the combination of `bic_current` | current |
| `stab_bic_relaxed` | stability selection at the combination of `bic_relaxed_current` | current |

- Current grid: alpha {0.5, 0.75, 0.9, 0.99} × theta {0.8, 1, 1.1, 1.3}; extended grid: alpha {0.5, 0.75, 0.9, 0.95, 0.99, 1} × theta {0.4, 0.6, 0.8, 1, 1.3}. Both with 50% random search, sampled as in the template (seed 234 + 10 i). alpha = 1 is the lasso without the group penalty.
- Relaxed RSS: for each distinct support on the path, RSS = y'y − b'X_s'y with b = (X_s'X_s)⁻¹X_s'y, from the Gram matrix X'X computed once per node (identical to `lm.fit` to 1e-14; QR fallback if X_s'X_s is singular). Locally (n = 1k, p = 50) it adds 0.2 s per node to the 0.5 s of the paths.
- Timing per variant = paths of its combinations (+ Gram matrix and relaxed RSS) (+ stability subsampling).

Files: `selection_test/run_selection_test.R` (one task = one dataset, three priors), `selection_test/run_selection_test.sbatch` (array 1–20), `selection_test/analyse_selection_test.R` (writes `selection_test/results/summary/`: `metrics.csv`, `F1_OR.csv`, `F1_AND.csv`, `time.csv`, `choices.csv`).

Run from the Simulation_analysis_v2 folder:
```
sbatch --array=1-20 selection_test/run_selection_test.sbatch
sbatch --mem=8G --time=00:30:00 --output=logs/seltest_analyse_%j.out \
  --wrap="source /center/healthds/singularity_functions; Rscript --vanilla selection_test/analyse_selection_test.R"
```
### Results (2026-11-06; 20 tasks, 5 replications per n; summary run locally because the cluster rejected the summary job)
All 20 tasks: `bic_current` reproduces the main study exactly. Tables in `selection_test/results/summary/`.

Seconds per node (mean over nodes, priors, replications; p = 100, q = 3, 16 cores):

| Variant | n = 1k | 5k | 10k | 20k |
|---|---|---|---|---|
| BIC, current grid | 6.0 | 11.5 | 17.1 | 26.3 |
| BIC relaxed, current grid | 7.1 | 12.7 | 18.3 | 27.3 |
| BIC, extended grid | 11.6 | 22.1 | 32.7 | 50.5 |
| BIC relaxed, extended grid | 13.7 | 24.4 | 34.9 | 52.3 |
| Stability at BIC hyperparameters | 63.9 | 67.0 | 96.6 | 128.0 |
| Stability at BIC relaxed hyperparameters | 65.5 | 67.9 | 97.4 | 128.2 |

(EBIC relaxed costs the same as BIC relaxed.)

F1, OR, n = 1k / 5k / 10k / 20k (selected rows; full tables in `metrics.csv`):

| Component, prior | BIC (current) | BIC relaxed | EBIC relaxed | Stability (BIC) | Stability (BIC relaxed) |
|---|---|---|---|---|---|
| Baseline, Perfect | 0.99 / 0.97 / 0.96 / 0.92 | 0.85 / 0.92 / 0.95 / 0.96 | 0.98 / 0.99 / 0.99 / 0.98 | 0.99 / 0.99 / 0.98 / 0.97 | 0.99 / 0.99 / 0.97 / 0.96 |
| Baseline, Noisy | 0.87 / 0.87 / 0.85 / 0.80 | 0.81 / 0.91 / 0.92 / 0.93 | 0.93 / 0.96 / 0.96 / 0.95 | 0.91 / 0.93 / 0.93 / 0.91 | 0.93 / 0.95 / 0.93 / 0.92 |
| Baseline, None | 0.74 / 0.80 / 0.79 / 0.77 | 0.74 / 0.91 / 0.92 / 0.93 | 0.82 / 0.96 / 0.96 / 0.96 | 0.80 / 0.93 / 0.92 / 0.91 | 0.80 / 0.96 / 0.95 / 0.93 |
| Binary, None | 0.48 / 0.80 / 0.71 / 0.61 | 0.54 / 0.84 / 0.87 / 0.87 | 0.65 / 0.94 / 0.94 / 0.91 | 0.46 / 0.93 / 0.92 / 0.85 | 0.53 / 0.91 / 0.91 / 0.85 |
| Continuous, None | 0.04 / 0.59 / 0.81 / 0.80 | 0.14 / 0.63 / 0.85 / 0.96 | 0.08 / 0.58 / 0.87 / 0.98 | 0.02 / 0.39 / 0.81 / 0.97 | 0.01 / 0.49 / 0.83 / 0.97 |

Findings:
1. The relaxed criterion removes the decline with n for all priors and components, at +1.0–1.2 s per node (+4–18%).
2. BIC relaxed over-selects at small n relative to P (P = 396): Perfect-prior baseline OR 0.85 at n = 1k against 0.99. EBIC relaxed corrects this and is the best or close to the best OR rule throughout. Under AND, BIC relaxed is better for the continuous effects (e.g. 0.59 against 0.37 at 5k, Perfect prior), EBIC relaxed slightly too sparse.
3. Extended grid: about twice the cost, gains ≤ 0.03 in most cells (largest: noisy-prior baseline OR, BIC relaxed 0.93 → 0.96 at 20k), sometimes lower for the binary effect. With the relaxed criterion the theta choices no longer pile at 0.8 (current grid: 0.8 / 1 / 1.1 / 1.3 = 31 / 31 / 30 / 8%); alpha = 0.99 remains the modal choice (~50%), and alpha = 1 is chosen in ~28–34% of nodes when available.
4. Stability at the relaxed-BIC hyperparameters costs the same as at the penalised-BIC ones (< 2% extra, the subsampling dominates), with small gains (≤ 0.03 baseline, up to +0.10 continuous at 5k).
5. EBIC relaxed matches or beats stability selection at about one fifth of the cost.
Caveat: 5 replications, p = 100, q = 3 only.

Overall metrics (all components of a fit pooled: F1 and FPR from the summed TP / FP / FN / TN, Pearson as the TP-weighted mean of the per-component correlations): `selection_test/overall_metrics.R` → `selection_test/results/summary/overall_{F1,FPR,Pearson}.csv`. EBIC relaxed is best or tied in 23 of 24 OR/AND × prior × n cells; it is the only rule whose FPR stays flat (~1/1000) as n grows; Pearson is 0.98–0.998 for all rules.

## Decisions on the selection rule [DECIDED 2026-11-06]
1. Add `bic_relaxed` and `ebic_relaxed` to the template and the package.
2. Stability selection with BIC hyperparameters takes them from the relaxed BIC (consistency; extra cost < 2%).
3. Keep the current grid (alpha {0.5, 0.75, 0.9, 0.99} × theta {0.8, 1, 1.1, 1.3}, 50% random search).
4. ~~`ebic_relaxed` is the default and the main rule of the paper~~ **[REVISED 2026-11-07, after Parts 1 and 2]: `bic_relaxed` with AND symmetrisation is the default (template, package) and the main rule of the paper; OR also reported.** It improves further as n grows, is best for n ≥ 5000 and for binary covariate effects, the focus of the case study (large n, binary covariates). Guideline: BIC relaxed with AND (default, n ≥ 5000), EBIC relaxed with OR for small samples (n ≈ 1000), where it is best (Part 1 overall F1 0.73 vs 0.71). Table and considerations in `Add_simulation_analysis/README.md` §8, "Final decision".
Then: rerun p = 50 / 100 with the new rule and submit p = 250 / 500. The p = 50 / 100 results fitted with BIC on the penalised RSS were moved to `results_bic_penalised/` (2026-11-06); the reruns write to `results/` (`run_task.R` defaults to `bic_relaxed` since 2026-11-07; OR and AND are both evaluated for every fit).

## Structure of the simulation results in the paper [DECIDED 2026-11-06]
- **Supplementary Part 1 — selection rules** (`Add_simulation_analysis/results_v3`, 100 replications, p = 100, q = 3, noisy prior, n ∈ {1k, 5k, 10k, 20k}). Rules (all with OLS refit): CV lambda.min, CV lambda.1se, BIC and EBIC on the penalised RSS, BIC relaxed, EBIC relaxed, stability selection at the CV / BIC / BIC-relaxed hyperparameters, lambda.min + pcor threshold, lambda.min + backbone. Logic: (i) overall F1, FPR, Pearson pooled over all components, by n and OR/AND → choose the symmetrisation; (ii) per component (baseline, binary effect, continuous effect, null): density, TPR, FPR, precision, F1, Pearson/Spearman on true positives, Frobenius error, for the chosen symmetrisation; (iii) computing cost per rule; (iv) guidelines: BIC relaxed with AND (default, n ≥ 5000, best for binary effects), EBIC relaxed with OR for small n (confirmed 2026-11-07).
- **Supplementary Part 2 — robustness across priors and grid** (`selection_test/results_100`, 100 replications, p = 100, q = 3, the main-study datasets, three priors, n ∈ {1k, 5k, 10k, 20k}). BIC relaxed and EBIC relaxed (BIC on the penalised RSS kept as reference), current vs extended grid. Logic: both behave well across priors; the extended grid does not improve performance and doubles the cost.
- **Main text — impact of the prior** (this study, `results/`): factorial n × p × q × prior with BIC relaxed, AND (main) and OR.

## Part 2 (selection test at 100 replications)
`selection_test/run_selection_test.R --stab 0` (no stability selection), `selection_test/run_part2.sbatch` (array 1–400: tasks 1–100 n = 20k, 101–200 n = 10k, 201–300 n = 5k, 301–400 n = 1k). Expected about 20 min per n = 20k task (16 cores). Checked locally (n = 1k, p = 50, rep 2): the current BIC reproduces the main study.
```
sbatch --array=1-400%20 selection_test/run_part2.sbatch
# when finished
sbatch --mem=8G --time=01:00:00 --output=logs/part2_analyse_%j.out \
  --wrap="source /center/healthds/singularity_functions; Rscript --vanilla selection_test/analyse_selection_test.R selection_test/results_100; Rscript --vanilla selection_test/overall_metrics.R selection_test/results_100"
```
### Results (job 45436977, 2026-11-07; 400 tasks, 100 replications per n, summary run locally → `selection_test/results_100/summary/`)
Overall F1 (all components pooled), n = 1k / 5k / 10k / 20k, current grid (extended grid within ±0.01 for BIC/EBIC relaxed everywhere; for penalised BIC up to +0.03 at n ≥ 10k, still well below the relaxed criteria):

| Prior | Rule | OR | AND |
|---|---|---|---|
| Perfect | BIC (penalised RSS) | 0.79 / 0.89 / 0.91 / 0.87 | 0.75 / 0.89 / 0.93 / 0.97 |
| Perfect | BIC relaxed | 0.58 / 0.81 / 0.88 / 0.93 | 0.76 / 0.91 / 0.96 / 0.99 |
| Perfect | **EBIC relaxed** | **0.79 / 0.92 / 0.96 / 0.98** | 0.77 / 0.91 / 0.96 / 0.99 |
| Noisy | BIC (penalised RSS) | 0.71 / 0.82 / 0.83 / 0.79 | 0.70 / 0.87 / 0.91 / 0.94 |
| Noisy | BIC relaxed | 0.59 / 0.82 / 0.88 / 0.92 | 0.73 / 0.90 / 0.95 / 0.98 |
| Noisy | **EBIC relaxed** | **0.75 / 0.90 / 0.94 / 0.96** | 0.70 / 0.89 / 0.94 / 0.98 |
| None | BIC (penalised RSS) | 0.62 / 0.77 / 0.77 / 0.71 | 0.63 / 0.86 / 0.90 / 0.91 |
| None | BIC relaxed | 0.55 / 0.81 / 0.88 / 0.91 | 0.66 / 0.90 / 0.95 / 0.98 |
| None | **EBIC relaxed** | **0.69 / 0.89 / 0.94 / 0.95** | 0.65 / 0.88 / 0.94 / 0.98 |

Overall FPR × 1000 (OR): EBIC relaxed 0.4–1.4 under all priors and all n; BIC relaxed falls from ~10–12 (n = 1k) to ~2; penalised BIC rises with n (Perfect 0.6 → 3.7, Noisy 1.8 → 6.9, None 3.3 → 10.4). Pearson on true positives 0.99–0.997 for all rules.

Seconds per node (n = 1k → 20k): penalised BIC 6.2 → 27.7 (current) vs 12.0 → 53.1 (extended); BIC/EBIC relaxed 7.4 → 28.8 (current) vs 14.2 → 54.9 (extended). The relaxed criteria add ~1 s per node; the extended grid doubles the cost.

Conclusions:
1. **The extended grid does not help**: F1 and FPR unchanged (±0.01) for the relaxed criteria under every prior, at twice the cost → keep the current grid.
2. **EBIC relaxed with OR is the best rule under all three priors** and every n; its FPR stays flat (~1/1000). The ranking is the same as in Part 1 (Noisy prior there), so the choice does not depend on the prior.
3. **With AND, BIC relaxed is equal or slightly better** (+0.01–0.03, mainly at n = 1k and on covariate effects) → guideline confirmed: EBIC relaxed with OR, BIC relaxed with AND.
4. BIC relaxed with OR over-selects at small n (FPR ~1% at n = 1k), as expected with P = 396 predictors per node.
5. The prior gain survives the new rule: with EBIC relaxed + OR, Perfect > Noisy > None at every n (e.g. n = 1k: 0.79 / 0.75 / 0.69), the gap shrinking with n. Same with BIC relaxed + AND (n = 1k: 0.76 / 0.73 / 0.66; n = 20k: 0.99 / 0.98 / 0.98).
6. **BIC relaxed + AND (final default) vs EBIC relaxed + OR, overall F1 n = 1k / 5k / 10k / 20k:** Perfect 0.76 / 0.91 / 0.96 / 0.99 vs 0.79 / 0.92 / 0.96 / 0.98; Noisy 0.73 / 0.90 / 0.95 / 0.98 vs 0.75 / 0.90 / 0.94 / 0.96; None 0.66 / 0.90 / 0.95 / 0.98 vs 0.69 / 0.89 / 0.94 / 0.95. As in Part 1, EBIC + OR is better at n = 1k under every prior, BIC + AND equal or better from n = 5k on, so the final decision holds under every prior.
