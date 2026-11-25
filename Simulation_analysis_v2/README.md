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
- 2026-11-14: **main study with `bic_relaxed` complete** (p = 50 / 100 / 250 / 500, analysis job 45606506, `results/summary/`). 3,198 / 3,200 tasks: n = 1000, p = 250, q = 3, reps 20 and 61 (array 620 and 661 of the p = 250 block, job 45451372) hit the 3 h limit (~275 s per node). See "Main study results" below.
- 2026-11-16: the two timed-out tasks (n = 1000, p = 250, q = 3, reps 20 and 61) rerun with `bic_relaxed`: main study complete, 100 replications in all 32 settings (analysis job 45671105). EBIC relaxed rerun at n = 1000 complete, 8 cells × 100 replications (analysis job 45671103, `results_ebic_relaxed/summary/`). **EBIC does not fix the P ≳ n cells: the cause is a degeneracy of the relaxed criteria near saturation** (see "EBIC relaxed at n = 1000: results").

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

## Main study results (`bic_relaxed`, 2026-11-14, analysis job 45606506)
Tables and figures in `results/summary/` (`table_<metric>_<symm>.csv`, `prior_gain_<symm>.csv`, `tuning.csv`, `computing_time.csv`). Complete except n = 1000, p = 250, q = 3 (98 replications; two tasks timed out).

F1, AND, mean over 100 replications, n = 1k / 5k / 10k / 20k (Noisy prior):

| p | q | Baseline | Binary effect | Continuous effect |
|---|---|---|---|---|
| 50 | 1 | 0.97 / 0.99 / 0.99 / 0.99 | 0.73 / 0.99 / 1.00 / 1.00 | — |
| 50 | 3 | 0.93 / 0.98 / 0.99 / 0.99 | 0.69 / 0.98 / 0.99 / 0.99 | 0.10 / 0.59 / 0.87 / 0.97 |
| 100 | 1 | 0.95 / 0.98 / 0.98 / 0.99 | 0.66 / 0.98 / 0.99 / 0.99 | — |
| 100 | 3 | 0.90 / 0.98 / 0.99 / 0.99 | 0.63 / 0.97 / 0.99 / 0.99 | 0.07 / 0.50 / 0.81 / 0.97 |
| 250 | 1 | 0.91 / 0.97 / 0.98 / 0.99 | 0.51 / 0.95 / 0.98 / 0.99 | — |
| 250 | 3 | 0.61 / 0.96 / 0.98 / 0.98 | 0.25 / 0.93 / 0.97 / 0.98 | 0.03 / 0.37 / 0.71 / 0.93 |
| 500 | 1 | 0.70 / 0.96 / 0.97 / 0.98 | 0.31 / 0.93 / 0.97 / 0.99 | — |
| 500 | 3 | 0.03 / 0.93 / 0.96 / 0.98 | 0.01 / 0.89 / 0.94 / 0.97 | 0.01 / 0.29 / 0.60 / 0.87 |

Findings:
1. **For n ≥ 5000 BIC relaxed + AND works at every p**: baseline F1 ≥ 0.90 and binary-effect F1 ≥ 0.82 under every prior, improving with n; continuous effects need n ≥ 10k–20k, more so at large p. AND is above OR for the baseline and binary effects at every n ≥ 5k (gap growing with p: e.g. p = 500, q = 3, n = 20k, binary 0.97 vs 0.79, Noisy); for continuous effects OR is better at n = 5k–10k (e.g. p = 50, q = 3, n = 5k: 0.73 vs 0.59), as in Part 1.
2. **BIC relaxed breaks down when the number of predictors per node P = (p − 1)(q + 1) is close to or above n** (n = 1000): P = 1996 (p = 500, q = 3): ~991 predictors selected per node, ~48,000 baseline edges against 748 true, FPR ≈ 0.38, F1 0.03; P ≈ 1000 (p = 500, q = 1; p = 250, q = 3): 3–2.4× too many edges, baseline F1 0.61–0.70; P = 498 (p = 250, q = 1): mild inflation (F1 0.87–0.91). As df → n the refit RSS → 0, so log(RSS) dominates the log(n) penalty. Part 1 / Part 2 (P = 396, n = 1k) did not reach this regime. EBIC is designed for P ≥ n (its penalty grows with log C(P, df)). These fits are also the slowest (n = 1k, p = 500, q = 1: 6.1 h mean, 8.3 h max per task).
3. **Prior gain is in the hard settings and on the baseline**: q = 3, n = 1000, Perfect − None baseline F1 +0.12 to +0.14 (TPR +0.21 to +0.24) for p ≤ 250, Noisy − None +0.07 to +0.08; it falls to ≤ 0.05 at n = 5k and vanishes at n = 20k, where all priors reach F1 ≥ 0.97. With q = 1 the baseline is already easy and the prior adds ≤ 0.03 (except p = 500).
4. **Covariate effects are hardly affected by the prior** (it penalises only the baseline): ±0.02 for p ≤ 100. At p ≥ 250 and n ≥ 5k the **Perfect prior lowers binary-effect F1 by 0.02–0.04** (TPR unchanged, more false positives: p = 500, q = 1, n = 5k, FPR 0.78 vs 0.24 per 1000 for Noisy), while the Noisy prior raises it (+0.02 to +0.06). The Perfect prior gives penalty factor 0 (unpenalised) to every true baseline edge; the selected λ is then lower (median 0.00064 vs 0.00102 Noisy, 0.00069 None) and more covariate terms enter. To be checked; a floor on the penalty factor (e.g. 1 − 0.9 · prior) would avoid fully unpenalised predictors.
5. Pearson on true positives ≥ 0.99 wherever edges are found.

Open: (a) how to handle P ≳ n — rerun the n = 1000 cells with `ebic_relaxed` and/or make the rule depend on P / n; (b) rerun the two timed-out tasks with a longer limit; (c) the Perfect-prior effect on covariate components.

### EBIC relaxed at n = 1000 (decided 2026-11-14)
Rerun of the n = 1000 cells of all p and q with `selection = "ebic_relaxed"` into `results_ebic_relaxed/`, to handle P ≳ n (finding 2) and to set the BIC / EBIC guideline in terms of P / n. `slurm/run_tasks.sbatch` takes `SELECTION` and `OUT`, `slurm/analyse.sbatch` takes `RES` (checked locally: `config.rds` records `ebic_relaxed`). Within each p block the n = 1000 tasks are array indices 601–800.
```
sbatch --array=601-800%50 --export=ALL,OFFSET=2400,SELECTION=ebic_relaxed,OUT=results_ebic_relaxed --time=00:30:00 --mem-per-cpu=2G slurm/run_tasks.sbatch   # p = 50
sbatch --array=601-800%50 --export=ALL,OFFSET=1600,SELECTION=ebic_relaxed,OUT=results_ebic_relaxed --time=01:00:00 --mem-per-cpu=3G slurm/run_tasks.sbatch   # p = 100
sbatch --array=601-800%50 --export=ALL,OFFSET=800,SELECTION=ebic_relaxed,OUT=results_ebic_relaxed  --time=04:00:00 --mem-per-cpu=3G slurm/run_tasks.sbatch   # p = 250
sbatch --array=601-800%50 --export=ALL,OFFSET=0,SELECTION=ebic_relaxed,OUT=results_ebic_relaxed    --time=10:00:00 --mem-per-cpu=4G slurm/run_tasks.sbatch   # p = 500
# when finished
sbatch --export=ALL,RES=results_ebic_relaxed slurm/analyse.sbatch
```

### EBIC relaxed at n = 1000: results (2026-11-16, analysis job 45671103)
F1 at n = 1000, Noisy prior, BIC relaxed → EBIC relaxed (P = (p − 1)(q + 1) predictors per node):

| P (p, q) | AND baseline | AND binary | OR baseline | OR binary | OR continuous |
|---|---|---|---|---|---|
| 98 (50, 1) | 0.97 → 0.99 | 0.73 → 0.63 | 0.86 → 0.94 | 0.73 → 0.79 | — |
| 196 (50, 3) | 0.93 → 0.91 | 0.69 → 0.53 | 0.86 → 0.92 | 0.69 → 0.71 | 0.21 → 0.11 |
| 198 (100, 1) | 0.95 → 0.99 | 0.66 → 0.53 | 0.79 → 0.93 | 0.61 → 0.71 | — |
| 396 (100, 3) | 0.90 → 0.90 | 0.63 → 0.41 | 0.80 → 0.91 | 0.58 → 0.64 | 0.13 → 0.06 |
| 498 (250, 1) | 0.91 → 0.99 | 0.51 → 0.34 | 0.65 → 0.94 | 0.41 → 0.54 | — |
| 996 (250, 3) | 0.61 → 0.61 | 0.26 → 0.11 | 0.13 → 0.13 | 0.04 → 0.02 | 0.01 → 0.01 |
| 998 (500, 1) | 0.70 → 0.83 | 0.31 → 0.15 | 0.39 → 0.70 | 0.20 → 0.27 | — |
| 1996 (500, 3) | 0.03 → 0.03 | 0.01 → 0.01 | 0.02 → 0.02 | 0.01 → 0.01 | 0.00 → 0.00 |

Findings:
1. **Where nothing saturates (P ≤ 498) the Part 1 guideline holds at every P**: EBIC relaxed + OR vs BIC relaxed + AND, binary effect 0.79 / 0.71 / 0.71 / 0.64 / 0.54 vs 0.73 / 0.69 / 0.66 / 0.63 / 0.51, baseline similar (0.91–0.94 vs 0.90–0.97). EBIC + AND is too sparse for covariate effects (binary −0.08 to −0.22); EBIC lowers continuous-effect F1 under both symmetrisations (TPR falls).
2. **EBIC does not fix P ≳ n; the problem is a degeneracy, not the BIC penalty.** Per-node selections are bimodal: either ≤ 10 predictors or saturated (994–996 of 996 predictors, λ ≈ 1e-6, the end of the path). p = 500, q = 3: ~496 of 500 nodes saturated under both rules; p = 250, q = 3: ~43 of 250 nodes (rep 1). With df ≈ P ≈ n the OLS refit nearly interpolates (n − df ≈ 4), RSS → 0 and n log(RSS / n) → −∞, which beats any finite penalty, BIC or EBIC. Checked on node 11 of n = 1000, p = 250, q = 3, rep 1 (debug copy of the template, 2026-11-16): the selected support has df = 996 = P, λ = 1.3e-6, RSS = 0.165 against y'y = 999 (Gram RSS identical to `lm.fit`, so not a numerical artefact), EBIC = −1827 (n log(RSS/n) = −8707 against the penalty log(n)·df = 6880). The support is chosen by the lasso path to fit y, so near saturation the refit RSS is far below its nominal (n − df)σ². The EBIC term 2γ log C(P, df) does not help because it returns to 0 as df → P (C(P, P) = 1). p = 500, q = 1 (P = 998) saturates only in some replications, which EBIC partly reduces.
3. **Proposed fix**: exclude from the relaxed criteria (BIC and EBIC) the supports with too many predictors, e.g. df ≤ n / 2 (new `selection_control` setting), in the template and the package; then rerun the n = 1000 cells with both rules. The true node degrees are ≤ 10, so the cap does not bind in any non-degenerate fit (n ≥ 5000 unaffected: there n − P ≥ 3000).

### Diagnostics on single nodes and on the prior (local, 2026-11-16)
**Stability selection on a saturated node** (n = 1000, p = 250, q = 3, rep 1, no prior; scratchpad script `node_test.R`). TP / FP per node:

| Rule | Node 11 baseline (4 true) | Node 11 covariates (3 true) | Node 1 baseline (5 true) | Node 1 covariates (5 true) | s per node |
|---|---|---|---|---|---|
| BIC relaxed | 4 / 245 | 3 / 744 | 3 / 0 | 0 / 3 | 63–73 |
| EBIC relaxed | 4 / 245 | 3 / 744 | 2 / 0 | 0 / 0 | 68–71 |
| Stability (BIC) | 2 / 0 | 0 / 0 | 2 / 0 | 0 / 0 | 90–100 |

Stability selection avoids the saturation (it only uses the part of the path with small models) but is conservative at n = 1000; cost ~1.3–1.4× BIC relaxed. Cluster pilot submitted: `stability_bic`, n = 1000, p = 250, q = 3, reps 1–10 → `results_stability/`.

**Where the false positives fall under each prior** (n = 5000, p = 250, q = 1, reps 1–2, nodes 1–100, `bic_relaxed`, node-wise supports before symmetrisation; scratchpad script `prior_fp.R`). Totals over 200 node fits:

| Prior | Baseline TP / FN | Baseline FP | of which on (false) prior edges | Binary TP / FN | Binary FP | of which on true baseline edges | median λ | mean predictors |
|---|---|---|---|---|---|---|---|---|
| Perfect | 275 / 0 | 100 | — | 178 / 0 | 123 | 0 | 0.00069 | 3.38 |
| Noisy | 275 / 0 | 69 | 36 (of 4,817 false prior pairs) | 176 / 2 | 56 | 3 | 0.00102 | 2.88 |
| None | 275 / 0 | 100 | — | 177 / 1 | 87 | 3 | 0.00079 | 3.20 |

- True positives are the same under all priors (TPR ≈ 1): at q = 1 and n = 5000 the prior can only change false positives.
- **Hypothesis rejected:** the extra binary-effect false positives of the Perfect prior are *not* on the unpenalised true baseline pairs (0 of 123); they are spread over other pairs.
- The differences follow the selected λ: Noisy selects the largest λ, Perfect the smallest, so Noisy has the fewest false positives everywhere (baseline off-prior 33 vs 100; binary 56 vs 123), at the price of 36 false positives on its own false-prior edges (0.75% of them).
- Likely mechanism (not verified): the prior changes the order in which predictors enter the λ path, and therefore which supports the relaxed BIC compares. With the Noisy prior, partially favoured false pairs enter early without improving the fit, and the minimum lands at a larger λ; with the Perfect prior all true baseline edges are in from the start and the minimum lands at a smaller λ, where more noise terms have entered. The effect is small in absolute terms (≤ 0.6 false positives per node).


**Interpretation for the publication.**
- At q = 1 and n ≥ 5000, the "messy" effect of the prior on baseline F1 (see `results/summary/`) is a false-positive effect, not a recovery effect. TPR is already about 1 without a prior (`Fig_TPR_AND.png`), so the prior can only move the selected λ, and with it the number of false positives.
- That the noisy prior has fewer false positives than the perfect prior is a side effect of how the prior interacts with BIC selection along the path. It is not evidence that a noisy prior is better than a perfect one, and the paper should not present it that way.
- The prior's real benefit is recovering true edges in the hard settings (q = 3, smaller n), where TPR without a prior is well below 1. That is the message for the main text.
- Suggested wording: *"When the edges are already recovered without prior information (q = 1, n ≥ 5000), the prior mainly shifts the λ selected by the relaxed BIC, and the resulting differences in F1 between priors are small (≤ 0.6 false positives per node) and not systematic."*
- Where to put it: as a short remark in the results, with this table in the Supplementary Material if a reviewer asks why the noisy prior can beat the perfect prior at q = 1.

### Stability selection pilot at n = 1000, p = 250, q = 3 (job 45679218, reps 1–10, analysed locally 2026-11-17)
`selection = "stability_bic"` (B = 50 half-samples, π = 0.75, PFER = 1, α and ϑ from the relaxed BIC), compared with BIC and EBIC relaxed on the same 10 datasets (P = 996, P/n ≈ 1). Means over 10 replications; FPR × 1000.

| Component | Prior | BIC relaxed AND: TPR / FPR / F1 | BIC relaxed OR: F1 | Stability AND: TPR / FPR / F1 | Stability OR: TPR / FPR / F1 |
|---|---|---|---|---|---|
| Baseline | Perfect | 0.99 / 20.4 / 0.55 | 0.09 | 0.91 / 0.04 / 0.95 | 0.99 / 0.32 / **0.98** |
| Baseline | Noisy | 0.90 / 20.0 / 0.51 | 0.09 | 0.82 / 0.24 / 0.89 | 0.90 / 1.60 / **0.89** |
| Baseline | None | 0.76 / 21.1 / 0.44 | 0.08 | 0.67 / 0.07 / 0.80 | 0.69 / 0.54 / **0.80** |
| Binary effect | Noisy | 0.45 / 18.7 / 0.15 | 0.02 | 0.09 / 0.02 / 0.16 | 0.22 / 0.22 / **0.34** |
| Continuous effect | Noisy | 0.06 / 18.9 / 0.02 | 0.01 | 0.005 / 0.02 / 0.01 | 0.005 / 0.09 / 0.01 |

EBIC relaxed is close to BIC relaxed (baseline F1 0.44–0.54 AND, 0.07–0.09 OR).

Findings:
1. **Stability selection removes the saturation.** Predictors selected per node: mean 2.4–3.1, max 12, no node above 900. BIC and EBIC relaxed: median 3–10 but 324–342 of 2,500 node fits saturated (996 of 996), mean 137.
2. **Baseline recovered at P ≈ n:** F1 0.80 / 0.89 / 0.98 (None / Noisy / Perfect) with OR, against 0.44–0.55 for BIC relaxed + AND. This is comparable to BIC relaxed in the non-saturated cells at n = 1000 (e.g. p = 100, q = 3: 0.82–0.96). Precision 0.87–0.99.
3. **The prior gain is clear** under stability selection (+0.09 Noisy, +0.18 Perfect vs None, baseline F1, OR), and it comes through TPR, as in the main study.
4. **Covariate effects remain hard:** binary F1 0.34–0.36 with OR (precision 0.73–0.91, TPR 0.22–0.24), better than BIC relaxed (0.15–0.16 AND). Continuous effects are not recovered by any rule at n = 1000 (F1 ≤ 0.03). They need n ≥ 10k even when P/n is small.
5. **With stability selection OR is better than AND** (the rule is already conservative), the opposite of the relaxed BIC. Guideline: when P/n ≥ 0.5 use `stability_bic` with OR.
6. **Cost:** 200 s per node vs 99 s for BIC relaxed (≈ 2×); about 50–60 min per prior per replication on 16 cores.

To feed `paper_figures.R` (FigS_n1000): `sbatch --export=ALL,RES=results_stability slurm/analyse.sbatch`.

### Stability selection in the other P/n ≥ 0.5 cells (planned 2026-11-17, to run when the cluster is free)
- Cells: n = 1000 with p = 500, q = 1 (P = 998) and p = 500, q = 3 (P = 1996). Replications 1–10, the same datasets as the BIC / EBIC relaxed fits. 10 replications rather than 100: the F1 standard deviations across replications are ≤ 0.07, so the standard error of a mean is ≤ 0.02, far below the differences between rules (≥ 0.3 in the pilot). Report as "10 replications (paired)" in S7.
- Array indices in `slurm/settings.csv` (OFFSET = 0, p = 500 block): 601–610 = n 1000, q = 3, reps 1–10; 701–710 = n 1000, q = 1, reps 1–10.
- Expected cost (unknown, hence step 1): BIC relaxed took 3.1 h (q = 3) and 6.1 h (q = 1, max 8.3 h) per task; stability took ≈ 2× BIC relaxed in the pilot.
```
# 1. Timing: 16 nodes, Noisy prior, rep 1 (writes nodes_timing.rds only)
# (2026-11-19: the first attempt used `sbatch --wrap`, which runs under /bin/sh; sourcing
#  singularity_functions then fails at line 43 (`aws-cli': not a valid identifier) and the job
#  exits before Rscript. Use the bash script instead.)
for Q in 1 3; do sbatch --output=logs/stabtime_q${Q}_%j.out slurm/stability_timing.sbatch $Q; done
# expected wall per task ≈ mean node time (in the log) × 500 nodes × 3 priors / 16 cores; set --time to ~1.5× that

# Timing result (2026-11-20; jobs 45930943 q = 3, 45930942 q = 1; 16 nodes, Noisy prior, n = 1000, p = 500):
#   q = 3: mean node time 314.5 s (16 nodes in 373 s) -> 314.5 x 500 x 3 / 16 = 8.2 h per task
#          (BIC relaxed 3.1 h: x2.6)
#   q = 1: mean node time 454.3 s (16 nodes in 502 s) -> 454.3 x 500 x 3 / 16 = 11.8 h per task
#          (BIC relaxed 6.1 h, max 8.3 h: x1.9)
#   --time set to about 1.5-1.7x the expected time to cover variation across replications
#   (BIC relaxed max/mean = 1.36 for q = 1).
# 2. Full run (replace --time with the estimate from step 1)
sbatch --array=601-610%10 --export=ALL,OFFSET=0,SELECTION=stability_bic,OUT=results_stability --time=14:00:00 --mem-per-cpu=4G slurm/run_tasks.sbatch   # p = 500, q = 3 (expected about 8.2 h)
sbatch --array=701-710%10 --export=ALL,OFFSET=0,SELECTION=stability_bic,OUT=results_stability --time=20:00:00 --mem-per-cpu=4G slurm/run_tasks.sbatch   # p = 500, q = 1 (expected about 11.8 h)
# 3. When finished (summarises all three stability cells), then rerun paper_figures.R
sbatch --export=ALL,RES=results_stability slurm/analyse.sbatch
```

### Stability selection in all P/n ≥ 0.5 cells: results (jobs 45975774 q = 3 and 45975775 q = 1, analysed locally 2026-11-22)

**Run.** All 20 tasks completed: 601–610 (p = 500, q = 3) and 701–710 (p = 500, q = 1). No time-outs.
- **Task wall time** (mean / max over replications): q = 3 539 / 623 min (estimate 8.2 h); q = 1 643 / 841 min (estimate 11.8 h).
- **Mean node time:** 342 s (q = 3) and 407 s (q = 1).

The summary was written with `Rscript --vanilla analyse_results.R results_stability`, which covers all three cells, 10 replications each.
`paper_figures/FigS_n1000.{png,pdf}` was regenerated from it; the pilot-only version is kept as `FigS_n1000_v1_pilot.*`.

**p = 500 cells**, means over 10 replications, FPR × 1000. BIC and EBIC are the relaxed rules.

| q (P, P/n) | Component | Prior | BIC AND F1 | EBIC AND F1 | Stab AND: TPR / FPR / F1 | BIC OR F1 | EBIC OR F1 | Stab OR: TPR / FPR / F1 |
|---|---|---|---|---|---|---|---|---|
| 1 (998, 1.0) | Baseline | None | 0.61 | 0.82 | 0.99 / 0.04 / 0.99 | 0.33 | 0.67 | 1.00 / 0.23 / 0.97 |
| | | Noisy | 0.70 | 0.83 | 0.99 / 0.20 / 0.97 | 0.39 | 0.70 | 1.00 / 1.05 / 0.88 |
| | | Perfect | 0.64 | 0.82 | 0.99 / 0.02 / 0.99 | 0.33 | 0.67 | 1.00 / 0.18 / 0.98 |
| | Binary | Noisy | 0.31 | 0.15 | 0.07 / 0.02 / 0.13 | 0.20 | 0.28 | 0.19 / 0.18 / 0.30 |
| 3 (1,996, 2.0) | Baseline | None | 0.03 | 0.03 | 0.64 / 0.07 / 0.78 | 0.02 | 0.02 | 0.68 / 0.49 / 0.77 |
| | | Noisy | 0.03 | 0.03 | 0.77 / 0.31 / 0.85 | 0.02 | 0.02 | 0.89 / 1.80 / 0.82 |
| | | Perfect | 0.03 | 0.03 | 0.85 / 0.03 / 0.92 | 0.02 | 0.02 | 0.98 / 0.44 / 0.96 |
| | Binary | Noisy | 0.01 | 0.01 | 0.01 / 0.00 / 0.01 | 0.01 | 0.01 | 0.06 / 0.04 / 0.12 |
| | Continuous | Noisy | 0.01 | 0.01 | 0.00 / 0.00 / 0.00 | 0.00 | 0.00 | 0.00 / 0.01 / 0.00 |

**Observations:**
1. **Stability selection rescues the baseline when P > n.** At P/n = 2, BIC and EBIC relaxed saturate: FPR ≈ 0.4–0.7 and F1 ≈ 0.03. Stability keeps FPR below 0.002 and reaches F1 0.78–0.96. At P/n = 1, F1 rises from 0.70 / 0.83 (BIC / EBIC AND) to 0.97, with much smaller variability across replications.
2. **Covariate effects remain hard at n = 1000.** Binary effect F1 is about 0.3 at P/n = 1, which is no better than the best of BIC and EBIC, and about 0.12 at P/n = 2. The continuous effect is not recovered by any rule (F1 ≈ 0). Stability is conservative: TPR for the binary effect is 0.06–0.19.
3. **AND vs OR under stability.** The baseline is similar under both rules (AND slightly better with the Noisy prior: 0.97 vs 0.88 and 0.85 vs 0.82). The binary effect is clearly better with OR (0.30 vs 0.13; 0.12 vs 0.01), as in the pilot (p = 250). The guideline "P/n ≥ 0.5 → `stability_bic` with OR" stands, because OR is the only rule that recovers covariate effects. AND may be used if only the baseline matters.
4. **The prior gain depends on P/n.** At P/n = 2 the prior helps (baseline AND F1: None 0.78, Noisy 0.85, Perfect 0.92). At P/n = 1 the baseline is already at the ceiling without a prior (0.99), and the Noisy prior lowers it slightly (0.97), through a few more false positives.
5. **Cost.** Mean task wall time vs BIC relaxed: 9.0 h vs 3.1 h for q = 3 (about 2.9×) and 10.7 h vs 6.1 h for q = 1 (about 1.8×).

**FigS_n1000_prior (added 2026-11-22).** New block in `paper_figures.R`, written to `paper_figures/FigS_n1000_prior.{pdf,png}`. It shows F1 (mean ± SD) at n = 1000 by prior (None / Noisy / Perfect): rows = rule × symmetrisation (BIC relaxed, EBIC relaxed, stability; AND / OR), columns = component.
What it shows:
- **The prior gain is concentrated in the baseline of the q = 3 settings.** For example, P = 396 with BIC AND: None 0.82, Noisy 0.90, Perfect 0.96. At P = 1,996 with stability AND: 0.78 → 0.85 → 0.92.
- **The q = 1 settings are already near the ceiling** without a prior.
- **Covariate effects (binary, continuous) are essentially insensitive to the prior** under every rule (differences within the SD; EBIC AND with the Perfect prior is slightly higher for the binary effect).
- **With stability selection the Noisy prior can fall below None:** P = 998, both rules, through a few extra false positives.

## Reproducibility experiment (designed 2026-11-23)
**Question.** How variable is the estimated network across independent datasets drawn from the *same* true network, for BIC relaxed, EBIC relaxed and stability selection, with and without the prior? The main study cannot answer this: each replication draws its own graph and only metrics were stored.

**Design.**
- p = 100, q = 3 (continuous effect x1, null covariate x2, binary effect x3); n ∈ {1k, 5k, 10k}; R = 50 independent datasets per n.
- **Fixed truth:** the baseline and covariate networks are drawn from a graph seed `repro_graph_seed(p, q)` = 5e7 + 10p + q, the same for every n and dataset (new argument `graph_seed` of `simulate_dataset()`; with `graph_seed = NULL` the main-study datasets are unchanged, checked on n = 1000, p = 100, q = 3, rep 1). The noisy prior is fixed too (seed graph seed + 7). Datasets differ only in covariates and subjects (seed `repro_data_seed()`, 6e7 + …, disjoint from `dataset_seed()`). True supports: baseline target 162 edges, x1 62, x3 54, x2 0.
- Rules: `bic_relaxed`, `ebic_relaxed`, `stability_bic` (B = 50, π = 0.75, PFER = 1); priors: noisy and none; OR and AND stored. Everything else as in the main study (grid, 50% random hyperparameter search with node-specific seeds, OLS refit). Note: the random hyperparameter subset and the stability subsample seeds are node-specific and therefore the same across datasets; the variability measured is sampling variability of the data.
- The estimated networks are stored (edge lists, `networks.rds`), not only the metrics.

**Metrics** (`analyse_repro.R`, per n × rule × prior × symmetrisation × component):
- pairwise Jaccard similarity of the edge sets over the 1,225 pairs of datasets;
- edge-selection frequency across datasets, for true edges (mean, share selected in ≥ 80% of datasets, never selected) and for non-edges (mean × 1000, number selected in ≥ 20% / ≥ 50% of datasets); instability index mean 4f(1 − f) over edges ever selected;
- number of estimated edges (density): mean, SD, CV; mean F1 / TPR / FPR for reference.
- Null covariate x2: density and recurrent false positives only.

**Files.** `run_repro_task.R`, `slurm/settings_repro.csv` (150 rows: tasks 1–50 n = 10k, 51–100 n = 5k, 101–150 n = 1k), `slurm/run_repro.sbatch`, `analyse_repro.R`. Outputs in `results_reproducibility/n<n>_p100_q3/rep<r>/` and `results_reproducibility/summary/` (`repro_summary.csv`, `repro_jaccard_pairs.csv`, `repro_edge_frequency.csv`, `Fig_repro_jaccard`, `Fig_repro_frequency`, `Fig_repro_density`).

**Cost.** From the main study and the selection test (p = 100, q = 3, 16 cores): BIC / EBIC relaxed ≈ 55 s (1k) to 140 s (10k) per prior for 100 nodes; stability ≈ 65–97 s per node. Expected ≈ 20 min (1k) to 35 min (10k) per task (3 rules × 2 priors), ≈ 1,200 core-hours in total. Local test (4 nodes, n = 1k): 7.5 s per node BIC/EBIC, 44 s stability. Pipeline tested end to end on p = 50, 3 datasets.

**Submission** (from Simulation_analysis_v2):
```
sbatch --array=1-150%30 slurm/run_repro.sbatch
# when finished
Rscript --vanilla analyse_repro.R results_reproducibility
```

### Reproducibility experiment: results (job 46038647, analysed locally 2026-11-23)
**Run.** All 150 tasks completed (50 datasets per n), no errors. Wall time per prior on 16 cores (mean): BIC and EBIC relaxed 1.0 / 1.7 / 2.3 min at n = 1k / 5k / 10k; stability 8.5 / 9.3 / 12.1 min (5–9× more). Outputs: `results_reproducibility/summary/` (`repro_summary.csv` holds every number below; figures `Fig_repro_jaccard`, `Fig_repro_frequency`, `Fig_repro_density`).

Reading key: J = mean pairwise Jaccard over the 1,225 pairs of datasets; stable TP = share of true edges selected in ≥ 80% of the 50 datasets; recurrent FP = non-edges selected in ≥ 20% of datasets. True edges: baseline 162, binary effect (x3) 54, continuous effect (x1) 62, null covariate (x2) 0. Pairs of datasets are not independent, so no tests are reported; the boxplots show the spread.

**Baseline, AND** (J / stable TP / recurrent FP; None → Noisy prior):

| Rule | n = 1k | n = 5k | n = 10k |
|---|---|---|---|
| BIC relaxed | 0.66 → 0.74 / 0.62 → 0.72 / 0 → 3 | 0.84 → 0.93 / 0.88 → 0.97 / 2 → 2 | 0.94 → 0.96 / 1.00 → 1.00 / 2 → 2 |
| EBIC relaxed | 0.89 → 0.88 / 0.62 → 0.69 / 0 → 0 | 0.84 → 0.93 / 0.75 → 0.92 / 1 → 2 | 0.94 → 0.98 / 0.98 → 0.99 / 2 → 2 |
| Stability | 0.93 → 0.88 / 0.62 → 0.72 / 0 → 1 | 0.81 → 0.91 / 0.69 → 0.89 / 0 → 2 | 0.92 → 0.96 / 0.95 → 0.98 / 1 → 2 |

**Baseline, OR:** J 0.43–0.54 (BIC relaxed), 0.76–0.81 (EBIC, stability) at 1k; 0.72–0.90 at 5k–10k. Recurrent FP: BIC relaxed 14–27, EBIC 3–22, stability 1–20 (increasing with n).

**Binary effect (x3), AND:** J 0.07–0.19 at 1k for every rule (stable TP ≤ 0.02); 0.83–0.89 at 5k and 0.93–0.98 at 10k (stable TP 0.94–1.00; no recurrent FP). OR: BIC relaxed J 0.53–0.66 at 5k–10k with 3–10 recurrent FP, EBIC and stability 0.75–0.86 with 0–6.

**Continuous effect (x1):** not reproducible at these n. AND J ≤ 0.004 at 1k, 0.03–0.14 at 5k, 0.17–0.39 at 10k; OR up to 0.40–0.48 at 10k (stable TP 0.10–0.58, highest for BIC relaxed + OR). Recurrent FP 0 for all rules.

**Null covariate (x2):** no recurrent false positives for any rule, n or symmetrisation (J ≤ 0.005). Mean edges per dataset: BIC relaxed AND 1–11, OR 8–43 (largest at 1k); EBIC ≤ 4; stability ≤ 0.5. The false positives of the null component are scattered, not reproducible.

**Density (number of estimated edges, CV across datasets):** baseline CV ≤ 5% for every rule and n (≤ 2.5% at n ≥ 5k with AND); binary effect CV 2–7% at n ≥ 5k, 17–44% at 1k; continuous effect CV 9–67% at n ≥ 5k and > 100% for EBIC / stability at 1k (mostly empty networks).

**Findings.**
1. **The default (BIC relaxed + AND, n ≥ 5000) is reproducible:** baseline J 0.84–0.96, binary-effect J 0.83–0.93, ≥ 88% of true baseline edges and ≥ 98% of true binary-effect edges selected in ≥ 80% of datasets, ≤ 2 recurrent false positives, edge-count CV ≤ 2.5%.
2. **Stability selection is not more reproducible than the relaxed criteria at n ≥ 5000** (baseline AND J 0.81–0.96 vs 0.84–0.98), at 5–9× the cost. Its advantage is at n = 1k: J 0.88–0.93 vs 0.66–0.74 for BIC relaxed. This supports the guideline (relaxed BIC by default, stability / EBIC when n is small relative to P).
3. **High Jaccard ≠ good recovery.** At n = 1k EBIC and stability have J ≈ 0.9 because they reproducibly select the same conservative subset: only 62–72% of true baseline edges are stable, the same as BIC relaxed. Reproducibility must be read with the stable-TP share (or F1).
4. **The prior improves reproducibility** where TPR is not saturated: baseline AND J +0.04–0.10 at n = 5k for all rules, stable TP +0.09–0.20 (5k) and +0.07–0.10 (1k); binary effect J +0.01–0.10 (OR). Exception: at 1k EBIC / stability J drops slightly (0.89 → 0.88, 0.93 → 0.88) because the prior adds true edges that are not yet selected in every dataset (stable TP 0.62 → 0.69–0.72). The prior does not change the continuous-effect or null components.
5. **AND over OR for the relaxed BIC is confirmed** from the reproducibility side: OR adds 10–27 recurrent false baseline edges and lowers J by 0.12–0.23. With OR, EBIC and stability are more reproducible than BIC relaxed.
6. **Covariate effects:** binary effects become reproducible between n = 1k and 5k; continuous effects need n > 10k (consistent with the main study, F1 ≤ 0.7 at 10k).
7. **The null covariate never produces recurring edges,** so spurious covariate-specific edges do not replicate across samples.

Suggested wording (Supplementary, with `Fig_repro_jaccard` and the table above): *"Across 50 independent datasets drawn from the same network, BIC relaxed with AND symmetrisation recovered nearly the same baseline and binary-effect networks at n ≥ 5000 (mean pairwise Jaccard 0.84–0.96; ≥ 88% of true edges selected in at least 80% of datasets; edge-count CV ≤ 2.5%). Stability selection was not more reproducible in this regime, at 5–9 times the cost, but was at n = 1000. The prior increased reproducibility where recovery was incomplete, and spurious edges on a covariate without effect never recurred across datasets."*

## Ablation: prior on all terms (set up 2026-11-24)
**Question.** Is it right to apply the prior to the baseline only? The naive alternative applies the same prior weights to every component. This is the experiment behind the paper's first contribution ("selective prior").

**Design.** The same setting as the selection-rule study: p = 100, q = 3, n ∈ {1k, 5k, 10k, 20k}, replications 1–100, the **same datasets as the main study** (`dataset_seed()`), `bic_relaxed`, OR and AND.
- Arms 1–3 are the main-study fits in `results/`: no prior; noisy prior on the baseline only (prCGGM); perfect prior on the baseline only (prCGGM).
- Arms 4–5 are new, in `results_ablation/`: noisy prior on all terms; perfect prior on all terms.
- **Prior on all terms:** the covariate-specific term h × node k gets penalty factor ϑ (1 − W_jk) instead of ϑ. Template option `prior_on_covariates = TRUE` (`ggReg_main_functions.R`, both copies; default FALSE) and `run_task.R --prior_on_covariates TRUE`.
- **Checks (local, n = 1k, rep 1, nodes 1–4):**
  - with the option off, the selections reproduce the main study exactly (n_selected and λ per node);
  - with the option on, the perfect-prior selections change (nodes 2 and 3), and the no-prior ones do not.
- **What we expected (before the run; see the results below, which differ).** The prior is built from T₀ = Δ₀ + mean(x1)·B1, so it also contains the continuous-effect edges. The binary-effect and null-covariate edges are not in it. With the prior on all terms we expect:
  - a small gain for the continuous effect;
  - binary-effect edges outside the prior penalised as before, those at prior pairs favoured;
  - null-covariate false positives concentrated at prior-supported pairs.

**Files.**
- `slurm/settings_ablation.csv` (400 rows: tasks 1–100 n = 20k, 101–200 10k, 201–300 5k, 301–400 1k);
- `slurm/run_ablation.sbatch` (noisy and perfect priors, 16 cores, 1 h);
- `analyse_ablation.R` → `results_ablation/summary/`:
  - `ablation_summary.csv`;
  - `ablation_paired.csv` (all terms minus baseline only, paired Wilcoxon);
  - `Fig_ablation` (AND) and `Fig_ablation_OR`.

**Cost.** Local: about 6 s per node at n = 1k (2 priors: 2.5 min on 8 cores). On the cluster: about 1–7 min per task, < 100 core-hours in total.

**Submission** (from Simulation_analysis_v2):
```
sbatch --array=1-400%50 slurm/run_ablation.sbatch
# when finished
Rscript --vanilla analyse_ablation.R results_ablation results
```

### Ablation: results (full run, 100 replications per n, analysed locally 2026-11-25)
`analyse_ablation.R results_ablation results`: all 400 tasks complete and paired with the main-study fits on the same datasets. All numbers below are AND means; the differences are all terms minus baseline only, with paired Wilcoxon p < 1e-9 unless stated.

| n | Arm | F1 baseline | F1 binary | F1 continuous | FPR×1000 binary | FPR×1000 continuous | FPR×1000 null |
|---|---|---|---|---|---|---|---|
| 1k | No prior | 0.821 | 0.616 | 0.084 | 1.76 | 1.82 | 1.84 |
| 1k | Noisy, baseline only | 0.905 | 0.627 | 0.071 | 1.35 | 1.60 | 1.56 |
| 1k | Noisy, all terms | 0.888 | **0.204** | 0.068 | 0.35 | 0.18 | 0.16 |
| 1k | Perfect, baseline only | 0.958 | 0.654 | 0.080 | 1.75 | 2.16 | 2.16 |
| 1k | Perfect, all terms | 0.928 | **0.067** | 0.016 | 0.17 | 0.10 | 0.08 |
| 5k | Noisy, baseline only / all terms | 0.979 / 0.983 | 0.970 / **0.824** | 0.502 / 0.576 | 0.47 / **2.48** | 0.48 / 0.39 | 0.39 / 0.11 |
| 5k | Perfect, baseline only / all terms | 0.990 / 0.976 | 0.967 / **0.608** | 0.557 / **0.321** | 0.63 / **4.30** | 0.77 / 0.63 | 0.68 / 0.03 |
| 20k | Noisy, baseline only / all terms | 0.987 / 0.996 | 0.991 / **0.866** | 0.969 / 0.934 | 0.18 / **3.12** | 0.10 / 0.74 | 0.10 / 0.06 |
| 20k | Perfect, baseline only / all terms | 0.997 / 0.997 | 0.991 / **0.633** | 0.974 / **0.724** | 0.18 / **11.55** | 0.20 / **5.22** | 0.21 / 0.00 |

(The 10k results lie between 5k and 20k; all cells are in `results_ablation/summary/ablation_summary.csv`. With OR the pattern is the same: e.g. perfect prior at 20k, binary F1 0.917 → 0.432 and continuous F1 0.911 → 0.559.)

**Findings**
1. **Covariate-specific networks are damaged by applying the prior to all terms, at every n.**
   - Binary-effect F1 falls by 0.12–0.42 with the noisy prior and by 0.33–0.59 with the perfect prior, and does not recover with n: about 0.63 for the perfect prior at 20k, against 0.99 with the selective prior.
   - Continuous-effect F1 falls by 0.24–0.25 with the perfect prior at n ≥ 5k. This is despite its true edges being in the prior (T₀ = Δ₀ + mean(x1)·B1). With the noisy prior it is mixed: +0.07 at 5k, −0.03 at 20k.
2. **The failure has two regimes.**
   - Small n (1k): the covariate terms are over-shrunk. Binary TPR 0.52 → 0.12 (noisy) and 0.58 → 0.04 (perfect), and FPRs go down for every covariate.
   - Large n (≥ 5k): **the baseline structure leaks into the covariate networks**. Binary FPR×1000 0.18 → 3.1 (noisy) and → 11.6 (perfect) at 20k, and continuous FPR×1000 0.20 → 5.2 (perfect). The perfect-prior binary network has 105 estimated edges against 50 true.
   - The leak grows with n and is worst with the most informative (perfect) prior. This is the signature of the prior favouring the selection of covariate terms at baseline-prior pairs where the covariate has no effect.
   - Not verified edge by edge: `evaluation.rds` does not store edge lists. The interpretation follows from the FPR rising with n under the perfect prior while TPR is ≈ 1.
   - The small-n over-shrinkage is consistent with the pf normalisation (sparsegl rescales the penalty factors to mean 1, so terms outside the prior get pf > 1).
3. **The baseline network is barely affected.** |ΔF1| ≤ 0.03, largest at n = 1k (perfect −0.030). With the noisy prior, all terms gives a slightly higher baseline F1 at n ≥ 5k (+0.004 to +0.009), with lower FPR. The prior's benefit for the baseline is therefore the same whichever way it is applied. The difference is entirely in the covariate-specific networks.
4. **The null covariate is not a discriminating test.** All terms gives fewer null false positives (2.16 → 0.08 ×1000 at 1k, perfect), but this is part of the general over-shrinkage of covariate terms and comes with the loss of the true covariate effects. Our expectation of null false positives at prior pairs was wrong for this design.
5. **The selective prior (prCGGM) is never worse than no prior** on the covariate networks (binary F1 0.627/0.654 vs 0.616 at 1k; ≥ 0.967 vs 0.961 at 5k and above). So it keeps the baseline gain at no cost to the covariate effects.

**Suggested wording (main text, Results, with Fig_ablation_main).** "Restricting the prior to the baseline is essential. When the same prior weights were applied to the covariate-specific terms, the baseline network was essentially unchanged (|ΔF1| ≤ 0.03), but the covariate-specific networks deteriorated at every sample size. At n = 1,000 the binary-effect F1 dropped from 0.65 to 0.07 (perfect prior) because true effects were over-shrunk. At n = 20,000 its false-positive rate increased 60-fold (0.18 to 11.6 per 1,000 pairs) as baseline structure leaked into the covariate networks, and the F1 plateaued at 0.63 instead of 0.99. A more accurate prior made the leakage worse."

**Figures.**
- `results_ablation/summary/Fig_ablation_main.{png,pdf}` (`fig_ablation_main.R`; 4 panels, AND: baseline F1, binary F1, continuous F1, binary FPR) for main Fig 2;
- `Fig_ablation` (7 panels, AND) and `Fig_ablation_OR` for the supplementary.
