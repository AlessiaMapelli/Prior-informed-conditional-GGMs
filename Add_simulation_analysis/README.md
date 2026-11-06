# Additional simulation analysis (exploratory)

Why does the cross-validated lasso produce very dense networks when n is large, and which selection rule should be used instead? This folder is self-contained. It uses copies of the original code in `code/` and never modifies `Simulation_analysis/` or `Computational_templates/`.

## Design

- **Data**: `simulate_setting()` follows the design of Zhang & Li (2023, Section 6). It uses p = 100 and q = 3 covariates: x1 and x2 are Uniform[0, 1] and x3 is binary. Two of the three, chosen at random, carry an effect; the third (the null covariate) has none. The prior is noisy (sensitivity 0.8, specificity 0.9).
  - **Population network**: power-law degree distribution with exponent 2.5 (`igraph::sample_fitness_pl`), p edges, weights ±U[0.35, 0.5].
  - **Covariate networks**: Erdős–Rényi with edge probability 0.01 (as in Zhang & Li), weights ±U[0.175, 0.25], zero diagonal.
  - **Diagonal**: a single constant d for all nodes and all covariate values. It is the smallest value that keeps Ω(u) = dI + B₀ + Σ uₕBₕ positive definite (minimum eigenvalue 0.1) over the whole covariate range. Ω(u) is linear in u, so its minimum eigenvalue is lowest at a corner of the covariate box, and checking the 2^q corners is exact.
  - **Model assumptions**: the conditional variances 1/d do not depend on the covariates, so the node-wise coefficients are linear in u, as the model assumes. In a typical replicate, d ≈ 1.5, the baseline |partial correlation| has median 0.28, and the covariate effects at u = 1 have median 0.14.
- **First run (`results/`, kept for reference)**: this used the paper's generator (`generate_simulation_dataset()`), which has two problems.
  - `igraph::sample_pa(power = 2.5)` is super-linear preferential attachment and produces a star: one hub linked to 88–99 of the 100 nodes. Zhang & Li's "power law with parameter 2.5" refers to the degree distribution instead.
  - Δ₀ + Σ uₕΔₕ is not positive definite for any subject. A subject-specific constant (+1.6 to +3.8) is therefore added to the diagonal, and the covariate matrices also carry a diagonal. Both shrink the partial correlations severalfold and make the conditional variances depend on u.
  The current run writes to `results_v2/`.
- **Sample sizes**: n ∈ {1k, 5k, 10k, 20k}, with 5 replications each. Each node regression has (q+1)(p−1) = 396 predictors, so n/predictors ranges from 2.5 to 50. In the UKB application this ratio is about 16 (Sex and Diab as covariates).
- **Model fit**: the same design matrix, penalty factors, hyperparameter grid (α ∈ {0.5, 0.75, 0.9, 0.99}, ϑ ∈ {0.8, 1, 1.1, 1.3}) and 50 % random search as the paper's simulation. Five-fold CV uses one shared fold assignment.
- **Selection rules**: each node and hyperparameter combination is fitted once, and every rule is computed from that same path.

| Rule | λ | Hyperparameters | Cost |
|---|---|---|---|
| `min` | CV, `lambda.min` (current default) | BIC across combinations | path + CV |
| `1se` | CV, `lambda.1se` | BIC across combinations | path + CV |
| `bic` | BIC along the path | joint minimum | path only |
| `ebic` | extended BIC (γ = 0.5) along the path | joint minimum | path only |
| `stab` | stability selection (Meinshausen & Bühlmann 2010): 50 half-subsamples, π = 0.75, per-node PFER ≤ 1 | those of `min` | path + CV + 50 subsample paths |
| `min+backbone` | `min`, then the disparity-filter backbone (α = 0.05, BH) as in the T2D report | | negligible |
| `min+pcor` | `min`, then \|partial correlation\| > 0.05 | | negligible |

- **Coefficients**: every rule is evaluated with the penalised coefficients (`pen`) and after an OLS refit on the selected predictors (`refit`). Both symmetrisation rules, OR and AND, are evaluated.
- **Targets**: the baseline is scored against Δ₀, as in the paper, and against the centred target Δ₀ + Σ x̄ₕΔₕ. The covariates are centred before fitting, so the centred target is what the baseline actually estimates. This difference grows with n.
- **Metrics**:
  - density, TPR, FPR, F1;
  - Spearman over all pairs (the paper's MagCor);
  - Spearman and Pearson over true-positive edges only (`*_TP`), which measure the edge weights rather than the support;
  - Frobenius error;
  - computing time per node and rule.
- **CV diagnostics**: where `lambda.min` and `lambda.1se` fall on the path, how often they hit its end, the degrees of freedom at each, and the relative CV gain from `1se` to `min`. The full CV curves are stored for nodes 1–5.

## Hyperparameter tuning and selection rules (v3 run; reference for the paper)

This section describes the fitting procedure of the 100-replication run (`results_v3/`, §8) in the form to be used in the Methods and Supplementary material. The v2 description above differs only in the points listed in §8: covariate encoding, baseline target, and penalised coefficients.

### Node-wise model
Each node j is fitted separately. Node j is regressed on:
- the other p − 1 = 99 nodes (baseline block);
- their interactions with each covariate (one block per covariate: x1, x2, x3).

This gives P = (q + 1)(p − 1) = 396 predictors. Continuous covariates are standardised; the binary covariate enters as a 0/1 dummy. The model is fitted by the sparse-group lasso (`sparsegl`, no intercept).

### Hyperparameters

| | Role | Values |
|---|---|---|
| λ | overall penalty strength | path of about 100 values computed by `sparsegl` |
| α | balance between the lasso penalty (single edges) and the group penalty (a whole covariate block) | {0.5, 0.75, 0.9, 0.99} |
| ϑ | penalty factor of the covariate terms relative to the baseline terms | {0.8, 1, 1.1, 1.3} |

Fixed, not tuned:
- **Prior:** baseline predictor k of node j has penalty factor 1 − W_jk, so edges supported by the prior are penalised less. Covariate terms have penalty factor ϑ.
- **Group penalty factors:** 0 for the baseline block and 1 for each covariate block.
- **Random search:** each node evaluates a random 50% of the 16 (α, ϑ) pairs, i.e. 8 pairs, with a node-specific seed as in the template.
- **Cross-validation:** 5 folds, with one fold assignment shared by all nodes and pairs.

### Shared first step
For each of the 8 (α, ϑ) pairs, the full λ path is fitted on all n subjects:

```
                ┌─ (α,ϑ)₁ ── λ path ──┐
node j ──8 ×──  ├─ (α,ϑ)₂ ── λ path ──┤──► 8 × ~100 candidate models
                └─ ...                ┘
```

The selection rules differ only in how one model (or, for stability selection, one edge set) is chosen from these candidates.

### Selection rules

```
                    λ within each (α,ϑ)           choice of (α,ϑ)              additional fits
───────────────────────────────────────────────────────────────────────────────────────────────
CV lambda.min       5-fold CV → λ.min             BIC at λ.min                 5 CV fits per (α,ϑ)
  (= published method / current template)
CV lambda.1se       5-fold CV → λ.1se             BIC at λ.1se                 5 CV fits per (α,ϑ)
BIC                 ── jointly: minimum of BIC over all (α,ϑ,λ) ──             none
EBIC                ── jointly: minimum of EBIC over all (α,ϑ,λ) ──            none
Stability (CV)      no single λ (see below)       (α,ϑ) of CV lambda.min       5 CV fits per (α,ϑ)
                                                                               + 50 subsample paths
Stability (BIC)     no single λ (see below)       (α,ϑ) of BIC                 50 subsample paths
lambda.min + pcor   CV lambda.min, then remove edges with |partial correlation| ≤ 0.05
lambda.min + bb     CV lambda.min, then disparity-filter backbone (α = 0.05, Benjamini–Hochberg)
```

**Criteria:**
- BIC = n·log(RSS/n) + log(n)·df, with df the number of non-zero coefficients.
- EBIC = BIC + 2γ·log C(P, df) with γ = 0.5 (Chen & Chen, 2008). It penalises larger models more when P is large.

**Stability selection** (Meinshausen & Bühlmann, 2010):
1. Fix (α, ϑ), from CV lambda.min or from BIC, and the λ path of the full-data fit.
2. Draw B = 50 random subsamples of size ⌊n/2⌋ and fit the whole λ path on each.
3. For every predictor and every λ, record the selection frequency over the subsamples.
4. Restrict to the high-λ region of the path, where the average number of selected predictors is at most q_Λ = √(PFER·(2π − 1)·P) ≈ 14. Here PFER = 1, an expected number of false selections of at most about one per node, and π = 0.75.
5. Select the predictors with selection frequency ≥ π = 0.75 at some λ in that region.

Stability selection does not choose a single λ, so it only needs (α, ϑ). The two versions differ only in where (α, ϑ) comes from: from the lambda.min fit, which needs the cross-validation, or from the BIC fit, which does not.

**Post-selection thresholds** (applied to the CV lambda.min network after the refit):
- *pcor:* remove edges with |ρ_jk| ≤ 0.05, where ρ_jk = −δ_jk·σ_j·σ_k.
- *Backbone:* disparity filter (Serrano et al., 2009). An edge is kept if it is significant for at least one endpoint at α = 0.05 after Benjamini–Hochberg correction; this is the procedure used in the T2D report.

### After selection (all rules)

```
selected predictors ─► OLS refit (no penalty) ─► β̂_jk, σ̂²_j ─► δ_jk = −β̂_jk / σ̂²_j
                                                                    │
            all p nodes ─► symmetrise each block (OR: kept if selected in either node's
                           regression; AND: kept only if selected in both)
                                                                    │
                 Baseline (network at the mean of x1, x2 and x3 = 0), x1, x2, x3 networks
```

- **OLS refit:** an ordinary least-squares fit on the selected predictors only, without penalty, for each node. It removes the lasso shrinkage from the edge weights without changing the selected edges. Residual variance: σ̂²_j = RSS / (n − |S_j|), where S_j is the set of predictors selected for node j.
- **Precision scale:** δ_jk = −β̂_jk / σ̂²_j, because β_jk = −ω_jk/ω_jj and σ²_j = 1/ω_jj.
- **Symmetrisation:** applied to each block (baseline, x1, x2, x3) with the template's `build_delta_dictionary()`. The OR/AND decision is taken on the regression coefficients.

### Computing cost per node
Measured on the cluster (16 cores per setting), n = 1k → 20k, mean over 100 replications:

| Rule | Fits | s/node |
|---|---|---|
| BIC / EBIC | 8 paths | 6.8 → 30.0 |
| CV lambda.min / 1se | 8 paths + 5-fold CV for each pair | 46.8 → 190.5 |
| Stability (CV) | 8 paths + CV + 50 subsample paths | 105.6 → 303.8 |
| Stability (BIC) | 8 paths + 50 subsample paths | 62.2 → 144.7 |
| BIC / EBIC relaxed | 8 paths + Gram matrix + relaxed RSS of each distinct support | pending (expected +1–2 s over BIC/EBIC) |
| Stability (BIC relaxed) | as BIC relaxed + 50 subsample paths | pending |

**Relaxed criteria (added 2026-11-06, third pass):** BIC relaxed = n·log(RSS_OLS/n) + log(n)·df, where RSS_OLS is the residual sum of squares of the OLS refit of the support at each λ (instead of the RSS of the penalised fit); EBIC relaxed adds 2γ·log C(P, df). The penalised RSS is inflated by the lasso shrinkage of the true edges by an amount that grows with n, so BIC on the penalised RSS keeps choosing smaller λ and admits noise edges as n grows (selected λ ∝ 1/n); the relaxed RSS removes this bias. RSS_OLS = y'y − b'X_s'y with b = (X_s'X_s)⁻¹X_s'y, from the Gram matrix X'X computed once per node (QR fallback if singular); identical to `lm.fit` to 1e-14. See `Simulation_analysis_v2/README.md` (selection test) for the discovery and the 5-replication comparison.
| OLS refit | 1 least-squares fit per rule | < 0.05 |

`cv.sparsegl` refits each pair's path on the 5 training folds. The 50 subsample paths are fitted only for the chosen pair.

### Relation to the published method
The published method (and the current template) corresponds to **CV lambda.min without refit**: λ is chosen by 5-fold CV within each (α, ϑ) pair, the pair is chosen by BIC, and the penalised coefficients are reported. The supplementary analysis compares it with the alternatives above. These differ only in how λ (or, for stability selection, the edge set) is chosen; the OLS refit is applied to all of them.

## Files

| File | Content |
|---|---|
| `code/` | Unmodified copies of `Sim_functions.R`, `ggReg_main_functions.R` and the original job scripts |
| `R/explore_functions.R` | Data generation wrapper, node exploration under all rules, assembly and evaluation |
| `run_setting.R` | Runs one setting (n, replication) with the nodes in parallel; saves each node as soon as it is done, then evaluates |
| `slurm/settings.csv` | The 20 settings (n × replication), largest n first |
| `slurm/test_timing.sbatch` | Timing test: one node per core for each sample size |
| `slurm/run_settings.sbatch` | Array job, one task per setting, 16 cores per task |
| `analyse_results.R` | Tables and figures in `results/summary/` (v1/v2 runs) |
| `slurm/settings_v3.csv`, `slurm/run_settings_v3.sbatch` | The 400 settings of the 100-replication run (v3) and its array job |
| `analyse_results_v3.R` | Tables (mean (SD) over replications) and boxplots of the v3 run in `results_v3/summary/` |
| `archive_v2_code/` | The code as used for `results_v2/` (before the v3 changes) |

## Running on the cluster

From this folder on the cluster:

```bash
mkdir -p logs
# 1. Timing test (8 cores, nodes 1-8, n = 1000, 5000, 20000; results_test/)
sbatch slurm/test_timing.sbatch
sacct -j <jobid> --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS
# 2. Full run (writes to results_v2/)
sbatch --array=1-20%10 slurm/run_settings.sbatch
# 3. Summaries
Rscript analyse_results.R results_v2
# v3: 100 replications (writes to results_v3/), then its summaries
sbatch --array=1-400%20 slurm/run_settings_v3.sbatch
# v3 second pass: stability selection at the BIC hyperparameters (no CV)
sbatch --array=1-400%20 slurm/run_settings_v3_stab_bic.sbatch
# v3 third pass: BIC / EBIC relaxed and stability at the relaxed-BIC hyperparameters
sbatch --array=1-400%20 slurm/run_settings_v3_relaxed.sbatch
Rscript analyse_results_v3.R results_v3
```

A task that runs out of time can be resubmitted with the same index (`sbatch --array=3 ...`), and the nodes already fitted are skipped. Packages needed: sparsegl, igraph, MASS, matrixcalc, reshape2, ggplot2, and dplyr/tidyr for the analysis.

Local reference (Apple M-series, one core, n = 1000): about 70 s per node with all rules, of which 3 s for the path, 21 s for CV over the 8 combinations, and 34 s for 50 stability subsamples.

## Findings log (for the manuscript revision)

### 1. Code versions used by the paper analyses (checked 2026-11-02)
- `Simulation_analysis/` and `T2D_application_code/` ran on template commit 4cfcecd. `T2D_application_prot/` (the new T2D analysis, AUC 0.78) ran on commit 330a36e. All three now contain the current `Computational_templates/ggReg_main_functions.R`.
- **T2D networks are unchanged by the code fixes**:
  - the old name-based predictor selection picks exactly the same predictors as the new positional one for all 366 proteins;
  - screening is off, both λ come from CV, and both covariates (Sex, Diab) are factors, so the screening and covariate-scaling fixes don't apply;
  - the downstream analyses only use `Dic_adj_matrics` and `Cov_effect`;
  - node 279 refitted with the current template matches the saved fit to about 1e-16 (same λ, α, ϑ, BIC and 766 non-zero coefficients).
- The current code names covariate networks `<variable>_<level>` (e.g. `Diab_Incident`), while the paper's run produced `DiabIncident`. `T2D_application_code/` (Rmd params, `enrichment_communities.py`, README) now uses `Diab_Incident`. `T2D_application_prot/` keeps `DiabIncident` because its saved results use that name.
- On 2026-10-28 a 10-variable test overwrote `T2D_application_prot/results/input_data_nodes.rda` and node files 1–10. The collected `PPI_T2D_app.RData` (2026-08-05) is intact.

### 2. Paper simulation regenerated with the δ fix (2026-11-03, job 45334610)
- `Simulation_analysis/Sim_regenerate_results.R` re-collects all 240 replications from the saved node files, with δ = −β/σ² and the OR/AND choice made on β. Nothing is refitted.
- TPR, FPR, F1 and Accuracy are identical to before. Seven of the 72 MagCor cells change in the third decimal (at most ±0.001). The Frobenius error changes slightly but is not reported in the paper.
- The regenerated values agree with an independent recompute to about 1e-16. The previous outputs are in `Simulation_results/archive_before_delta_fix/`.

### 3. Problems found in the simulation data generator (`Sim_functions.R`, also used by the paper)
- **The baseline network is a star.** `igraph::sample_pa(power = 2.5)` is super-linear preferential attachment. In the paper's simulation one hub is linked to 46 of 50 nodes at p = 50; in the first exploratory run, to 88–99 of 100 nodes. Zhang & Li draw the *degree distribution* from a power law with exponent 2.5.
- **No subject's precision matrix is positive definite.** Δ₀ + Σ uₕΔₕ has minimum eigenvalue −1.7 to −3.5, so `simulate_subject_data_robust()` adds a subject-specific constant to the diagonal. Over 100% of subjects this averaged +1.6 to +3.1 in the paper's configurations and +3.8 in the first exploratory run, against a baseline diagonal of 0.6. On top of that, the covariate matrices Δₕ have a diagonal of about 0.8–1.2.
- **Consequences**:
  - the true partial correlations shrink severalfold (baseline about 0.09, covariate effects 0.02–0.05), which explains why covariate networks were hard to detect in the paper (F1 about 0.1);
  - the conditional variances depend on u, violating the model's assumption (Zhang & Li, Section 2: Ω(u)_jj = σ_jj; their Supplement S1.1 tolerates σ_jj(u) ∈ [1, 2]), so the node-wise coefficients are not linear in u.
- **Diagonal schemes compared** on a power-law baseline (|ρ| medians, baseline / covariate):
  - (A) Zhang & Li's row-rescaling with σ_jj = 1: 0.07 / 0.07 at edge probability 0.05, and 0.17 / 0.19 at 0.01. Positive definiteness is not guaranteed: it failed on a star network, because averaging after rescaling can break diagonal dominance.
  - (B) a common diagonal, the smallest that is PD at the covariate corners: 0.20 / 0.10, and 0.24 / 0.12 at 0.01.
  - (C) node-specific row dominance: 0.11 / 0.06.
  - **Chosen: (B)**, with a power-law baseline (exponent 2.5) and covariate edge probability 0.01 (as Zhang & Li).
- The package's `simulate_prcggm_data()` has the same two issues: attachment probability ∝ degree^2.5, and a per-subject correction for 46% of subjects with the defaults. **To fix after the analyses are settled.**
- The covariates are centred before fitting, so the estimated baseline is Δ₀ + Σ x̄ₕΔₕ, not Δ₀. In the paper's simulation this explains only part of the baseline false positives (covariate edges are selected in the baseline about 6% of the time vs 3% for true zeros), but it grows with n. The exploratory study therefore scores the baseline against both targets.

### 4. Density diagnosis in the T2D application (`T2D_application_prot`)
- Median selected λ_prec ≈ 1.6e-6 (range 4.6e-8 to 1.4e-4).
- Network density: baseline 0.97, Sex 0.92, Diab 0.90.
- ϑ = 1.5, the largest grid value, is chosen for 263 of 366 nodes.

### 5. First exploratory run (`results/`, star generator; superseded, kept for reference)
- `lambda.min` never reaches the end of the λ path, but moves deeper into it as n grows (median position 0.13 → 0.40). Its baseline F1 (OR) falls from 0.60 to 0.11 between n = 1k and 20k, and density rises from 0.05 to 0.35 (truth 0.02). The CV gain over `lambda.1se` is under 1%. It also puts 19% density on the null covariate at n = 20k.
- `1se`, BIC, EBIC and stability selection never select edges on the null covariate. Their baseline F1 (OR) degrades more slowly with n: `1se` 0.86 → 0.63, EBIC 0.83 → 0.62, BIC and stability selection about 0.79 → 0.51.
- AND strongly helps BIC and EBIC (EBIC F1 0.94–0.98 across n), but hurts `1se` and stability selection.
- The partial-correlation threshold and the backbone after `lambda.min` reached F1 ≈ 0.99, which is likely an artefact of the star network and the uniformly small partial correlations.
- Computing time per node (n = 1k → 20k): about 100 → 290 s for the CV rules, 14 → 44 s for BIC/EBIC, and 240 → 480 s for stability selection. The refit is negligible.
- Timing test on the cluster (cnode01, 8 cores, job 45315243): about 215 s per node per core at n = 1k–5k and 425 s at n = 20k, with about 2 GB per worker at n = 20k.

### 6. Exploratory run with the corrected generator (`results_v2/`, 2026-11, 20 settings × 100 nodes, all complete)
Setting: p = 100, q = 3 (two covariates with an effect, one null), noisy prior, 5 replications per n. Results are means over replications. F1 for the baseline is computed against the **centred target** Δ₀ + Σ x̄ₕΔₕ, which is what the baseline estimates because the covariates are centred. Against Δ₀ alone, every rule loses F1 as n grows because it correctly detects the shifted edges.

**CV diagnostics.** λ never reaches the end of the path. `lambda.min` moves deeper into the path as n grows (median position 0.13 → 0.35, median non-zero coefficients per node 6 → 38), while the CV error gained over `lambda.1se` stays at 0.9–1.8%.

**F1, OR symmetrisation, penalised coefficients** (n = 1k / 5k / 10k / 20k):

| Rule | Baseline (centred) | Covariate with effect | Density on the null covariate |
|---|---|---|---|
| min | 0.60 / 0.56 / 0.47 / 0.36 | 0.16 / 0.34 / 0.28 / 0.18 | 0.01 / 0.03 / 0.04 / 0.08 |
| 1se | 0.56 / 0.87 / 0.93 / 0.95 | 0.00 / 0.35 / 0.37 / 0.65 | 0 |
| bic | 0.69 / 0.86 / 0.83 / 0.82 | 0.05 / 0.70 / 0.79 / 0.78 | ≤ 0.002 |
| ebic | 0.69 / 0.90 / 0.87 / 0.85 | 0.02 / 0.62 / 0.77 / 0.86 | ≤ 0.001 |
| stab | 0.71 / 0.94 / 0.95 / 0.95 | 0.06 / 0.71 / 0.85 / 0.97 | 0 |
| min + pcor > 0.05 | 0.70 / 0.83 / 0.87 / 0.92 | 0.03 / 0.23 / 0.17 / 0.13 | 0 |
| min + backbone | 0.08 / 0.38 / 0.58 / 0.71 | 0.00 / 0.23 / 0.58 / 0.86 | 0 |

**With AND**:
- BIC, EBIC and stability selection reach baseline F1 of 0.94–0.99 for n ≥ 5k. BIC with AND gives covariate F1 of 0.55 / 0.69 / 0.95 at n = 5k / 10k / 20k.
- AND lowers the power for covariate effects for most rules, and `1se` with AND is too conservative.

**Refit (OLS on the selected predictors).**
- It doesn't change the support, so it doesn't change F1.
- It improves the edge weights on true positives. Pearson with the truth rises from 0.95–0.97 to 0.96–0.98 on the baseline, and from 0.77–0.96 to 0.96–0.99 on the covariate effects. The largest gains are for `1se` and EBIC, which shrink the most.
- Its cost is under 0.03 s per node.
- The paper's MagCor (Spearman over all pairs, including zeros) mostly reflects the support and doesn't change with the refit.

**Computing time per node** (n = 1k / 5k / 10k / 20k, seconds on one cluster core):
- CV rules (`min`, `1se`): 34 / 64 / 103 / 174.
- BIC and EBIC (path only): 5 / 10 / 16 / 27, about **6–7× cheaper** than CV.
- Stability selection (50 subsamples, on top of CV): 79 / 106 / 167 / 273.

**n = 1k** (n / predictors ≈ 2.5): the baseline is recovered (F1 about 0.7), but no rule detects the covariate effects reliably.

### 7. Draft guidelines (to be confirmed by the paper simulations)
1. **Do not use `lambda.min`.** It over-selects more as n grows (here baseline density 0.06 → 0.18 against a true 0.04), puts spurious edges on covariates without any effect, and improves the CV error by under 2%. This is the mechanism behind the dense T2D networks.
2. **Default: λ by BIC along the path, chosen jointly with α and ϑ.** It needs no cross-validation (6–7× faster), and it is consistent with the BIC already used to tune α and ϑ. It gives near-zero false positives on the null covariate and good power for covariate effects. EBIC is the more conservative alternative when p is large relative to n.
3. **Stability selection** is the most accurate overall, especially for covariate effects at large n. It costs about 1.6× CV and 10× BIC, so offer it as an option for final analyses when compute allows.
4. **`lambda.1se`** is accurate for the baseline when n is large but has low power for covariate effects. Not recommended as the default.
5. **Symmetrisation:** OR by default with BIC/EBIC/stability selection. AND increases baseline precision but reduces power for covariate effects.
6. **Refit:** always refit the selected edges by OLS to report edge weights. It is essentially free and removes the lasso shrinkage.
7. **Post-hoc thresholds are not useful.** The partial-correlation threshold removes real covariate effects, and the disparity backbone is unstable.
8. **Interpretation:** with centred covariates, the baseline network is the network at the average covariate profile, not at u = 0. The paper should state this explicitly, and simulations should score against that target.

Caveats: one setting (p = 100, q = 3, noisy prior, ER edge probability 0.01). The paper simulations (step 3) should vary the prior type, p and the effect sizes.

### 8. Supplementary sensitivity analysis: 100 replications (`results_v3/`)
Purpose: supplementary material justifying the choices made for the paper simulations (selection rule, symmetrisation, refit, no post-hoc thresholds). Same setting as v2 (p = 100, q = 3, noisy prior, n ∈ {1k, 5k, 10k, 20k}, generator of §3 with scheme B), 100 replications per n (seed 1000·rep + n). Changes with respect to v2 (code of v2 in `archive_v2_code/`):
- **Covariate encoding as in the template and the package:** continuous covariates (x1, x2) are standardised, the binary one (x3) enters as an unscaled 0/1 dummy (reference level 0). In v2 x3 was centred as well.
- **Baseline target and prior:** the estimated baseline is the network at the mean of the continuous covariates and at the reference level of the binary ones, T₀ = Δ₀ + x̄₁Δ₁ + x̄₂Δ₂. The baseline is scored against T₀ only, and the prior (sensitivity 0.8, specificity 0.9) is generated from the support of T₀, so that it describes the network the baseline estimates. The generator itself is unchanged.
- **Truth on the scale of the fit:** the proteins are standardised (precision S Ω S, S = marginal SDs) and the continuous covariates too (component × sd(x_h)); the true components are expressed on this scale, so edge weights and Frobenius errors are directly comparable. Checked: slope of the refitted estimates on the truth over the true positives 1.00 for the baseline.
- **OLS refit only:** all rules (CV lambda.min, lambda.1se, BIC, EBIC, stability selection) use the OLS refit on the selected predictors; lambda.min + partial-correlation threshold (|ρ| > 0.05) and lambda.min + disparity backbone are applied to the refitted lambda.min network. OR and AND.
- **Metrics:** density, TPR, FPR, precision, F1, Spearman and Pearson correlation with the truth on the true-positive edges, Frobenius error (off-diagonal). MagCor and the CV diagnostics are dropped; computing times are kept.
- **Storage:** node files and simulated data are deleted after each evaluation; each `rep*/` keeps `evaluation.rds`, `times.rds` and `truth.rds` (true components, prior, configuration).

**Run:** job 45365482 (array 1–400), all 400 settings evaluated, no failed nodes. True baseline T₀ has on average 162–166 edges (100 from Δ₀ plus the edges of the effective continuous covariates); every replication has two effective covariates (2 effect components and 1 null component). Outputs in `results_v3/summary/` (tables `table_<metric>_<OR|AND>.csv`, `covariate_effects_by_type.csv`, `computing_time.csv`, figures `Fig_*.pdf/png`).

**F1, mean (SD) over 100 replications, n = 1k / 5k / 10k / 20k**

| Method | Baseline OR | Baseline AND | Covariate effect OR | Covariate effect AND |
|---|---|---|---|---|
| CV lambda.min | 0.70 / 0.61 / 0.50 / 0.38 | 0.87 / 0.89 / 0.82 / 0.71 | 0.22 / 0.39 / 0.30 / 0.20 | 0.14 / 0.62 / 0.70 / 0.58 |
| CV lambda.1se | 0.77 / 0.90 / 0.94 / 0.94 | 0.42 / 0.71 / 0.82 / 0.90 | 0.02 / 0.27 / 0.36 / 0.58 | 0.00 / 0.13 / 0.24 / 0.38 |
| BIC | 0.87 / 0.90 / 0.88 / 0.83 | 0.90 / 0.97 / 0.97 / 0.95 | 0.17 / 0.56 / 0.77 / 0.80 | 0.07 / 0.43 / 0.66 / 0.92 |
| EBIC | 0.90 / 0.93 / 0.91 / 0.86 | 0.88 / 0.97 / 0.98 / 0.96 | 0.10 / 0.49 / 0.72 / 0.84 | 0.03 / 0.36 / 0.54 / 0.88 |
| Stability selection | 0.90 / 0.93 / 0.93 / 0.92 | 0.90 / 0.97 / 0.98 / 0.98 | 0.16 / 0.60 / 0.84 / 0.94 | 0.07 / 0.42 / 0.68 / 0.92 |
| lambda.min + pcor | 0.73 / 0.97 / 0.99 / 0.99 | 0.85 / 0.96 / 0.98 / 0.99 | 0.22 / 0.48 / 0.45 / 0.42 | 0.14 / 0.47 / 0.45 / 0.43 |
| lambda.min + backbone | 0.00 / 0.01 / 0.23 / 0.57 | 0.00 / 0.00 / 0.00 / 0.07 | 0.00 / 0.00 / 0.00 / 0.08 | 0.00 / 0.00 / 0.00 / 0.00 |

SDs: baseline ≤ 0.05 for all rules except the backbone; covariate effects 0.02–0.25, largest at n = 5k–10k (see below).

**Covariate effects, binary vs continuous (F1, OR / AND).** The large SDs of the covariate effects come from the covariate type: the binary covariate (SD 0.5) carries more signal per unit effect than the U[0,1] ones (SD 0.29).
- Binary: BIC 0.49 / 0.87 / 0.83 / 0.74 (OR) and 0.22 / 0.94 / 0.97 / 0.96 (AND); stability selection 0.48 / 0.94 / 0.94 / 0.91 (OR) and 0.21 / 0.95 / 0.99 / 0.98 (AND); lambda.min 0.47 / 0.44 / 0.32 / 0.21 (OR).
- Continuous: BIC 0.02 / 0.41 / 0.74 / 0.83 (OR) and 0.00 / 0.16 / 0.51 / 0.90 (AND); stability selection 0.02 / 0.42 / 0.79 / 0.96 (OR) and 0.00 / 0.15 / 0.53 / 0.89 (AND); lambda.min 0.10 / 0.36 / 0.30 / 0.20 (OR) and 0.02 / 0.50 / 0.67 / 0.57 (AND).

**Other metrics**
- *False positives:* lambda.min (OR) baseline FPR 0.03 → 0.11 and density 0.06 → 0.14 (true 0.033); null-covariate density 0.01 → 0.07. All other rules: null-covariate density 0.00 at every n; baseline FPR ≤ 0.01.
- *Precision vs recall:* with OR, BIC/EBIC lose some precision at large n (baseline precision 0.84 → 0.71 for BIC, recall 1.00), which is why their baseline F1 decreases slightly at 20k; AND restores precision (0.91–0.97). Stability selection keeps precision 0.85–0.88 (OR) at all n.
- *Edge weights after the OLS refit:* on the true positives, Pearson 0.99–1.00 for the baseline and 0.94–0.99 for covariate effects with BIC, EBIC and stability selection (Spearman 0.98–1.00 baseline, 0.81–0.94 covariate effects, increasing with n). Once an edge is selected, the refitted weight is accurate for every rule; the rules differ in selection only.
- *pcor threshold* (|ρ| > 0.05 on lambda.min): excellent for the baseline at n ≥ 5k (F1 0.96–0.99), but caps covariate-effect recall at about 0.40–0.47: the covariate effects are smaller than the threshold on the partial-correlation scale, so a fixed threshold cannot be calibrated for both.
- *Disparity backbone:* fails on these sparse networks (baseline recall 0.00–0.40). With 3–5 edges per node, the null model of the disparity filter rarely rejects; it is not suitable for sparse estimated networks.
- *n = 1000* (n / predictors per node ≈ 2.5): baseline F1 0.87–0.90 with BIC/EBIC/stability; covariate effects essentially undetectable for continuous covariates (F1 ≤ 0.10) and partially for the binary one (≤ 0.49).

**Computing time** (s per node, mean (SD); setting wall time on 16 cores)

| n | BIC/EBIC | CV rules (path + CV) | Stability selection (path + CV + 50 subsamples) | OLS refit | Setting wall time (min) |
|---|---|---|---|---|---|
| 1,000 | 6.8 (0.7) | 46.8 (4.8) | 105.6 (11.1) | < 0.05 | 12.0 (1.2) |
| 5,000 | 13.0 (1.4) | 87.8 (9.4) | 148.7 (16.1) | < 0.05 | 16.7 (1.8) |
| 10,000 | 17.9 (2.1) | 116.4 (13.7) | 195.8 (23.8) | < 0.05 | 22.0 (2.6) |
| 20,000 | 30.0 (3.3) | 190.5 (20.8) | 303.8 (33.8) | < 0.05 | 34.2 (3.8) |

(Stability selection uses the hyperparameters chosen by CV lambda.min here, hence the CV cost; running it at the BIC hyperparameters would cost about path + subsamples.)

**Second pass: stability selection at the BIC hyperparameters** (`--mode stab_bic`, `slurm/run_settings_v3_stab_bic.sbatch`). Stability selection uses the whole λ path on each subsample, so it needs only α and ϑ; in the first pass these came from the lambda.min fit, which requires the cross-validation. The second pass takes them from the BIC fit instead (no cross-validation), on the same 400 datasets regenerated from the seed (checked against `truth.rds`; the BIC/EBIC results of the second pass must equal the first pass, which `analyse_results_v3.R` checks). Outputs `evaluation_stab_bic.rds` and `times_stab_bic.rds` per replication; the summary reports both versions ("Stability selection (CV hyperparameters)" and "(BIC hyperparameters)"). Completed (job 45378868, all 400; BIC rows identical to the first pass, 3,200 rows). F1, OR, n = 1k / 5k / 10k / 20k: baseline 0.90 / 0.94 / 0.94 / 0.93 (CV hyperparameters 0.90 / 0.93 / 0.93 / 0.92); covariate effects 0.16 / 0.54 / 0.80 / 0.94 (CV 0.16 / 0.60 / 0.84 / 0.94); AND baseline identical to the CV version. Time per node 62 → 145 s against 106 → 304 s.

**Third pass: relaxed criteria** (`--mode relaxed`, `slurm/run_settings_v3_relaxed.sbatch`, 2026-11-06). Same 400 datasets (checked against `truth.rds`), no cross-validation: BIC/EBIC on the penalised RSS (must equal the first pass), BIC relaxed, EBIC relaxed, and stability selection at the hyperparameters of BIC relaxed. Outputs `evaluation_relaxed.rds`, `times_relaxed.rds`. `analyse_results_v3.R` now also: splits the covariate effects into binary and continuous in all tables and figures; writes the overall metrics (`overall_{F1,FPR,Pearson}_{OR,AND}.csv`: all components of a fit pooled, F1 and FPR from the summed TP / FP / FN / TN, Pearson as the TP-weighted mean of the per-component correlations); adds the relaxed rules to the cost table. Checked locally (n = 1k, rep 1): data identical to `truth.rds`, BIC/EBIC identical to the first pass (16 rows). Results: pending. This is Part 1 of the supplementary sensitivity analysis (paper structure in `Simulation_analysis_v2/README.md`).

**Conclusions for the paper (update of §7)**
1. lambda.min should not be used: it is the only rule whose false positives grow with n, on the baseline and on null covariates, and its baseline F1 falls from 0.70 to 0.38.
2. BIC along the path is the recommended default: no cross-validation (6–7× faster than CV), near-zero false positives on null covariates, baseline F1 0.87–0.97 and covariate-effect F1 up to 0.80 (OR) / 0.92 (AND) at 20k.
3. Stability selection is the most accurate overall (covariate effects up to 0.94 with OR) and the most stable across n, at about 10× the cost of BIC: recommended as an option for final analyses.
4. lambda.1se has low power for covariate effects (F1 ≤ 0.58) and is not recommended.
5. Symmetrisation: AND gives the best baseline (BIC/EBIC/stability 0.95–0.98 at n ≥ 5k) and the best covariate effects at 20k, but lower power for continuous effects at 5k–10k; OR gives higher recall for covariate effects at moderate n. The choice depends on whether precision or recall is preferred; in this run AND is better for n ≥ 20k and for the baseline at every n ≥ 5k.
6. OLS refit: always, since with it the edge weights of the selected edges are very accurate (Pearson ≥ 0.94).
7. Post-hoc thresholds and backbone: not recommended (the threshold removes covariate effects; the backbone fails on sparse networks).
8. Design point: covariate effects need roughly n ≥ 5k for a binary covariate and n ≥ 10k for a continuous one (at these effect sizes, p = 100).

