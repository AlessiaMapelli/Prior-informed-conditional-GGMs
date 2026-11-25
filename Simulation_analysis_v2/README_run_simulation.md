# Reproducing the simulation results of the paper

This file explains how to regenerate every simulation figure and table of the paper (main text and Supplementary) from scratch. The design, the decisions and the full results are in `README.md` (this folder). The exploratory analyses that are not reported in the paper are in `../Add_simulation_analysis/`.

## Software

- **Estimator:** the single template `../Computational_templates/ggReg_main_functions.R`, sourced by every script.
  - It implements the same estimator as the R package prCGGM ≥ 0.2.0. `../Computational_templates/check_template_vs_package.R` checks this.
  - It also has two template-only options, used only in the Supplementary: the post-selection filters `post_filter = "pcor"` / `"backbone"`, and `prior_on_covariates = TRUE`. See `../Computational_templates/README.md`, "Template and package".
- **R ≥ 4.3** with:
  - estimation: `sparsegl`, `glmnet`, `parallel`, `doParallel`, `foreach`;
  - data generation and evaluation: `igraph`, `MASS`, `matrixcalc`;
  - figures: `ggplot2`, `reshape2`, `dplyr`, `tidyr`.
- **Selection-rules comparison:** also uses the generator helpers of the preprint simulation (`../Preprint_simulation_analysis/Sim_functions.R`: `generate_prior_knowledge()`, `calculate_binary_metrics()`).
- **Cluster:** the scripts were run on a SLURM cluster with R in a Singularity container (`source /center/healthds/singularity_functions` in the `.sbatch` files; replace that line with your cluster's R setup). Each array task uses 16 cores and runs the node regressions in parallel with forked workers (`parallel::mclapply`).
- **Folders:** every command below is run from `Simulation_analysis_v2/` unless stated otherwise. Results go to `results*/` folders, which git ignores.

## Seeds

All datasets are regenerated from seeds, so no data files are needed.

| Study | Seed |
|---|---|
| Main study, EBIC / stability reruns, selection test Part 2, ablation | `dataset_seed(n, p, q, rep)` in `R/sim_design.R` (the same dataset for every rule and prior) |
| Reproducibility | graph `5e7 + 10p + q` (fixed network); data `6e7 + 1e6·index(n) + 1000p + 100q + rep`; prior: graph seed + 7 |
| Selection rules (Part 1) | `1000·rep + n` (`selection_rules/run_setting.R`) |
| Within a fit | random hyperparameter search `234 + 10·node`; stability subsamples `5000 + node` |

## Paper item → how to reproduce

The figure and table files are written by `paper_figures.R` (→ `paper_figures/`) unless stated otherwise.

| Paper item | Content | Fits (cluster) | Summaries | Figure / table |
|---|---|---|---|---|
| **Main Fig. 2** (prior) | F1 / TPR / FPR per component, q = 3, AND, three priors | **A** | `analyse_results.R` | `Fig2_F1_q3`, `Fig3_TPR_FPR_q3` |
| **Main Fig. 2**, ablation panel | Prior on the baseline only vs on all terms | **A** + **F** | `analyse_ablation.R results_ablation results` | `fig_ablation_main.R` → `results_ablation/summary/Fig_ablation_main` |
| **Main Fig. 3** (selection rules) | Overall F1 / FPR of the rules (Part 1) + BIC vs EBIC relaxed (and stability) at n = 1000 | **D**, **B**, **C** | `selection_rules/analyse_results_v3.R`; `analyse_results.R` on `results_ebic_relaxed`, `results_stability` | `FigS_rules_overall`, `FigS_n1000` |
| **Main Table 1** (guideline) | Which rule when | from the summaries of A–E | — | written in the text |
| Suppl.: computing time | Time per node by n, p, q | **A** | `analyse_results.R` | `FigS_time` |
| Suppl.: full main study | q = 1; other q = 3 metrics; AND vs OR; Pearson on TPs; all metrics | **A** | `analyse_results.R` | `FigS_q1`, `FigS_q3_other`, `FigS_AND_OR`, `FigS_Pearson`, `TableS_main_full.csv`, `TableS_main_compact.csv` |
| Suppl.: selection rules, per component (Part 1) | All rules incl. CV, stability, pcor threshold, backbone | **D** | `selection_rules/analyse_results_v3.R` | `FigS_rules_components` |
| Suppl.: robustness across priors and grid (Part 2) | BIC / EBIC relaxed under three priors, current vs extended grid | **E** | `selection_test/analyse_selection_test.R`, `selection_test/overall_metrics.R` | `TableS_part2.csv` |
| Suppl.: operating range at n = 1000 | BIC relaxed vs EBIC relaxed vs stability, effect of the prior | **A**, **B**, **C** | as above | `FigS_n1000`, `FigS_n1000_prior` |
| Suppl.: reproducibility | Jaccard, edge-selection frequency, density variation across independent datasets | **G** | `analyse_repro.R results_reproducibility` | `results_reproducibility/summary/Fig_repro_*` |
| Suppl.: ablation (full) | All metrics, AND and OR | **A** + **F** | `analyse_ablation.R` | `results_ablation/summary/Fig_ablation`, `Fig_ablation_OR` |

## Cluster runs

The core-hours are measured on 16-core tasks; the details are in `README.md`. A task that runs out of time can be resubmitted with the same index, and the parts already done are skipped. Logs go to `logs/`, so create the folder first (`mkdir -p logs selection_rules/logs`).

**A. Main study** (`run_task.R`, `slurm/settings.csv`; 3,200 tasks = n ∈ {1k, 5k, 10k, 20k} × p ∈ {50, 100, 250, 500} × q ∈ {1, 3} × 100 replications; each dataset fitted under the three priors, BIC relaxed, OR and AND evaluated; about 43,000 core-hours, 70% of it at p = 500):
```
sbatch --array=1-800%50 --export=ALL,OFFSET=2400 --time=00:30:00 --mem-per-cpu=2G slurm/run_tasks.sbatch   # p = 50
sbatch --array=1-800%50 --export=ALL,OFFSET=1600 --time=01:00:00 --mem-per-cpu=3G slurm/run_tasks.sbatch   # p = 100
sbatch --array=1-800%50 --export=ALL,OFFSET=800  --time=03:00:00 --mem-per-cpu=3G slurm/run_tasks.sbatch   # p = 250
sbatch --array=1-800%50 --export=ALL,OFFSET=0    --time=10:00:00 --mem-per-cpu=4G slurm/run_tasks.sbatch   # p = 500
sbatch slurm/analyse.sbatch                       # -> results/summary/
```
Optional timing test first: `sbatch --array=1-5 slurm/timing_test.sbatch`.

**B. EBIC relaxed at n = 1000** (all p and q; the n = 1000 tasks are indices 601–800 of each block):
```
sbatch --array=601-800%50 --export=ALL,OFFSET=2400,SELECTION=ebic_relaxed,OUT=results_ebic_relaxed --time=00:30:00 --mem-per-cpu=2G slurm/run_tasks.sbatch
sbatch --array=601-800%50 --export=ALL,OFFSET=1600,SELECTION=ebic_relaxed,OUT=results_ebic_relaxed --time=01:00:00 --mem-per-cpu=3G slurm/run_tasks.sbatch
sbatch --array=601-800%50 --export=ALL,OFFSET=800,SELECTION=ebic_relaxed,OUT=results_ebic_relaxed  --time=04:00:00 --mem-per-cpu=3G slurm/run_tasks.sbatch
sbatch --array=601-800%50 --export=ALL,OFFSET=0,SELECTION=ebic_relaxed,OUT=results_ebic_relaxed    --time=10:00:00 --mem-per-cpu=4G slurm/run_tasks.sbatch
sbatch --export=ALL,RES=results_ebic_relaxed slurm/analyse.sbatch
```

**C. Stability selection where P/n ≥ 0.5** (n = 1000: p = 250, q = 3; p = 500, q = 3; p = 500, q = 1; replications 1–10, paired with A and B):
```
sbatch --array=601-610%10 --export=ALL,OFFSET=800,SELECTION=stability_bic,OUT=results_stability --time=06:00:00 --mem-per-cpu=3G slurm/run_tasks.sbatch    # p = 250, q = 3 (pilot; command reconstructed: tasks 1401-1410)
sbatch --array=601-610%10 --export=ALL,OFFSET=0,SELECTION=stability_bic,OUT=results_stability --time=14:00:00 --mem-per-cpu=4G slurm/run_tasks.sbatch      # p = 500, q = 3
sbatch --array=701-710%10 --export=ALL,OFFSET=0,SELECTION=stability_bic,OUT=results_stability --time=20:00:00 --mem-per-cpu=4G slurm/run_tasks.sbatch      # p = 500, q = 1
sbatch --export=ALL,RES=results_stability slurm/analyse.sbatch
```

**D. Selection rules, Part 1** (`selection_rules/`; p = 100, q = 3, noisy prior, n ∈ {1k, 5k, 10k, 20k}, 100 replications, 400 tasks).
- Each node is fitted once per hyperparameter combination and every rule is read off the same paths (`selection_rules/R/explore_functions.R`).
- There are three passes over the same datasets:
  - CV lambda.min / 1se, BIC / EBIC, stability at the CV hyperparameters, and the pcor and backbone filters on the lambda.min networks;
  - stability at the BIC hyperparameters;
  - BIC / EBIC relaxed and stability at the relaxed-BIC hyperparameters.
```
cd selection_rules
sbatch --array=1-400%20 slurm/run_settings_v3.sbatch            # pass 1 (2 h per task)
sbatch --array=1-400%20 slurm/run_settings_v3_stab_bic.sbatch   # pass 2 (1 h)
sbatch --array=1-400%20 slurm/run_settings_v3_relaxed.sbatch    # pass 3 (1 h)
Rscript --vanilla analyse_results_v3.R results_v3               # -> results_v3/summary/
cd ..
```

**E. Selection test, Part 2** (`selection_test/`; the main-study datasets, p = 100, q = 3, three priors, 400 tasks; about 20 min per n = 20k task):
```
sbatch --array=1-400%20 selection_test/run_part2.sbatch
Rscript --vanilla selection_test/analyse_selection_test.R selection_test/results_100
Rscript --vanilla selection_test/overall_metrics.R selection_test/results_100
```

**F. Ablation of the selective prior** (`run_task.R --prior_on_covariates TRUE`; the A datasets at p = 100, q = 3, noisy and perfect prior, 400 tasks, < 100 core-hours):
```
sbatch --array=1-400%50 slurm/run_ablation.sbatch
Rscript --vanilla analyse_ablation.R results_ablation results
Rscript --vanilla fig_ablation_main.R results_ablation
```

**G. Reproducibility** (`run_repro_task.R`; a fixed true network, 50 independent datasets per n ∈ {1k, 5k, 10k}, BIC relaxed, EBIC relaxed and stability selection, noisy prior and no prior, 150 tasks, about 1,200 core-hours):
```
sbatch --array=1-150%30 slurm/run_repro.sbatch
Rscript --vanilla analyse_repro.R results_reproducibility
```

**Figures and tables** (after A–E):
```
Rscript --vanilla paper_figures.R paper_figures 0.5
```
The second argument greys out the settings with P/n ≥ 0.5, which are outside the operating range of BIC relaxed.

## Quick local check

A single task runs on a laptop, for example p = 50, q = 1, n = 1000, replication 1, three priors, with nodes in parallel on 4 cores (about 20 seconds; checked 2026-11-25 with the single template):
```
Rscript --vanilla run_task.R --n 1000 --p 50 --q 1 --rep 1 --ncores 4 --out results_local
```
