# ggReg: Prior-Informed Conditional Gaussian Graphical Regression

**Associated publication:** *Prior-Informed Conditional Gaussian Graphical Models for Personalized Protein Interaction Network Reconstruction* 

> **The method is now available as the R package [`prCGGM`](https://CRAN.R-project.org/package=prCGGM), which is the recommended way to use it** (`install.packages("prCGGM")`; see `prCGGM_example.R` for the example below run with the package). In the package, `GGReg_full_estimation()` becomes `prcggm()` and `predict_personalized_network()` becomes `predict()`. On the example data, both implementations give the same networks. This folder documents the standalone code associated with the paper.

---

## Overview

`ggReg` is an R pipeline for estimating prior-infused conditional gaussian graphical regression models — particularly personalized biological networks from omics — in high-dimensional settings. It extends the standard Gaussian graphical model (GGM) framework in two key ways:

1. **Prior knowledge integration**: confidence scores from biological databases (e.g., STRING) are used to penalize edges with low prior support less heavily, while standard edges are penalized more.
2. **Covariate-conditional estimation**: network structure is allowed to vary across individuals as a function of external covariates (e.g., age, sex, clinical variables), enabling truly *personalized* network reconstruction.

The result is a subject-specific partial correlation network — and associated precision matrix — that reflects both population-level biology and individual characteristics. These personalized networks can be used directly for downstream analyses such as biomarker discovery, differential network analysis, and disease prediction.

---

## Method Summary

The estimation follows a three-step procedure:

- **Step 1 – Mean estimation**: a node-wise Lasso regression removes the effect of covariates on the mean expression of each feature, producing residuals `Z`.
- **Step 2 – Precision matrix estimation**: node-wise penalized regressions on `Z` estimate how network edges vary with covariates. The penalty incorporates prior knowledge through a weighted Lasso: edges with strong prior support (high `W_ij`) receive reduced penalization. A sparse-group Lasso structure separates population-level edge effects from covariate-modulated effects.
- **Step 3 – Personalized network prediction**: given a new individual's covariate vector, their personalized precision matrix is reconstructed by linearly combining the estimated population-level and covariate-specific network components.

For each node, the penalty `λ` and the hyperparameters (the sparsity-group balance `α` and the covariate penalty scaling `ϑ`) are selected by the rule in `selection`:

| `selection` | `λ` | `α`, `ϑ` | Cost |
|---|---|---|---|
| `"ebic_relaxed"` (default) | minimum extended BIC (γ = 0.5) along the `λ` path, computed on the RSS of the OLS refit of each support | jointly with `λ` | `λ` paths + one Gram matrix per node |
| `"bic_relaxed"` | as above with the BIC | jointly with `λ` | as above |
| `"bic"` / `"ebic"` | minimum BIC / extended BIC on the RSS of the penalised fit (over-select as n grows; kept for comparison) | jointly with `λ` | `λ` paths only |
| `"cv_min"` / `"cv_1se"` | K-fold cross-validation, `lambda.min` / `lambda.1se` | BIC at the selected `λ` | paths + CV |
| `"stability_bic"` | stability selection (50 half-subsamples, selection probability ≥ 0.75, ≤ 1 expected false selection per node) | those of `"bic_relaxed"` | paths + 50 subsample paths |
| `"stability_cv"` | as above | those of `"cv_min"` | paths + CV + 50 subsample paths |

With `refit = TRUE` (default) the selected predictors of each node are refitted by OLS, which removes the shrinkage of the edge weights; the selected edges are unchanged. The default (`"ebic_relaxed"` with refit) and the alternatives were compared in simulation studies (see `Add_simulation_analysis/` and `Simulation_analysis_v2/selection_test/`): cross-validation with `lambda.min`, the rule of the first version of this template, and BIC on the penalised RSS both select increasingly dense networks as the sample size grows; the relaxed criteria do not, and `"ebic_relaxed"` was as accurate as stability selection at about one fifth of its cost.

---

## Repository Contents

| File | Description |
|---|---|
| `ggReg_main_functions.R` | Core estimation and utility functions |
| `ggReg_node_job.R` | Single-node estimation script for SLURM array jobs |
| `slurm_ggReg_node.sbatch` | SLURM batch script template for HPC execution |
| `prCGGM_example.R` | The same example run with the `prCGGM` R package |

---

## Quick Start

### Installation

Install the required R packages:

```r
install.packages(c("glmnet", "sparsegl", "doParallel", "foreach"))
# Optional, for network visualization:
install.packages("corrplot")
```

Then source the main functions file:

```r
source("ggReg_main_functions.R")
```

---

### Basic Usage

Test data can be loaded via: 
```r
load("Example_data.RData")
```

The entire pipeline is accessible through a single entry-point function: `GGReg_full_estimation()`.

**Minimal call — standard GGM, no covariates or prior:**
```r
results <- GGReg_full_estimation(
  x = expression_data  # data.frame: subjects × features
)
```

**With covariates and prior knowledge (recommended):**
```r
results <- GGReg_full_estimation(
  x           = expression_data,   # data.frame: subjects × features
  known_ppi   = ppi_weight_matrix, # matrix [0,1]: features × features
  covariates  = covariate_data,    # data.frame: subjects × covariates
  use_slurm   = FALSE              # set TRUE on HPC (see below)
)
```

The function automatically selects the computation strategy based on the number of features `p`:
- `p ≤ 150`: sequential
- `p > 150`: R parallel (`foreach`/`doParallel`)
- `use_slurm = TRUE`: SLURM array jobs (one job per node)

---

### Output

`GGReg_full_estimation()` returns a list with two components:

- **`results`**: the primary outputs for downstream use
  - `Cov_effect`: estimated covariate effects on the mean (p × (q+1) matrix)
  - `Dic_adj_matrics`: named list of partial-correlation adjacency matrices, one per covariate effect plus `"Baseline"`

- **`additional_info`**: diagnostic and intermediate outputs. Of note:
  - `Dic_Delta_hat`: symmetrized precision matrices per covariate effect
  - `Sigma_hat`: estimated node-wise residual variances
  - `optimal_params`: per node, the selected `λ`, `(α, ϑ)`, the BIC of the reported model, the number of selected predictors, the selection rule and whether the coefficients were refitted
  - `scaling_params`: scale parameters for the continuos varaibles
  - `dummy_params`: dummy coding for the categorical variables

---

### Predicting a Personalized Network

Once the model is estimated, generate a subject-specific network using `predict_personalized_network()`:

```r
new_subjects <- data.frame(x1=c(0,1))
new_subjects$x1 <- as.factor(new_subjects$x1)

personal_net <- predict_personalized_network(
  Dic_adj_matrics        = results$results$Dic_adj_matrics,
  new_subject_covariates = new_subjects,
  scaling_params = results$additional_info$scaling_params,
  dummy_params = results$additional_info$dummy_params
)
```

The function handles covariate scaling and dummy encoding automatically, matching the preprocessing applied during training and return the personalized adjacency matrix for each subject included in the new_subjects data.frame.

---

### Visualization

Two utility functions are provided for network exploration:

```r
# Visualize a single network component
plot_personalized_network(
  network   = results$results$Dic_adj_matrics$Baseline,
)

# Visualize a single personalized network
plot_personalized_network(
  network   = personal_net$Subject_1,
)

# Compare two networks (e.g., two subjects or conditions)
plot_network_difference(network1 = personal_net$Subject_1, network2 = personal_net$Subject_2)

# Compute a summary comparison metric
compare_networks(network1 = personal_net$Subject_1, network2 = personal_net$Subject_2, method = "frobenius")  # also: "difference", "correlation"
```

---

### HPC Execution (SLURM)

For large-scale problems (`p > 1000`), node-wise estimation is reccomended to be parallelized via SLURM array jobs. Before running, edit `slurm_ggReg_node.sbatch` to set your email, partition name, and any cluster-specific module loading. Then call `GGReg_full_estimation()` with:

```r
results <- GGReg_full_estimation(
  x                 = expression_data,
  known_ppi         = ppi_weight_matrix,
  covariates        = covariate_data,
  use_slurm         = TRUE,
  slurm_script_path = "./slurm_ggReg_node.sbatch",
  output_path       = "./results/",
  name_output       = "my_run"
)
```

Each SLURM array task runs `ggReg_node_job.R` independently for a single node, and results are automatically collected and aggregated once all jobs complete.

---

## Key Parameters Reference

| Parameter | Default | Description |
|---|---|---|
| `known_ppi` | `NULL` | Prior confidence matrix `[0,1]`; `NULL` disables prior |
| `covariates` | `NULL` | Covariate data frame; `NULL` estimates a standard GGM |
| `scr` | `FALSE` | Pre-screening to reduce problem size (correlation-based or via `scr_matrix`) |
| `gamma` | `NULL` | Correlation screening threshold; `NULL` uses the 20th percentile. Ignored if `scr_matrix` is provided |
| `scr_matrix` | `NULL` | Custom screening matrix `{0,1}`; entry `(j,k)=0` excludes feature `k` from node `j`'s regression. Overrides correlation-based screening when provided. Prior-supported edges (`known_ppi > 0`) are never screened out regardless |
| `tune_hyperparams` | `TRUE` | Tune `α` and `ϑ` via BIC. If `FALSE` take `asparse_grid[1]` and `weight_grid[1]` as parameters|
| `asparse_grid` | `c(0.75, 0.5, 0.9, 0.99)` | Grid for the group/element sparsity balance `α` |
| `weight_grid` | `c(1.0, 0.8, 1.1, 1.3)` | Grid for the covariate penalty scaling `ϑ` |
| `random_hyper_search` | `FALSE` | Use random search instead of full grid for hyperparameter tuning |
| `selection` | `"ebic_relaxed"` | Rule selecting `λ` (and `α`, `ϑ`): `"ebic_relaxed"`, `"bic_relaxed"`, `"bic"`, `"ebic"`, `"cv_min"`, `"cv_1se"`, `"stability_bic"`, `"stability_cv"` (see Method Summary) |
| `refit` | `TRUE` | OLS refit of the selected predictors (always `TRUE` for stability selection) |
| `selection_control` | `list()` | Settings of the rules: `ebic_gamma` (0.5), `stab_B` (50), `stab_fraction` (0.5), `stab_pi` (0.75), `stab_pfer` (1) |
| `K` | `5` | Folds of the cross-validation rules |
| `lambda_prec` | `NULL` | Fixed `λ` for all nodes (the hyperparameters are then chosen by BIC) |
| `symm_method` | `"OR"` | Symmetrization rule: `"OR"` (less conservative) or `"AND"` |
| `use_slurm` | `FALSE` | Use SLURM for HPC parallelization |

---

## Citation

If you use this code, please cite:

> Mapelli, A. *Incorporating prior knowledge into Conditional GGM: an application to Protein Interaction Networks*. Draft, Politecnico di Milano.

---

## Dependencies

- [`glmnet`](https://cran.r-project.org/package=glmnet) — Lasso mean estimation
- [`sparsegl`](https://cran.r-project.org/package=sparsegl) — Sparse group Lasso precision estimation
- [`doParallel`](https://cran.r-project.org/package=doParallel) / [`foreach`](https://cran.r-project.org/package=foreach) — R-based parallelization
- [`corrplot`](https://cran.r-project.org/package=corrplot) *(optional)* — Network visualization
