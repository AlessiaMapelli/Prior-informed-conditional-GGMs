# ggReg: Prior-Informed Conditional Gaussian Graphical Regression

**Associated publication:** *Incorporating prior knowledge into Conditional GGM: an application to Protein Interaction Networks* 

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

Hyperparameters (the sparsity-group balance `α` and the prior weight scaling `ϑ`) are tuned via BIC over a grid or random search. Model selection for `λ` uses cross-validation via `glmnet` and `sparsegl`.

---

## Repository Contents

| File | Description |
|---|---|
| `ggReg_main_functions.R` | Core estimation and utility functions |
| `ggReg_node_job.R` | Single-node estimation script for SLURM array jobs |
| `slurm_ggReg_node.sbatch` | SLURM batch script template for HPC execution |

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
  - `Cov_effect`: estimated covariate effects on the mean (p × q matrix)
  - `Dic_adj_matrics`: named list of partial-correlation adjacency matrices, one per covariate effect plus `"Baseline"`

- **`additional_info`**: diagnostic and intermediate outputs
  - `z`: residuals from mean estimation
  - `Dic_Delta_hat`: symmetrized precision matrices per covariate effect
  - `Sigma_hat`: estimated node-wise residual variances
  - `optimal_params`: best `(α, ϑ)` and BIC per node

---

### Predicting a Personalized Network

Once the model is estimated, generate a subject-specific network using `predict_personalized_network()`:

```r
new_subject <- data.frame(age = 52, sex = "Female", bmi = 27.4)

personal_net <- predict_personalized_network(
  Dic_Delta_hat        = results$additional_info$Dic_Delta_hat,
  training_covariates  = covariate_data,
  new_subject_covariates = new_subject
)
```

The function handles covariate scaling and dummy encoding automatically, matching the preprocessing applied during training.

---

### Visualization

Two utility functions are provided for network exploration:

```r
# Visualize a single personalized network
plot_personalized_network(
  network   = personal_net,
  labels    = protein_names,
  threshold = 0.05,
  main      = "Subject-specific PPI network"
)

# Compare two networks (e.g., two subjects or conditions)
plot_network_difference(network1 = net_A, network2 = net_B, labels = protein_names)

# Compute a summary comparison metric
compare_networks(net_A, net_B, method = "frobenius")  # also: "difference", "correlation"
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
| `scr` | `TRUE` | Correlation-based pre-screening to reduce problem size |
| `gamma` | `NULL` | Screening threshold; `NULL` uses the 20th percentile |
| `tune_hyperparams` | `TRUE` | Tune `α` and `ϑ` via BIC |
| `asparse_grid` | `c(0.5, 0.75, 0.9, 0.95)` | Grid for the group/element sparsity balance `α` |
| `weight_grid` | `c(0.8, 1.0, 1.1, 1.3, 1.5)` | Grid for the prior weight scaling `ϑ` |
| `random_hyper_search` | `FALSE` | Use random search instead of full grid for hyperparameter tuning |
| `symm_method` | `"OR"` | Symmetrization rule: `"OR"` (less conservative) or `"AND"` |
| `use_slurm` | `TRUE` | Use SLURM for HPC parallelization |

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
