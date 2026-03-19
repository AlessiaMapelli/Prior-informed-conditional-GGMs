# ggReg: Prior-Informed Conditional Gaussian Graphical Models

This repository contains all code to reproduce the results of the paper:

> **Prior-informed conditional gaussian graphical models for personalized protein interaction network reconstruction**
>
> Alessia Mapelli<sup>1,2,∗</sup>, Michela Carlotta Massi<sup>2</sup>, Gianmauro Cuccuru<sup>2</sup>, Emanuele Di Angelantonio<sup>2</sup>, Francesca Ieva<sup>1,2</sup>
>
> <sup>1</sup> MOX, Department of Mathematics, Politecnico di Milano, Milan, Italy
> <sup>2</sup> Health Data Science Research Centre, Human Technopole, Milan, Italy

---

## Overview

**ggReg** is an R pipeline for network estimation from multivariate omics data. It is built around a node-wise penalized regression framework and supports a spectrum of models depending on what inputs are available:

- **Lasso GGM** — standard sparse network estimation with no covariates or prior, using only expression data.
- **Weighted GGM** — incorporates edge-specific confidence scores from biological databases (e.g., STRING): edges with strong prior support are penalized less, letting the data confirm or override prior information.
- **Conditional GGM** — network structure is allowed to vary across individuals as a function of external covariates (age, sex, clinical variables), producing a baseline population-level network plus covariate-specific perturbation terms.
- **Prior-informed Conditional GGM** (full model) — combines both extensions: prior knowledge informs the penalization while covariates drive individual-level variation, enabling truly *personalized* network reconstruction from a single shared model.

Given any new individual's covariates, their personalized precision matrix is recovered by linear combination of the estimated population-level and covariate-specific network components. Hyperparameters are tuned via BIC; the regularization path for `λ` is selected by cross-validation. For large-scale problems, node-wise jobs are parallelized on HPC clusters via SLURM array jobs.

---

## Repository Structure

```
.
├── Computational_templates/    # Core algorithm implementation and quick-start template
├── Simulation_analysis/        # Simulation study scripts
└── T2D_application_code/       # Real-data application to UKB T2D proteomics
```

---

### `Computational_templates/`

The self-contained implementation of the ggReg pipeline and the recommended starting point for new users. It contains:

- **`ggReg_main_functions.R`** — all core estimation, hyperparameter tuning, symmetrization, personalized network prediction, and visualization functions, accessible through the single entry point `GGReg_full_estimation()`.
- **`ggReg_node_job.R`** and **`slurm_ggReg_node.sbatch`** — scripts for parallelizing node-wise estimation on HPC clusters via SLURM array jobs.
- **Example data** (`Example_data.RData`, `expression_data.csv`, `covariate_data.csv`, `ppi_weight_matrix.txt`) — a small simulated dataset (1000 subjects, 10 proteins, 1 binary covariate, perfect prior) to test and demonstrate the pipeline locally before scaling up.

See `Computational_templates/README.md` for full documentation, including installation, usage examples, a parameter reference, and visualization utilities.

---

### `Simulation_analysis/`

Scripts to reproduce the simulation study evaluating ggReg's performance across controlled settings. The study crosses sample sizes (500, 1000), number of nodes (50, 100), prior quality (perfect, noisy, absent), and number of covariates (1 or 3), with 10 replications per configuration. Networks are evaluated on FPR, F1, accuracy, and edge-weight preservation.

Key scripts:

- **`Sim_functions.R`** — all helper functions for network generation (scale-free and Erdős-Rényi topologies), data simulation, prior construction, and result evaluation.
- **`Sim_pipeline_run.R`** — the main orchestration script, covering a test run, sequential execution, HPC parallel execution (via `generate_input_datasets_simulation()` + SLURM), result collection, and figure production for the paper.
- **`Simulation_sequential.sh`** — Bash script that submits one SLURM array job per simulation configuration, reading the job list produced by `generate_input_datasets_simulation()`, with a cap of 100 concurrent jobs.
- **`ggReg_main_functions.R`**, **`ggReg_node_job.R`**, **`slurm_ggReg_node.sbatch`** — identical copies of the files in `Computational_templates/`.
- **`Example_data/Simulate_example_data.R`** — generates the example dataset included in `Computational_templates/` using the same simulation infrastructure but with a minimal configuration.

See `Simulation_analysis/README.md` for details.

---

### `T2D_application_code/`

Scripts to reproduce the real-data application of ggReg to UK Biobank (UKB) plasma proteomics for Type 2 Diabetes (T2D) network reconstruction. The analysis runs in three sequential steps:

1. **`Network_estimation_case_study_T2D_UKB.R`** — split the UKB cohort (366 proteins, covariates: Sex, Age, T2D status) into a network-estimation set, a matched case-control training set, and a test set, and runs `GGReg_full_estimation()` with SLURM parallelization including a preprocessed STRING based prior weight matrix. Results are saved to `./results/PPI_T2D_app.RData`.

2. **`Network_description_case_study_T2D_UKB.Rmd`** — a parameterized R Markdown report performing all downstream network analysis: mean effect extraction, network visualization, community detection, node centrality analysis, identification of top influential nodes, differential protein expression, and logistic risk prediction evaluated on the held-out test set. All figures and tables are saved to `./results/Network_description/`.

3. **`enrichment_communities.py`** — functional enrichment analysis of detected communities with results annotated using `kegg_hierarchy.csv` and summarized as a KEGG subcategory × community heatmap.

See `T2D_application_code/README.md` for details.

---

## Citation

If you use this code, please cite:

> Mapelli, A. *Prior-informed conditional gaussian graphical models for personalized protein interaction network reconstruction*. Draft, Politecnico di Milano.
