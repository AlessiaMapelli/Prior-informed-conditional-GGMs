# prCGGM: prior-informed conditional Gaussian graphical models

This repository contains the code to reproduce the results of the paper:

> **Prior-informed conditional Gaussian graphical models: an application to protein interaction network reconstruction**
>
> Alessia Mapelli<sup>1,2,∗</sup>, Michela Carlotta Massi<sup>2</sup>, Gianmauro Cuccuru<sup>2</sup>, Emanuele Di Angelantonio<sup>2</sup>, Francesca Ieva<sup>1,2</sup>
>
> <sup>1</sup> MOX, Department of Mathematics, Politecnico di Milano, Milan, Italy
> <sup>2</sup> Health Data Science Research Centre, Human Technopole, Milan, Italy
>
> arXiv:[2606.31805](https://arxiv.org/abs/2606.31805)

---

## The method

A prior-informed conditional GGM estimates a Gaussian graphical model whose mean and precision matrix depend linearly on subject-level covariates:

Ω(X) = Δ₀ + Σₕ Xₕ Δₕ

- **Δ₀** is the population-level (baseline) network.
- **Δₕ** is the change of the network associated with covariate h (for example sex, age or disease status).
- **Prior knowledge** on the edges, for example STRING confidence scores, enters through a weighted sparse-group lasso penalty that is applied **only to Δ₀**. The covariate-specific changes stay entirely data-driven.
- **Personalised networks** of new individuals are obtained from their covariates.

The model is fitted by node-wise regressions:
- **Penalty and hyperparameters:** chosen per node by the BIC of the least-squares refit (relaxed BIC).
- **Edges:** the selected edges are refitted by least squares, and an edge is kept only if both of its node-wise regressions select it (AND).
- **Special cases:** without covariates or without a prior, the model reduces to a lasso GGM, a prior-weighted GGM or a conditional GGM.

## Installation

The method is implemented in the R package [`prCGGM`](https://CRAN.R-project.org/package=prCGGM). It is the recommended way to use the method on new data:

```r
install.packages("prCGGM")   # version >= 0.2.0
library(prCGGM)
vignette("prCGGM")
```

Use **prCGGM ≥ 0.2.0**. It has the selection rules, refit and defaults used in the paper; version 0.1.0 used cross-validation with `lambda.min`. The package entry point is `prcggm()`, and `predict()` gives the personalised networks. `Computational_templates/prCGGM_example.R` runs the example of the template with the package.

---

## Repository structure

```
.
├── Computational_templates/        # standalone estimator (single copy) and example
├── Simulation_analysis_v2/         # simulation study of the paper (main text and Supplementary)
│   └── selection_rules/            # Supplementary Part 1: comparison of selection rules
├── Add_simulation_analysis/        # exploratory simulations that led to the design choices
├── T2D_application_code/           # UK Biobank application (code only, no data)
├── Preprint_simulation_analysis/   # code of the preprint's simulation (superseded)
├── Preprint_T2D_application_code/  # code of the preprint's application (superseded)
└── update_in_dev_version.md        # log of the changes since the preprint
```

### `Computational_templates/`

The standalone R implementation, **`ggReg_main_functions.R`**, with entry point `GGReg_full_estimation()`.
- **Single copy:** it is the only copy in the repository, and every simulation script sources it.
- **Agreement with the package:** it implements the same estimator as prCGGM ≥ 0.2.0. `check_template_vs_package.R` checks that the two give the same networks.
- **Template-only options:** it also has options used only in the Supplementary of the paper, which are not in the package:
  - the post-selection filters `post_filter = "pcor"` and `"backbone"`;
  - the ablation option `prior_on_covariates`.
- **Other files:**
  - SLURM scripts for node-wise array jobs;
  - a small simulated example dataset (1000 subjects, 10 proteins, 1 binary covariate);
  - `prCGGM_example.R`.

See `Computational_templates/README.md`, including the section "Template and package".

### `Simulation_analysis_v2/`

The simulation study of the paper. The design crosses:

| Factor | Levels |
|---|---|
| Sample size n | 1000, 5000, 10,000, 20,000 |
| Number of proteins p | 50, 100, 250, 500 |
| Number of covariates q | 1, 3 |
| Prior | perfect, noisy, none |
| Replications | 100 per setting |

Networks are evaluated per component (baseline, binary effect, continuous effect) on F1, TPR, FPR and edge weights.

The folder also holds:
- the comparison of selection rules (`selection_rules/`, Part 1, and `selection_test/`, Part 2);
- the reruns at n = 1000 with the extended BIC and stability selection;
- the reproducibility experiment;
- the ablation of the prior on the covariate terms;
- the paper figures (`paper_figures.R`, `fig_ablation_main.R`).

**`README_run_simulation.md`** maps every simulation figure and table of the paper to its cluster job, summary script and output, with the seeds and the software. `README.md` records the design, the decisions and the results.

### `Add_simulation_analysis/`

Exploratory simulations, not reported in the paper. They showed why the cross-validated lasso of the preprint selects very dense networks when n is large, and they led to the relaxed BIC with refit. Its README documents these analyses.

### `T2D_application_code/`

The real-data application: 366 cardiometabolic Olink proteins of UK Biobank participants without diabetes at baseline, with sex, age and 10-year incident type 2 diabetes as covariates and a STRING v12.0 prior.

It covers:
- the preprocessing (data fields, phenotypes, cohort, STRING prior);
- the network fit with prCGGM (SLURM array);
- the permutation null of the T2D-specific network;
- the downstream analyses (topology, enrichment, replication in held-out data, external interaction databases);
- the prediction of incident T2D.

The analysis uses UK Biobank data under application 102297. **UK Biobank data cannot be redistributed, so the folder contains code only.** The scripts are run from a private folder that holds the data, with paths set by environment variables. See `T2D_application_code/README.md`.

### `Preprint_*`

The code of the preprint's version of the paper, kept for the record:
- `Preprint_simulation_analysis/`: the simulation of the preprint;
- `Preprint_T2D_application_code/`: the application of the preprint (5-year endpoint, cross-validated selection).

The revised analyses above supersede both.

---

## Paper items

| Paper item | Folder | Details |
|---|---|---|
| Simulation study, main text and Supplementary | `Simulation_analysis_v2/` | `README_run_simulation.md` |
| Selection rules (Supplementary Parts 1 and 2) | `Simulation_analysis_v2/selection_rules/`, `Simulation_analysis_v2/selection_test/` | `README_run_simulation.md` |
| T2D application: network, permutation null, downstream analyses, prediction | `T2D_application_code/` | `T2D_application_code/README.md`, "Paper items" |
| Software | prCGGM ≥ 0.2.0 (CRAN); `Computational_templates/` | `Computational_templates/README.md` |

---

## Citation

If you use this code or the package, please cite:

> Mapelli, A., Massi, M. C., Cuccuru, G., Di Angelantonio, E. and Ieva, F. (2026). Prior-informed conditional Gaussian graphical models: an application to protein interaction network reconstruction. arXiv:2606.31805. doi:[10.48550/arXiv.2606.31805](https://doi.org/10.48550/arXiv.2606.31805)

In R: `citation("prCGGM")`.
