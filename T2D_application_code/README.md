# T2D application: UK Biobank proteomics and incident type 2 diabetes

This folder holds the code of the paper's real-data application. A prior-informed conditional GGM is fitted with the R package prCGGM (≥ 0.2.0) to 366 cardiometabolic Olink proteins (UK Biobank Pharma Proteomics Project) of participants without diabetes at baseline. The covariates are sex, age and 10-year incident type 2 diabetes (T2D), and the prior comes from STRING v12.0.

**Data.** The analysis uses UK Biobank data under application **102297**. UK Biobank data cannot be redistributed, so **this folder contains code only**.
- Researchers can apply for access at <https://www.ukbiobank.ac.uk>.
- The input data, the participant IDs of the splits and all results stay in a private project folder.
- The scripts are run from that private folder and read the code from here.
- `.gitignore` blocks data files (`*.RData`, `*.rds`, `*.csv`, `*.tsv`, `*.txt`, `results*/`, ...) from being committed here.

The code that produced the preprint's version of this application is in `../Preprint_T2D_application_code/`.

## Requirements

- **R ≥ 4.3** and **prCGGM ≥ 0.2.0** (`install.packages("prCGGM")`; until 0.2.0 is on CRAN, install it from the `v.0.2` branch of the package repository).
- R packages:
  - modelling: `glmnet`, `sparsegl`, `Matrix`, `pROC`, `caret`;
  - data handling and figures: `dplyr`, `tidyr`, `ggplot2`, `ggraph`, `ggrepel`, `cowplot`, `RColorBrewer`, `pheatmap`, `knitr`, `kableExtra`;
  - networks and enrichment: `igraph`, `msigdbr`;
  - DisGeNET download only: `disgenet2r`.
- A SLURM cluster for the network fit (one array task per protein, 366 tasks) and for the permutation null (101 tasks). The downstream analyses run on a workstation.

## Set-up

Lay the folders out like this:
```
<project>/
├── T2D_application_code/   # this folder (public)
└── <private folder>/       # UK Biobank data and results (never published)
```
Then run every command from the private folder, with these environment variables:

| Variable | Meaning | Default |
|---|---|---|
| `T2D_CODE_DIR` | this folder | `../T2D_application_code` |
| `T2D_INPUT_RDATA` | input built by `preprocessing/` | `T2D_app_input.RData` |
| `T2D_ENDPOINT` | `10y` (primary) or `5y` (sensitivity analysis) | `10y` |
| `T2D_RES` | results folder | `results_bic_relaxed_10y` |
| `PRCGGM_LIB` | optional R library holding prCGGM ≥ 0.2.0 | — |
| `STRING_DIR` | folder with the STRING v12.0 files | `external/STRING` |
| `R_SETUP` | optional script sourced by the `.sbatch` files to set up R on the cluster | — |
| `UKB_WORK_DIR`, `UKB_RAW_DIR`, `UKB_PPP_DIR`, `UKB_APP_ID`, `UKB_CPR_FILE` | preprocessing only: locations of the extracted UK Biobank tables and of the UKB-PPP proteomics files, and the application ID used in their file names | see `preprocessing/` |

## Steps

### 0. Data preparation (`preprocessing/`)

1. **`extract_fields/*.yml`** list the UK Biobank data fields extracted (participant characteristics, biomarkers including HbA1c, medication and primary-care records, follow-up). They are configuration files of the extraction tool used at the authors' institute. The field IDs document the variables used; with other tools, extract the same fields.
2. **`Phenotype_preprocessing.R`** derives the phenotypes:
   - diabetes status at baseline and incident T2D:
     - incident T2D is a T2D diagnosis after the baseline visit;
     - prevalent T2D is excluded, and so are participants flagged as undiagnosed diabetes (baseline HbA1c > 48 mmol/mol) or treated diabetes (antidiabetic prescription at baseline);
   - dates of diagnosis, death and loss to follow-up;
   - covariates;
   - the QDiabetes-2018 variables and the polygenic risk score.
   - Writes `full_pheno_data_*.csv` and `QDiab2018_data_*.csv`.
3. **`Case_specific_data_preprocessing_T2D_app_case.R`** builds the cohort and saves the input data (`T2D_INPUT_RDATA`):
   - participants without diabetes at baseline (type 1 diabetes excluded);
   - 10-year and 5-year endpoints;
   - follow-up time censored at death or loss to follow-up;
   - participants with a diagnosis or end of follow-up within 30 days of baseline excluded;
   - merge with the Olink proteomics;
   - the 366 cardiometabolic proteins.
4. **`STRING_preprocessing.R`** builds the prior W from STRING v12.0 (`9606.protein.links.v12.0.txt` and `9606.protein.info.v12.0.txt` from <https://string-db.org>). The combined scores are rescaled to [0, 1]; pairs without a STRING link get 0. Writes `STRING_W`, `STRING_W_rownames.txt` and `STRING_W_colnames.txt`.

### 1. Network estimation

`Network_estimation_case_study_T2D_UKB.R` (with `helper_functions.R` and `slurm_prcggm_node.sbatch`).

- **Splits:**
  - selection set: 70%, stratified by T2D × sex;
  - held-out set: 30%, used for replication and prediction.
- **Model:**
  - covariates sex, age (centred) and incident T2D;
  - p = 366, so P = 1,460 predictors per node;
  - `selection = "bic_relaxed"`, `symm_method = "AND"`, OLS refit, full 4 × 4 grid of (α, ϑ).

The script runs in two passes:
```
Rscript --vanilla ../T2D_application_code/Network_estimation_case_study_T2D_UKB.R   # 1: split, mean step, prints the sbatch command
sbatch --array=1-366 ../T2D_application_code/slurm_prcggm_node.sbatch <input.rds> <node dir>
Rscript --vanilla ../T2D_application_code/Network_estimation_case_study_T2D_UKB.R   # 2: collects the fit -> $T2D_RES/PPI_T2D_app.RData
```
For the 5-year sensitivity analysis set `T2D_ENDPOINT=5y`; results go to `results_bic_relaxed/`.

### 2. Permutation null of the T2D-specific network

`T2D_permutation_null.R`, `slurm_permutation_null.sbatch`.

- **What is permuted:** the incident-T2D labels, within sex × age-quintile strata of the selection set.
- **What is refitted:** the mean step and every node, at the hyperparameters selected on the real data, followed by AND / OR symmetrisation.
- **Size:** 100 permutations, plus permutation 0 with the real labels, which checks that the procedure reproduces the fit.
```
Rscript --vanilla ../T2D_application_code/T2D_permutation_null.R prepare 100
sbatch --array=1-101%50 ../T2D_application_code/slurm_permutation_null.sbatch
Rscript --vanilla ../T2D_application_code/T2D_permutation_null.R collect
Rscript --vanilla ../T2D_application_code/T2D_permutation_hub_degree.R   # hub test
Rscript --vanilla ../T2D_application_code/T2D_permutation_hub_fdr.R      # FDR of the hub classes
Rscript --vanilla ../T2D_application_code/T2D_permutation_maxT_edges.R   # edge-level maxT test
```

### 3. Downstream analysis (`downstream/`)

All scripts source `downstream/00_common.R` and write to `$T2D_RES/downstream/`. They are run in this order (`Rscript --vanilla ../T2D_application_code/downstream/<script>.R`):

| Script | Content |
|---|---|
| `A1_global_null.R` | Global test: number of T2D-specific edges against the permutation null |
| `A2_topology_DE.R` | Topology of the T2D OR graph (degree against the null, communities, centrality, hub classes); differential expression on the selection set |
| `A2b_enrichment.R` | Over-representation analysis of communities and protein sets (msigdbr: KEGG Legacy, Reactome, GO BP; background: the 366 panel proteins) |
| `download_DisGeNET.R` | DisGeNET curated T2D genes (UMLS C0011860). The API key is read from `DISGENET_API_KEY` and never stored |
| `A2c_A4c_disease_context.R` | DisGeNET T2D genes and diabetes pathways in the T2D graph |
| `A3_A4_TIMP1_replication_STRING.R` | TIMP1 neighbourhood; replication of the T2D edges in the held-out data; STRING support against the null |
| `A4b_replication_permutation.R` | Permutation calibration of the held-out replication |
| `A4d_STRING_physical.R` | STRING v12.5 physical interactions against the null (`9606.protein.physical.links.v12.5.txt.gz` in `$T2D_RES/downstream/external/`) |
| `A4e_HPA_interactions.R` | Human Protein Atlas interaction consensus against the null (`interaction_consensus.tsv` in `external/`) |
| `A5_edge_table.R` | Descriptive table of all T2D OR edges |
| `B_prediction.R` | Prediction of 10-year incident T2D in the held-out data: protein sets from the network against size-matched differential-expression and lasso sets, alone and added to age + sex, the QDiabetes-2018 variables and QDiabetes + HbA1c; DeLong and held-out likelihood-ratio tests |
| `B_paper_table_figure.R` | Prediction table and figure of the paper (`B_table_paper.tex`, `figures/Fig4b_prediction_by_base`) |
| `Figures.R` | Main and supplementary figures (`figures/`) |

The external reference files are downloaded by the user into `$T2D_RES/downstream/external/`:
- STRING v12.5 physical links;
- the HPA interaction consensus;
- the DisGeNET export;
- the protein list of García-Currás et al. (2026), Supplementary Table S6.

### Paper items

| Paper item | Script → file |
|---|---|
| T2D evidence against the permutation null (main) | `Figures.R` → `Fig1_T2D_validation` |
| T2D network, hubs against differential expression, TIMP1 (main) | `Figures.R` → `Fig2_T2D_network`, `Fig3_hubs_DE_TIMP1` |
| Prediction with likelihood-ratio tests (main) | `B_paper_table_figure.R` → `Fig4b_prediction_by_base`; table `B_table_paper.tex` (Supplementary) |
| Permutation nulls, enrichment, hub test, external references (Supplementary) | `Figures.R` → `FigS1`–`FigS4`; tables from A1–A5 |
| 5-year sensitivity analysis (Supplementary) | steps 1–2 with `T2D_ENDPOINT=5y`, `T2D_RES=results_bic_relaxed`; `A1_global_null.R` |
