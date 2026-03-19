# T2D Application Code

This folder contains all scripts needed to reproduce the real-data application of the prior-informed covariate-dependent Gaussian Graphical Model (GGM) to the UK Biobank (UKB) proteomic dataset for Type 2 Diabetes (T2D) research.

The analysis is organized in three sequential steps.

---

## File Overview

### `ggReg_main_functions.R`, `ggReg_node_job.R`, `slurm_ggReg_node.sbatch`

These three files are identical copies from the `Computational_templates` folder. They implement the core algorithm and its HPC execution logic. Refer to that folder for documentation.

---

## Step 1 — Network Estimation: `Network_estimation_case_study_T2D_UKB.R`

The main R script that prepares the data and estimates the prior-informed covariate-dependent protein interaction network.

**Cohort preparation.** The script loads the UKB proteomic dataset (`UKB_proteo_T2D.RData`), a data frame of 366 protein expression measurements alongside three covariates: Sex, Age, and T2D status (`Diab`). The cohort is split into three non-overlapping sets:

- **Selection set (80%)** — used for network estimation, stratified by T2D status.
- **Training set** — a matched case-control subset drawn from the remaining 20%, used for downstream risk prediction. Incident T2D cases are matched 1:4 to non-diabetic controls on Sex and Age.
- **Test set** — the remaining participants from the 20% holdout, used to evaluate risk prediction.


**Prior knowledge.** A prior preprocessed weight matrix from the STRING protein–protein interaction database is loaded, providing edge confidence scores that inform the network estimation.

**Network estimation.** `GGReg_full_estimation()` is called on the selection set with Sex and T2D status as covariates. The model jointly estimates:
- Mean effects of each covariate on each protein (via penalized regression).
- A covariate-dependent precision matrix encoding the conditional independence structure of the protein network.

Hyperparameters (sparsity `alpha`, prior weight `w`) are tuned via 5-fold cross-validation over a random search grid. Node-level computations are parallelized on the HPC cluster via SLURM (`slurm_ggReg_node.sbatch`).

Results are saved to `./results/PPI_T2D_app.RData`, together with the three cohort data frames.

---

## Step 2 — Network Description: `Network_description_case_study_T2D_UKB.Rmd`

A parameterized R Markdown report that performs all downstream network analysis. It is rendered for the T2D incidence outcome (`DiabIncident`) and writes all outputs to `./results/Network_description/DiabIncident/`. The report covers the following sections:

1. **Mean effects analysis.** Extracts proteins with a non-zero T2D mean effect from the estimated coefficient matrix and reports their R-squared values across covariates.

2. **Network construction and visualization.** Builds an `igraph` object from the estimated precision matrix and reports basic network statistics (number of nodes and edges, edge density, diameter). Produces two network plots: one showing the basic structure, one highlighting proteins influenced by T2D in the mean. A publication-quality figure using `ggraph` with edge weights encoded by width, colour, and transparency is also saved.

3. **Community detection.** Applies the fast-greedy algorithm (`cluster_fast_greedy`) to identify protein communities. Community membership is visualized with a colour-coded network plot.

4. **Node centrality analysis.** Computes degree, eigenvector centrality, PageRank, betweenness, and hub score for each node. Scores are normalized and averaged into a composite centrality measure. The top-10% nodes by both average centrality and hub score are identified; their intersection defines the set of *top influential nodes*, which are visualized as a labelled subgraph. The most important interactions (edges in the top 1% by weight within this subgraph) are also extracted.

5. **Node topology export.** The full per-node summary table (community membership, all centrality measures, important interactions) is saved as `Summary_node_topology_DiabIncident.csv`. This file is the direct input for Step 3.

6. **Differential expression analysis.** For each protein, a Wilcoxon rank-sum test compares expression between T2D groups. P-values are adjusted with the Benjamini–Hochberg procedure and results are displayed as a volcano plot.

7. **Risk prediction.** A logistic regression model is trained on the matched training set using the top influential nodes as predictors, and evaluated on the test set. Performance is reported as AUC with 95% confidence intervals, and ROC curves for both sets are plotted.

---

## Step 3 — Enrichment Analysis: `enrichment_communities.py`

A Python script that performs functional enrichment analysis on the protein communities identified in Step 2 and produces a summary heatmap.

**Input.** Reads `Summary_node_topology_DiabIncident.csv` (produced by Step 2) and groups proteins by their community assignment.

**Enrichment.** Runs over-representation analysis via `gseapy.enrichr` against three databases:
- KEGG 2021 Human
- GO Biological Process 2023
- Human Phenotype Ontology

Results are filtered at adjusted p-value ≤ 0.05 and saved as CSV files.

**KEGG hierarchy annotation.** KEGG terms are mapped to higher-level subcategories using `kegg_hierarchy.csv` (included in this folder), allowing biologically related pathways to be grouped for visualization. Related subcategories (e.g. different types of infectious disease, different branches of metabolism) are recoded into broader groups.

**Output.** A heatmap (`kegg_subcategory_by_community_heatmap.png`) showing the number of enriched KEGG subcategories per community (capped at 5 for visual clarity), saved to `./results/Network_description/`.
