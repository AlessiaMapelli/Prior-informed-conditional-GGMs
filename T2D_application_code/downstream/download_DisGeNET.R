#################################################
## Download the DisGeNET curated gene-disease associations for type 2 diabetes
## (UMLS C0011860) used in A2c / A4.3. Academic accounts can access CURATED sources only.
## The API key is read from the environment, never stored:
##   DISGENET_API_KEY=... Rscript --vanilla ../T2D_application_code/downstream/download_DisGeNET.R
## Run from the private project folder. Output: results_bic_relaxed_10y/downstream/external/
#################################################
suppressMessages(library(disgenet2r))
key <- Sys.getenv("DISGENET_API_KEY"); if (key == "") stop("Set DISGENET_API_KEY")
r <- disease2gene(disease = "UMLS_C0011860", database = "CURATED", api_key = key, verbose = FALSE, n_pags = 100)
d <- r@qresult[, c("gene_symbol", "score", "normalized_score", "evidence_index", "numPMIDs", "yearInitial", "yearFinal")]
out <- "results_bic_relaxed_10y/downstream/external"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
write.csv(d, file.path(out, "DisGeNET_T2D_C0011860_curated.csv"), row.names = FALSE)
cat("Saved", nrow(d), "associations (", format(Sys.Date()), ")\n")
