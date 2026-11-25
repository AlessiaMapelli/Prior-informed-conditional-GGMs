#################################################
## Configuration shared by all scripts of the T2D application
#################################################
# The scripts are run from the PRIVATE project folder (which holds the UK
# Biobank data and receives the results); the code is read from this folder.
# Environment variables:
#   T2D_CODE_DIR     this folder (default ../T2D_application_code, relative to
#                    the private folder)
#   T2D_INPUT_RDATA  input data built by preprocessing/ (default
#                    T2D_app_input.RData in the private folder)
#   T2D_RES          results folder (default results_bic_relaxed_10y)
#   PRCGGM_LIB       optional R library holding prCGGM >= 0.2.0
#   STRING_DIR       folder with the STRING v12.0 files (protein.info,
#                    protein.physical.links) used by preprocessing/ and A4d
CODE_DIR <- Sys.getenv("T2D_CODE_DIR", "../T2D_application_code")
INPUT_RDATA <- Sys.getenv("T2D_INPUT_RDATA", "T2D_app_input.RData")
STRING_DIR <- Sys.getenv("STRING_DIR", "external/STRING")

#' Load prCGGM (>= 0.2.0: selection = "bic_relaxed", symm_method = "AND")
#' from PRCGGM_LIB if set, otherwise from the default library
use_prcggm <- function() {
  lib <- Sys.getenv("PRCGGM_LIB")
  if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(prCGGM))
  if (!"bic_relaxed" %in% eval(formals(prCGGM::prcggm)$selection)) {
    stop("prCGGM ", packageVersion("prCGGM"), " has no selection = \"bic_relaxed\": install prCGGM >= 0.2.0")
  }
  invisible(packageVersion("prCGGM"))
}
