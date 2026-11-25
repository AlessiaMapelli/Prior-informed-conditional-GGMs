######################################################
## STRING prior for the T2D application
######################################################
# Builds the prior matrix W (proteins x proteins, values in [0, 1]) from the
# STRING v12.0 human interaction scores (https://string-db.org, files
# 9606.protein.links.v12.0.txt and 9606.protein.info.v12.0.txt, downloaded into
# STRING_DIR). The combined score is rescaled to [0, 1] by min-max over all
# STRING links; pairs without a STRING link get 0. Proteins are matched by gene
# symbol (protein.info column 2) to the Olink assay names of the proteomics
# data.
# Run from the private project folder:
#   Rscript --vanilla ../T2D_application_code/preprocessing/STRING_preprocessing.R
# Outputs (in the private project folder, read by
# Network_estimation_case_study_T2D_UKB.R): STRING_W (Matrix Market),
# STRING_W_rownames.txt, STRING_W_colnames.txt, STRING_database.txt.
rm(list=ls(all=TRUE))
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "config.R"))

library(dplyr)
library(Matrix)

# Protein names: all Olink assays of the proteomics data (column names of
# prot_data in the input built by preprocessing/; the original run used an
# earlier build of the same proteomics file)
load(INPUT_RDATA)
prot_names <- colnames(prot_data)[-1]

string_db_path = file.path(STRING_DIR, "9606.protein.links.v12.0.txt")
string_info_path = file.path(STRING_DIR, "9606.protein.info.v12.0.txt")

string_db <- read.table(string_db_path, sep=' ', header=TRUE)

# Normalize the 'combined_score' column
normalize_column <- function(column) {
  (column - min(column)) / (max(column) - min(column))
}
string_db$combined_score <- normalize_column(string_db$combined_score)

# Read STRING database protein information data
string_db_info <- read.table(string_info_path, sep= "\t", quote = "", header=F)

# Check the number of common proteins between the dataset and STRING database
prot_common <- sum(string_db_info$V2 %in% prot_names)
cat("The number of common proteins between the dataset and STRING database: ", prot_common, "\n")

# Filter STRING database protein information for common proteins
string_db_info <- string_db_info[string_db_info$V2 %in% prot_names, ]

# Merge interaction data with protein information data for protein1
string_db <- string_db %>%
  inner_join(string_db_info, by=c("protein1" = "V1")) %>%
  select(protein1, protein2, combined_score, V2) %>%
  rename(Target1 = V2)

# Merge interaction data with protein information data for protein2
string_db <- string_db %>%
  inner_join(string_db_info, by=c("protein2" = "V1")) %>%
  select(protein1, protein2, combined_score, Target1, V2) %>%
  rename(Target2 = V2)

string_db <- string_db[, c("Target1","Target2","combined_score")]

write.table(string_db, file="STRING_database.txt")

# Create an empty adjacency matrix
adjacency_matrix <- Matrix(0, nrow = length(prot_names), ncol = length(prot_names), 
                           dimnames = list(prot_names, prot_names))

# Fill in the matrix with combined scores where applicable
for (i in 1:nrow(string_db)) {
  prot1 <- string_db[i, "Target1"]
  prot2 <- string_db[i, "Target2"]
  score <- string_db[i, "combined_score"]
  
  adjacency_matrix[prot1, prot2] <- score
  
}

isSymmetric(adjacency_matrix)

write(colnames(adjacency_matrix), file = "STRING_W_colnames.txt")
write(rownames(adjacency_matrix), file = "STRING_W_rownames.txt")
writeMM(adjacency_matrix, "STRING_W")

