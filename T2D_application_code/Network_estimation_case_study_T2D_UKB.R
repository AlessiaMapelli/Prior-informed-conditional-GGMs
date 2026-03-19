setwd("T2D_application_code")
rm(list=ls(all=TRUE))


source("ggReg_main_functions.R")

# PREPARE THE DATA
load("UKB_proteo_T2D.RData")
#Load a dataframe called "study_cohort" with 366 protein expression data and 3 covariates ("Sex", "Age", "Diab")

str(study_cohort)
colSums(is.na(study_cohort))
study_cohort <- na.omit(study_cohort)

# Cohort selection
library(caret)
library(dplyr)
library(tidyr)

set.seed(123)

# Step 1: Split the dataset into 80% and 20% maintaining the proportion of Diab and Sex
index_1 <- createDataPartition(study_cohort$Diab, p = 0.8, list = FALSE, times = 1)
selection_data <- study_cohort[index_1, ]
remaining_data <- study_cohort[-index_1, ]

table(selection_data$Diab, selection_data$Sex)
table(study_cohort$Diab, study_cohort$Sex)


# Step 2: From the remaining 20%, divide into two parts matching the training set on sex and age
incident_participants <- remaining_data %>% filter(Diab == "Incident")
non_diabetic_participants <- remaining_data %>% filter(Diab == "Non diabetic")

set.seed(123)
index_2 <- createDataPartition(incident_participants$Sex, p = 0.5, list = FALSE, times = 1)
incident_training_data <- incident_participants[index_2, ]
incident_testing_data <- incident_participants[-index_2, ]
train_data <- data.frame()
for (i in 1:nrow(incident_training_data)) {
  incident_row <- incident_training_data[i, ]
  set.seed(123)
  matching_non_diabetic <- non_diabetic_participants %>%
    filter(Age == incident_row$Age, Sex == incident_row$Sex) %>%
    sample_n(4, replace = FALSE)
  train_data <- rbind(train_data, incident_row, matching_non_diabetic)
  non_diabetic_participants <- anti_join(non_diabetic_participants, matching_non_diabetic)
}

test_data <- anti_join(remaining_data, train_data)


# Load the prior matrix procedded from STRING database
STRING_prior <- readMM("STRING_W")
colnames_W <- read.table("STRING_W_colnames.txt")$V1
rownames_W <- read.table("STRING_W_rownames.txt")$V1
dimnames(STRING_prior)<- list(rownames_W,colnames_W)


expression_data <- selection_data[,-(1:4)]
covariates <- selection_data[,c(2,4)]


# 
# weighted conditional GGMReg - mean regressed and covariance over the covriate
#
results <- GGReg_full_estimation(
  x=expression_data,
  known_ppi = as.matrix(STRING_prior),
  covariates = covariates,
  scr = FALSE,
  gamma = NULL,
  mean_estimation = TRUE,
  lambda_mean = NULL,
  lambda_mean_type = "min",
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = TRUE,
  p.rand.hyper = 0.8,
  K = 5,
  use_slurm = TRUE,
  slurm_script_path = "slurm_ggReg_node.sbatch",
  output_path = "./results/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = TRUE)

output_path = "./results/"
selection_data=selection_data[,-3]
save(results, selection_data, train_data, test_data, file = paste0(output_path, "PPI_T2D_app.RData"))

