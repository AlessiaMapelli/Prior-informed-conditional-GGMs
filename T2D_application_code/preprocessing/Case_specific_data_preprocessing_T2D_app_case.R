# Paths (set as environment variables; see ../README.md):
#   UKB_WORK_DIR     working folder: holds the extracted UK Biobank tables
#                    (outputs of extract_fields/*.yml) and receives the outputs
#   UKB_RAW_DIR      folder with the centrally extracted UK Biobank tables
#                    (population characteristics, diabetes, GP prescriptions)
#   UKB_PPP_DIR      folder with the UKB-PPP Olink proteomics files
#   UKB_APP_ID       UK Biobank application ID that appears in the names of the
#                    extracted files
UKB_WORK_DIR <- Sys.getenv("UKB_WORK_DIR", ".")
UKB_RAW_DIR  <- Sys.getenv("UKB_RAW_DIR", "raw_data")
UKB_PPP_DIR  <- Sys.getenv("UKB_PPP_DIR", "ukb-ppp")
APP_ID       <- Sys.getenv("UKB_APP_ID", "APPID")
setwd(UKB_WORK_DIR)
rm(list=ls())
graphics.off()


library(dplyr)
library(tidyr)

# Import full phenotype data
pheno_data <- read.csv("full_pheno_data_20250911.csv")[,-1]
clinical_model_df <- read.csv("QDiab2018_data_20250911.csv")[,-1]

pheno_data_diab <- pheno_data[pheno_data$outcome_T2D_diagnosed.HbA1c.drug %in% c("Non diabetic", "Incident"),]
pheno_data_diab <- pheno_data_diab[!pheno_data_diab$E10,]

pheno_data_diab$outcome_T2D_diagnosed.HbA1c.drug <- as.factor(pheno_data_diab$outcome_T2D_diagnosed.HbA1c.drug) 
table(pheno_data_diab$outcome_T2D_diagnosed.HbA1c.drug, useNA = "ifany")
pheno_data_diab$Assessment_date <- as.Date(pheno_data_diab$Assessment_date)
pheno_data_diab$T2D_diagnosed_date <- as.Date(pheno_data_diab$T2D_diagnosed_date)
pheno_data_diab$Date_lost_to_follow.up <- as.Date(pheno_data_diab$Date_lost_to_follow.up)
pheno_data_diab$Date_of_death <- as.Date(pheno_data_diab$Date_of_death)

pheno_data_diab$ten_y_end_of_study = as.Date(pheno_data_diab$Assessment_date) + 3650

pheno_data_diab$end_of_study <-  if_else(!is.na(pheno_data_diab$Date_lost_to_follow.up) & (pheno_data_diab$ten_y_end_of_study > pheno_data_diab$Date_lost_to_follow.up),
  as.Date(pheno_data_diab$Date_lost_to_follow.up), as.Date(pheno_data_diab$ten_y_end_of_study))
pheno_data_diab$LTFU = !is.na(pheno_data_diab$Date_lost_to_follow.up) & (pheno_data_diab$ten_y_end_of_study > pheno_data_diab$Date_lost_to_follow.up)

pheno_data_diab$end_of_study = if_else(!is.na(pheno_data_diab$Date_of_death ) & (pheno_data_diab$end_of_study > pheno_data_diab$Date_of_death ), as.Date(pheno_data_diab$Date_of_death) , as.Date(pheno_data_diab$end_of_study ))

head(pheno_data_diab[, c("Date_of_death","ten_y_end_of_study","Date_lost_to_follow.up", "end_of_study","LTFU")], 50)

pheno_data_diab$surv_time = if_else( (is.na(pheno_data_diab$T2D_diagnosed_date) | pheno_data_diab$T2D_diagnosed_date > pheno_data_diab$end_of_study), as.numeric(pheno_data_diab$end_of_study - pheno_data_diab$Assessment_date), as.numeric(pheno_data_diab$T2D_diagnosed_date - pheno_data_diab$Assessment_date))
pheno_data_diab$surv_status = if_else( (!is.na(pheno_data_diab$T2D_diagnosed_date) & pheno_data_diab$T2D_diagnosed_date <= pheno_data_diab$end_of_study), 1,0)
summary(pheno_data_diab$surv_time)

pheno_data_diab$early_onset <- pheno_data_diab$surv_time<30
write.csv(pheno_data_diab, paste0("complete_DIAB_analysis_", APP_ID, "_20260226_early onset_T2D_app.csv"))

pheno_data_diab = subset(pheno_data_diab, select = -c(early_onset) )

pheno_data_diab <- pheno_data_diab[pheno_data_diab$surv_time>30, ]

prot_path <- file.path(UKB_PPP_DIR, paste0("proteomics_wide_filtered_imputed_assay_", APP_ID, ".csv"))
prot_data <- read.table(prot_path, sep = ",", header=T)
colnames(prot_data)[1]= "id"

data <- merge(pheno_data_diab, prot_data, by = "id")
str(data)

five_y_end_of_study = as.Date(pheno_data_diab$Assessment_date) + 1825
pheno_data_diab$LTFU_five_y = !is.na(pheno_data_diab$Date_lost_to_follow.up) & (five_y_end_of_study > pheno_data_diab$Date_lost_to_follow.up)
pheno_data_diab$surv_status_5y = if_else( (!is.na(pheno_data_diab$T2D_diagnosed_date) & pheno_data_diab$T2D_diagnosed_date <= five_y_end_of_study), 1,0)
pheno_data_diab$surv_status_5y = if_else(pheno_data_diab$LTFU_five_y,NA, pheno_data_diab$surv_status_5y)

colnames(clinical_model_df)
str(clinical_model_df)
colnames(clinical_model_df)[9] <- "PRS_T2D"

clinical_model_df[clinical_model_df$psy == "Flag", "psy"] <- NA
clinical_model_df[clinical_model_df$gdm == "Flag", "gdm"] <- NA
clinical_model_df[clinical_model_df$cvd == "Flag", "cvd"] <- NA

colSums(is.na(clinical_model_df))

clinical_model_df$psy  <- as.logical(clinical_model_df$psy)
clinical_model_df$gdm  <- as.logical(clinical_model_df$gdm)
clinical_model_df$cvd  <- as.logical(clinical_model_df$cvd)
clinical_model_df$ethn <- as.factor(clinical_model_df$ethn)

clinical_model_df<- clinical_model_df %>%
  filter(!is.na(ethn)) %>%
  dplyr::select(-outcome)

clinical_var_biomarkers_PRS <- colnames(clinical_model_df)[-1]
clinical_var_biomarkers <- clinical_var_biomarkers_PRS[clinical_var_biomarkers_PRS != "PRS_T2D"]
clinical_var <- clinical_var_biomarkers[clinical_var_biomarkers != "hba1c"]
prot_conv_path <- file.path(UKB_PPP_DIR, "uniprot_panel.rtf")
prot_coversion <- read.csv(prot_conv_path, sep = "\t", header=T)
prot_names <- colnames(prot_data)[-1]
prot_coversion <- prot_coversion[prot_coversion$Assay %in% prot_names,]
table(prot_coversion$Panel)

data <- merge(pheno_data_diab, prot_data, by = "id")
str(data)

Cardiometabolic_prot <- prot_coversion[prot_coversion$Panel == "Cardiometabolic",]
full_omic_var <- Cardiometabolic_prot$Assay
model_data <- merge(data, clinical_model_df, by =c("id", "sex"))
logistic_outcome <- "surv_status_5y"
surv_status="surv_status"
surv_time = "surv_time"

save(pheno_data,pheno_data_diab ,prot_data, data,clinical_model_df,model_data, clinical_var,clinical_var_biomarkers, clinical_var_biomarkers_PRS, full_omic_var,logistic_outcome, surv_status, surv_time , file = Sys.getenv("T2D_INPUT_RDATA", "T2D_app_input.RData"))







