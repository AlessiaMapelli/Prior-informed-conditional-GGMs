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

# Define diabetes
tsv_file <- file.path(UKB_RAW_DIR, "md_population_characteristics_19_12_2022_15_54_37.tsv")
population_df <- read.table(tsv_file, header = TRUE, sep = '\t')
colnames(population_df) <- c("id","sex","YOB","MOB","Etnic_group_visit_0",
                             "Etnic_group_visit_1","Etnic_group_visit_2","Age_at_recruitment",
                             "Gen_ethnic_group", "Date_of_death_1","Date_of_death_2")
population_df <- mutate(population_df, Date_of_death = pmin(as.Date(population_df$Date_of_death_1),as.Date(population_df$Date_of_death_2) , na.rm = TRUE))
population_df <- population_df[,-c(10,11)]

tsv_file <- paste0("md_population_characteristics_AM _", APP_ID, "_MD_20_02_2024_16_58_06.tsv")
assesment_date_df <- read.table(tsv_file, header = TRUE, sep = '\t')
assesment_date_df <- assesment_date_df[,c("f.eid" ,"Date_of_attending_assessment_centre.0.0")]
colnames(assesment_date_df)<- c("id","Assessment_date")


population_df <- merge(population_df, assesment_date_df, on="id")

tsv_file <- file.path(UKB_RAW_DIR, "md_diabetes_29_11_2022_10_38_01.tsv")
DIAB_df <- read.table(tsv_file, header = TRUE, sep = '\t')
colnames(DIAB_df) <- c("id","E10_date","E10_source","E11_date","E11_source",
                       "E12_date","E12_source","E13_date","E13_source","E14_date","E14_source")
DIAB_df <- merge(population_df, DIAB_df, on="id")
colnames(DIAB_df)

index_select <- DIAB_df$E10_date %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #1
DIAB_df <- DIAB_df[!index_select,]
index_select <- DIAB_df$E11_date %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #3
DIAB_df <- DIAB_df[!index_select,]
index_select <- DIAB_df$E12_date %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- DIAB_df$E13_date %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- DIAB_df$E14_date %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
sum(index_select)

var_considered <- c(12,14,16,18,20)
DIAB_dicotomy <-data.frame(matrix(nrow=dim(DIAB_df)[1], ncol=length(var_considered)))
colnames(DIAB_dicotomy) <- lapply(colnames(DIAB_df)[var_considered], function(x) substr(x, 1, 3))

for (i in 1:length(var_considered)){
  var <- DIAB_df[,var_considered[i]]
  DIAB_dicotomy[,i] <- !(var == "")
}

DIAB_dicotomy <- cbind(DIAB_dicotomy, DIAB_dicotomy[,2] | DIAB_dicotomy[,3] | DIAB_dicotomy[,4] | DIAB_dicotomy[,5])
colnames(DIAB_dicotomy)[6] <- "T2D_diagnosed"
sum(DIAB_dicotomy$T2D_diagnosed)


DIAB_dicotomy <- cbind(DIAB_dicotomy, lapply(DIAB_df[,c(11,var_considered)], as.Date))
colnames(DIAB_dicotomy)
DIAB_dicotomy  <- mutate(DIAB_dicotomy , T2D_diagnosed_date = pmin(DIAB_dicotomy[,9],DIAB_dicotomy[,10],DIAB_dicotomy[,11],DIAB_dicotomy[,12], na.rm = TRUE))
T2D_diagnosed_prevalent <- DIAB_dicotomy$Assessment_date > DIAB_dicotomy$T2D_diagnosed_date
DIAB_dicotomy <- cbind(DIAB_dicotomy,T2D_diagnosed_prevalent)
pheno_data <- cbind(DIAB_df[,-c(11,var_considered)],DIAB_dicotomy )

pheno_ids <- pheno_data$id

# Add other covariates
tsv_file <- paste0("md_population_characteristics_AM _", APP_ID, "_MD_20_02_2024_16_58_06.tsv")
assesment_date_df <- read.table(tsv_file, header = TRUE, sep = '\t')
assesment_date_df <- assesment_date_df[assesment_date_df$f.eid %in% pheno_ids,]

add_pheno_path = paste0("md_additional_info_AM_", APP_ID, "_MD_12_03_2024_10_23_36.tsv")
add_pheno_data <- read.csv(add_pheno_path, sep = "\t")
add_pheno_data <- add_pheno_data[add_pheno_data$f.eid %in% pheno_ids,]

add_pheno_path_2 = paste0("md_additional_info_AM_", APP_ID, "_MD_12_03_2024_18_08_22.tsv")
add_pheno_data_2 <- read.csv(add_pheno_path_2, sep = "\t")
add_pheno_data_2 <- add_pheno_data_2[add_pheno_data_2$f.eid %in% pheno_ids,]


add_pheno_path_3 <- paste0("md_additional_info_AM_", APP_ID, "_MD_29_03_2024_16_06_02.tsv")
add_pheno_data_3 <- read.table(add_pheno_path_3, header = TRUE, sep = '\t')
add_pheno_data_3 <- add_pheno_data_3[add_pheno_data_3$f.eid %in% pheno_ids,]



pheno_data$BMI <- add_pheno_data$Body_mass_index_.BMI..0.0
pheno_data$Alcohol_drinker_status <- add_pheno_data$Alcohol_drinker_status.0.0
pheno_data$Alcohol_intake_frequency <- add_pheno_data$Alcohol_intake_frequency..0.0
pheno_data$Smoking_status <- add_pheno_data$Smoking_status.0.0
pheno_data$Blood_type <- add_pheno_data$Blood.type_haplotype.0.0
pheno_data$Ethnic_background <- add_pheno_data$Ethnic_background.0.0
pheno_data[!is.na(pheno_data$Etnic_group_visit_1), "Ethnic_background"] <- pheno_data[!is.na(pheno_data$Etnic_group_visit_1), "Etnic_group_visit_1"] 
pheno_data[!is.na(pheno_data$Etnic_group_visit_2), "Ethnic_background"] <- pheno_data[!is.na(pheno_data$Etnic_group_visit_2), "Etnic_group_visit_2"] 
pheno_data$HbA1c <- add_pheno_data$Glycated_haemoglobin_.HbA1c..0.0
pheno_data$HbA1c_follow_up <- add_pheno_data$Glycated_haemoglobin_.HbA1c..1.0
pheno_data$Assessment_date_follow_up <- assesment_date_df$Date_of_attending_assessment_centre.1.0

str(pheno_data)
colnames(pheno_data)
pheno_data <- pheno_data[,-c(5,6,7)]
pheno_data$sex <- recode(as.factor(pheno_data$sex),
                         "0"  = "Female",
                         "1"  = "Male",
                         .default = NA_character_)
pheno_data$Gen_ethnic_group <-recode(as.factor(pheno_data$Gen_ethnic_group),
                                     "1"  = "Caucasian",
                                     .default = NA_character_)

pheno_data$Alcohol_drinker_status <- recode(as.factor(pheno_data$Alcohol_drinker_status),
                                            "-3" =  NA_character_,
                                            "0"  = "Never",
                                            "1"  = "Previous",
                                            "2"  = "Current",
                                            .default = NA_character_)

pheno_data$Alcohol_intake_frequency <- recode(as.factor(pheno_data$Alcohol_intake_frequency),
                                              "1" = "Daily or almost daily",
                                              "2" = "Three or four times a week",
                                              "3" = "Once or twice a week",
                                              "4" = "One to three times a month",
                                              "5" = "Special occasions only",
                                              "6" = "Never",
                                              "-3" = NA_character_,
                                              .default = NA_character_)

pheno_data$Smoking_status  <- recode(as.factor(pheno_data$Smoking_status ),
                                     "0" = "Never",
                                     "1" = "Previous",
                                     "2" = "Current",
                                     "-3" = NA_character_,
                                     .default = NA_character_)
pheno_data$Blood_type <- as.factor(pheno_data$Blood_type)
pheno_data$Assessment_date_follow_up <- as.Date(pheno_data$Assessment_date_follow_up)
pheno_data <- pheno_data %>%
  mutate(
    Ethnic_background = as.factor(Ethnic_background),
    Ethnic_background = case_when(
      Ethnic_background %in% c(-3, -1) ~ NA_character_,
      Ethnic_background %in% c(1, 1001, 1002, 1003) ~ "White",
      Ethnic_background %in% c(2, 2001, 2002, 2003, 2004) ~ "Mixed",
      Ethnic_background %in% c(3, 3001, 3002, 3003, 3004) ~ "Asian",
      Ethnic_background %in% c(4, 4001, 4002, 4003) ~ "Black",
      Ethnic_background %in% c(5) ~ "Chinese",
      Ethnic_background %in% c(6) ~ "Other",
      TRUE ~ as.factor(Ethnic_background)
    ),
    Ethnic_background = as.factor(Ethnic_background)
  )

str(pheno_data)

pheno_data$Pancreas_volume <- add_pheno_data$Pancreas_volume.2.0
pheno_data$Pancreas_PDFF <- add_pheno_data$Pancreas_PDFF_.fat_fraction..2.0
pheno_data$Abdominal_fat_ratio <- add_pheno_data$Abdominal_fat_ratio.2.0

index_select <- add_pheno_data$Date_E66_first_reported_.obesity..0.0  %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #9
pheno_data$Obesity_diag_prev <-ifelse(add_pheno_data$Date_E66_first_reported_.obesity..0.0 == "", NA,
                                      pheno_data$Assessment_date > as.Date(add_pheno_data$Date_E66_first_reported_.obesity..0.0))
pheno_data[index_select,"Obesity_diag_prev"] <- "Flag"
pheno_data$Obesity_date <- add_pheno_data$Date_E66_first_reported_.obesity..0.0
pheno_data[index_select,"Obesity_date"] <- "Flag"

add_pheno_path_4 = paste0("md_additional_info_AM_", APP_ID, "_MD_10_07_2024_16_36_39.tsv")
add_pheno_data_4 <- read.csv(add_pheno_path_4, sep = "\t")
add_pheno_data_4 <- add_pheno_data_4[add_pheno_data_4$f.eid %in% pheno_ids,]
add_pheno_data_4 <- merge(add_pheno_data_4, add_pheno_data[,c(1,11)])

add_pheno_path_5 <- paste0("md_additional_info_AM_", APP_ID, "_MD_11_07_2024_11_00_15.tsv")
add_pheno_data_5 <- read.csv(add_pheno_path_5, sep = "\t")
add_pheno_data_5 <- add_pheno_data_5[add_pheno_data_5$f.eid %in% pheno_ids,]
add_pheno_data_4 <- merge(add_pheno_data_4, add_pheno_data_5)


index_select <- add_pheno_data_4$Date_I11_first_reported_.hypertensive_heart_disease..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_4$Date_I12_first_reported_.hypertensive_heart_disease..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_4$Date_I13_first_reported_.hypertensive_heart_disease..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_4$Date_I10_first_reported_.hypertensive_heart_disease..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_4$Date_I15_first_reported_.secondary_hypertension..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
sum(index_select)

tot_hypertention  <-pmin(as.Date(add_pheno_data_4[,2]),as.Date(add_pheno_data_4[,3]),as.Date(add_pheno_data_4[,4]),as.Date(add_pheno_data_4[,5]), as.Date(add_pheno_data_4[,11]), na.rm = TRUE)

pheno_data$Hypertension_diag_prev <-ifelse(is.na(tot_hypertention), NA,
                                      pheno_data$Assessment_date > tot_hypertention)
pheno_data$Hypertension_date <- tot_hypertention

pheno_data$Physical_activity_days.week <- add_pheno_data$Number_of_days.week_of_moderate_physical_activity_10._minutes.0.0
pheno_data[pheno_data$Physical_activity_days.week %in% c(-3,-1),"Physical_activity_days.week"] <- NA
pheno_data$Physical_activity_vigorous_days.week <- add_pheno_data$Number_of_days.week_of_vigorous_physical_activity_10._minutes.0.0
pheno_data[pheno_data$Physical_activity_vigorous_days.week %in% c(-3,-1),"Physical_activity_vigorous_days.week"] <- NA

gen_pheno_data <- add_pheno_data[,c(45,46,47,48,49)]
colnames(gen_pheno_data)<- c("genPC1","genPC2","genPC3","genPC4","genPC5")
pheno_data <- cbind(pheno_data, gen_pheno_data)

pheno_data$Body_fat_percentage <- add_pheno_data$Body_fat_percentage.0.0
pheno_data$Whole_body_fat_mass <- add_pheno_data$Whole_body_fat_mass.0.0

pheno_data$Deprivation_index <- add_pheno_data_2$Townsend_deprivation_index_at_recruitment.0.0
pheno_data$Hip_circumference <- add_pheno_data_2$Hip_circumference.0.0


determine_obesity <- function(BMI, Ethnic_background) {
  case_when(
    Ethnic_background == "White" & BMI >= 30 ~ "Obese",
    Ethnic_background == "Mixed" & BMI >= 30 ~ "Obese",
    Ethnic_background == "Black" & BMI >= 28.1 ~ "Obese",
    Ethnic_background == "Chinese" & BMI >= 26.9 ~ "Obese",
    Ethnic_background == "Asian" & BMI >= 23.9 ~ "Obese",
    Ethnic_background == "Other" & BMI >= 30 ~ "Obese",
    TRUE ~ "Not Obese"
  )
}

pheno_data <- pheno_data %>%
  mutate(Obesity_BMI_based = determine_obesity(BMI, Ethnic_background),
         Obesity_BMI_based = as.factor(Obesity_BMI_based))

# Data_description.R
pheno_data$outcome_T2D_diagnosed <- pheno_data$T2D_diagnosed
pheno_data$outcome_T2D_diagnosed<- ifelse(pheno_data$outcome_T2D_diagnosed==FALSE, "Non diabetic", pheno_data$outcome_T2D_diagnosed)
table(pheno_data$T2D_diagnosed, pheno_data$T2D_diagnosed_prevalent)
prevalent_index <- ifelse(is.na(pheno_data$T2D_diagnosed_prevalent),FALSE, pheno_data$T2D_diagnosed_prevalent )
pheno_data[prevalent_index,"outcome_T2D_diagnosed"] <- "Prevalent"
incident_index <- ifelse(is.na(pheno_data$T2D_diagnosed_prevalent),TRUE, pheno_data$T2D_diagnosed_prevalent )
pheno_data[!incident_index,"outcome_T2D_diagnosed"] <- "Incident"
table(pheno_data$outcome_T2D_diagnosed,useNA ="ifany")


pheno_data$T2D_HbA1c <- pheno_data$HbA1c > 48
table(pheno_data$T2D_HbA1c,useNA ="ifany")
table(T2D_diagnosed= pheno_data$outcome_T2D_diagnosed,T2D_HbA1c=pheno_data$T2D_HbA1c,useNA ="ifany")

pheno_data$outcome_T2D_diagnosed.HbA1c <- as.character(pheno_data$outcome_T2D_diagnosed)

undiag_indexed <- (pheno_data$outcome_T2D_diagnosed  == "Incident" & pheno_data$T2D_HbA1c) | (pheno_data$outcome_T2D_diagnosed  == "Non diabetic" & pheno_data$T2D_HbA1c)
undiag_indexed[is.na(undiag_indexed)] <- FALSE

pheno_data[undiag_indexed,"outcome_T2D_diagnosed.HbA1c"] <- "Undiagnosed"
  
table(pheno_data$outcome_T2D_diagnosed.HbA1c,useNA ="ifany" )

# USE ANTIDIABETICS IN THE DEF
tsv_file <- paste0("gpp_drugs_antidiabetics_AM_", APP_ID, "_GP-P_11_07_2024_19_03_08.tsv")
antidiab_df <- read.table(tsv_file, header = TRUE, sep = '\t')

antidiab_df <- antidiab_df[, c("eid","issue_date")]
head(antidiab_df)

antidiab_df$issue_date <- as.Date(antidiab_df$issue_date, format = "%d/%m/%Y")

earliest_dates <- antidiab_df %>%
  group_by(eid) %>%
  summarise(earliest_date = min(issue_date, na.rm = TRUE))
head(earliest_dates)

latest_dates <- antidiab_df %>%
  group_by(eid) %>%
  summarise(latest_date = max(issue_date, na.rm = TRUE))
head(latest_dates)


antidiab_dates <- merge(earliest_dates,latest_dates, by="eid") 
head(antidiab_dates)
antidiab_dates$use_antidiab <- rep (TRUE,dim(antidiab_dates)[1] )

sum(antidiab_dates$eid %in% pheno_ids)
table(pheno_data[pheno_data$id %in% antidiab_dates$eid, "outcome_T2D_diagnosed.HbA1c" ])

colnames(antidiab_dates) <- c("id", "use_antidiab_earliest_date", "use_antidiab_latest_date", "use_antidiab")
pheno_data <- merge(pheno_data,  antidiab_dates, by="id", all.x=T)
pheno_data$use_antidiab_earliest_date <- as.Date(pheno_data$use_antidiab_earliest_date)
sum(is.infinite(pheno_data$use_antidiab_earliest_date))
sum(is.infinite(pheno_data$use_antidiab_latest_date))


pheno_data$outcome_T2D_antidiab <- ifelse(is.na(pheno_data$use_antidiab_earliest_date), "No Antidiab",
                                  ifelse(pheno_data$Assessment_date > pheno_data$use_antidiab_earliest_date,
                                         ifelse(pheno_data$Assessment_date - pheno_data$use_antidiab_latest_date > 365, "No Antidiab","Antidiab_at_Ass_date"),
                                         "Antidiab_after_Ass_date"))

table(pheno_data$outcome_T2D_antidiab, useNA = "ifany")
table(pheno_data$outcome_T2D_antidiab, pheno_data$outcome_T2D_diagnosed.HbA1c, useNA = "ifany")

pheno_data$outcome_T2D_diagnosed.HbA1c.drug <- pheno_data$outcome_T2D_diagnosed.HbA1c

add_prevalent_indexed <- pheno_data$outcome_T2D_diagnosed.HbA1c  == "Non diabetic" & pheno_data$outcome_T2D_antidiab == "Antidiab_at_Ass_date"
table(pheno_data[add_prevalent_indexed,"E10"], useNA = "ifany")
pheno_data[add_prevalent_indexed,"outcome_T2D_diagnosed.HbA1c.drug"] <- "Flag_ND"

add_prevalent_indexed <- pheno_data$outcome_T2D_diagnosed.HbA1c  == "Incident" & pheno_data$outcome_T2D_antidiab == "Antidiab_at_Ass_date"
table(pheno_data[add_prevalent_indexed,"E10"], useNA = "ifany")
pheno_data[add_prevalent_indexed,"outcome_T2D_diagnosed.HbA1c.drug"] <- "Flag_I"

add_incident_indexed <- pheno_data$outcome_T2D_diagnosed.HbA1c  == "Non diabetic" & pheno_data$outcome_T2D_antidiab == "Antidiab_after_Ass_date"
table(pheno_data[add_incident_indexed,"E10"], useNA = "ifany")
pheno_data[add_incident_indexed,"outcome_T2D_diagnosed.HbA1c.drug"] <- "Flag_ND"

table(pheno_data$outcome_T2D_diagnosed.HbA1c.drug, useNA = "ifany")

rm(antidiab_df)
# write.csv(pheno_data, paste0("info_DIAB_", APP_ID, "_add_cov_20240326.csv"), row.names=TRUE)

# Basic_prediction_models.R
Age_at_recruitment.2 <- pheno_data$Age_at_recruitment*pheno_data$Age_at_recruitment
pheno_data <- cbind(pheno_data, Age_at_recruitment.2)

additional_covar = read.csv(paste0("md_additional_info_AM_", APP_ID, "_MD_15_03_2024_11_47_34.tsv"), sep ="\t")
additional_covar = additional_covar[,c(1,2)]
colnames(additional_covar) = c("id", "Number_of_treatments")
pheno_data <- left_join(pheno_data, additional_covar,by="id")

colnames(pheno_data)

# LOAD ADDITIONAL COVARIATES
colnames(add_pheno_data_3)[1] <- "id"
colnames(add_pheno_data_3)
colSums(is.na(add_pheno_data_3))
add_pheno_data_3$Assessment_date <- pheno_data$Assessment_date


clinical_model_df <- pheno_data[,c(1,2,5,27,32,30,52,33)]
colnames(clinical_model_df) <- c("id","sex","age","bmi","ethn","smoke", "tds", "hba1c")
colnames(add_pheno_data_3)[2] <- "PRS_T2D.UKB"
clinical_model_df$PRS_T2D.UKB <- add_pheno_data_3$PRS_T2D.UKB
pheno_data$PRS_T2D.UKB <- add_pheno_data_3$PRS_T2D.UKB

index_select <- add_pheno_data_3$Date_E28_first_reported_.ovarian_dysfunction..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
sum(index_select)

pheno_data$Ovarian_dysfunction_date <- add_pheno_data_3$Date_E28_first_reported_.ovarian_dysfunction..0.0
pheno_data$Ovarian_dysfunction_diag_prev <-ifelse(add_pheno_data_3$Date_E28_first_reported_.ovarian_dysfunction..0.0 == "", FALSE,
                                                  as.Date(add_pheno_data_3$Assessment_date) > as.Date(add_pheno_data_3$Date_E28_first_reported_.ovarian_dysfunction..0.0))

clinical_model_df$pcos <- ifelse(is.na(pheno_data$Ovarian_dysfunction_diag_prev), FALSE, pheno_data$Ovarian_dysfunction_diag_prev)

colnames(add_pheno_data_5)[1] <- "id"
add_pheno_data_3 <- merge(add_pheno_data_3,add_pheno_data_5)

index_select <- add_pheno_data_3$Date_F25_first_reported_.schizoaffective_disorders..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_F24_first_reported_.induced_delusional_disorder..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_F20_first_reported_.schizophrenia..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_F21_first_reported_.schizotypal_disorder..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_F22_first_reported_.persistent_delusional_disorders..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
sum(index_select)

tot_schizophrenia <-pmin(as.Date(add_pheno_data_3[,144]),as.Date(add_pheno_data_3[,145]),as.Date(add_pheno_data_3[,146]),as.Date(add_pheno_data_3[,147]), as.Date(add_pheno_data_3[,4]), na.rm = TRUE)

pheno_data$Schizophrenia_diag_prev <-ifelse(is.na(tot_schizophrenia), NA,
                                           pheno_data$Assessment_date > tot_schizophrenia)
pheno_data$Schizophrenia_date <- tot_schizophrenia

index_select <- add_pheno_data_3$Date_F30_first_reported_.manic_episode..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_F31_first_reported_.bipolar_affective_disorder..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #1
sum(index_select)

tot_bipolar <-pmin(as.Date(add_pheno_data_3[,148]),as.Date(add_pheno_data_3[,5]), na.rm = TRUE)

pheno_data$Bipolar_diag_prev <-ifelse(is.na(tot_bipolar), NA,
                                            pheno_data$Assessment_date > tot_bipolar)
pheno_data$Bipolar_date <- tot_bipolar

clinical_model_df$psy <- ifelse(is.na((pheno_data$Schizophrenia_diag_prev | pheno_data$Bipolar_diag_prev)), FALSE, (pheno_data$Schizophrenia_diag_prev | pheno_data$Bipolar_diag_prev))

pheno_data[index_select,"Bipolar_diag_prev"] <- "Flag"
pheno_data$Bipolar_date <- as.character(pheno_data$Bipolar_date)
pheno_data[index_select,"Bipolar_date"] <- "Flag"
clinical_model_df[index_select,"psy"] <- "Flag"

index_select <- add_pheno_data_3$Date_O24_first_reported_.diabetes_mellitus_in_pregnancy..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #2
sum(index_select)

pheno_data$Gestational_diab_diag_prev <-ifelse(add_pheno_data_3$Date_O24_first_reported_.diabetes_mellitus_in_pregnancy..0.0 == "", FALSE,
                                                   add_pheno_data_3$Assessment_date > as.Date(add_pheno_data_3$Date_O24_first_reported_.diabetes_mellitus_in_pregnancy..0.0))
pheno_data$Gestational_diab_date <-add_pheno_data_3$Date_O24_first_reported_.diabetes_mellitus_in_pregnancy..0.0
pheno_data[index_select,"Gestational_diab_diag_prev"] <- "Flag"
pheno_data[index_select,"Gestational_diab_date"] <- "Flag"

clinical_model_df$gdm <- ifelse(is.na(pheno_data$Gestational_diab_diag_prev), FALSE, pheno_data$Gestational_diab_diag_prev)

clinical_model_df$htn <- ifelse(is.na(pheno_data$Hypertension_diag_prev), FALSE, pheno_data$Hypertension_diag_prev)

last_col_path <- paste0("md_additional_info_AM_", APP_ID, "_MD_02_04_2024_19_05_22.tsv")
last_col <- read.csv(last_col_path,sep = '\t')
last_col <- last_col[last_col$f.eid %in% pheno_ids,]
colnames(last_col)[1] <- "id"
add_pheno_data_3 <- merge(add_pheno_data_3,last_col)

index_select <- add_pheno_data_3$Date_I23_first_reported_.certain_current_complications_following_acute_myocardial_infarction..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_I21_first_reported_.acute_myocardial_infarction..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_I60_first_reported_.subarachnoid_haemorrhage..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_I22_first_reported_.subsequent_myocardial_infarction..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_I62_first_reported_.other_nontraumatic_intracranial_haemorrhage..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_I61_first_reported_.intracerebral_haemorrhage..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_I63_first_reported_.cerebral_infarction..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #0
index_select <- add_pheno_data_3$Date_I64_first_reported_.stroke._not_specified_as_haemorrhage_or_infarction..0.0 %in% c("1900-01-01","1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07") #2
sum(index_select)

tot_cardiovascular_disease <-pmin(as.Date(add_pheno_data_3[,"Date_I23_first_reported_.certain_current_complications_following_acute_myocardial_infarction..0.0"]),
                         as.Date(add_pheno_data_3[,"Date_I21_first_reported_.acute_myocardial_infarction..0.0"]),
                         as.Date(add_pheno_data_3[,"Date_I60_first_reported_.subarachnoid_haemorrhage..0.0"]),
                         as.Date(add_pheno_data_3[,"Date_I22_first_reported_.subsequent_myocardial_infarction..0.0"]),
                         as.Date(add_pheno_data_3[,"Date_I62_first_reported_.other_nontraumatic_intracranial_haemorrhage..0.0"]),
                         as.Date(add_pheno_data_3[,"Date_I61_first_reported_.intracerebral_haemorrhage..0.0"]),
                         as.Date(add_pheno_data_3[,"Date_I63_first_reported_.cerebral_infarction..0.0"]),
                         as.Date(add_pheno_data_3[,"Date_I64_first_reported_.stroke._not_specified_as_haemorrhage_or_infarction..0.0"]),
                         na.rm = TRUE)

pheno_data$Cardiovascular_disease_diag_prev <-ifelse(is.na(tot_cardiovascular_disease), NA,
                                            pheno_data$Assessment_date > tot_cardiovascular_disease)
pheno_data$Cardiovascular_disease_date <- tot_cardiovascular_disease

pheno_data[index_select,"Cardiovascular_disease_diag_prev"] <- "Flag"
pheno_data$Cardiovascular_disease_date <- as.character(pheno_data$Cardiovascular_disease_date)
pheno_data[index_select,"Cardiovascular_disease_date"] <- "Flag"

clinical_model_df$cvd <- ifelse(is.na(pheno_data$Cardiovascular_disease_diag_prev), FALSE, pheno_data$Cardiovascular_disease_diag_prev)
clinical_model_df$outcome <- pheno_data$outcome_T2D_diagnosed.HbA1c.drug

colnames(clinical_model_df)

# Father diab
add_pheno_data_3 <- add_pheno_data_3 %>%
  mutate(across(11:20, ~ifelse(is.na(.), -50, .)))

father_diab_dicotomy <- ( add_pheno_data_3[,11]== 9 | 
                            add_pheno_data_3[,12]== 9 |
                            add_pheno_data_3[,13]== 9 |
                            add_pheno_data_3[,14]== 9 |
                            add_pheno_data_3[,15]== 9 |
                            add_pheno_data_3[,16]== 9 |
                            add_pheno_data_3[,17]== 9 |
                            add_pheno_data_3[,18]== 9 |
                            add_pheno_data_3[,19]== 9 |
                            add_pheno_data_3[,20]== 9)
table(father_diab_dicotomy, useNA = "ifany")

pheno_data$Hist_father_diab <- father_diab_dicotomy

# Mother diab
add_pheno_data_3 <- add_pheno_data_3 %>%
  mutate(across(51:61, ~ifelse(is.na(.), -50, .)))

m_diab_dicotomy <- ( add_pheno_data_3[,51]== 9 | 
                            add_pheno_data_3[,52]== 9 |
                            add_pheno_data_3[,53]== 9 |
                            add_pheno_data_3[,54]== 9 |
                            add_pheno_data_3[,55]== 9 |
                            add_pheno_data_3[,56]== 9 |
                            add_pheno_data_3[,57]== 9 |
                            add_pheno_data_3[,58]== 9 |
                            add_pheno_data_3[,59]== 9 |
                            add_pheno_data_3[,60]== 9 |
                       add_pheno_data_3[,61]== 9)

table(m_diab_dicotomy, useNA = "ifany")

pheno_data$Hist_mother_diab <- m_diab_dicotomy

# sibilings diab
add_pheno_data_3 <- add_pheno_data_3 %>%
  mutate(across(95:106, ~ifelse(is.na(.), -50, .)))

s_diab_dicotomy <- ( add_pheno_data_3[,95]== 9 | 
                       add_pheno_data_3[,96]== 9 |
                       add_pheno_data_3[,97]== 9 |
                       add_pheno_data_3[,98]== 9 |
                       add_pheno_data_3[,99]== 9 |
                       add_pheno_data_3[,100]== 9 |
                       add_pheno_data_3[,101]== 9 |
                       add_pheno_data_3[,102]== 9 |
                       add_pheno_data_3[,103]== 9 |
                       add_pheno_data_3[,104]== 9 |
                       add_pheno_data_3[,105]== 9 |
                       add_pheno_data_3[,106]== 9)

table(s_diab_dicotomy, useNA = "ifany")

pheno_data$Hist_siblings_diab <- s_diab_dicotomy

pheno_data$Hist_family_diab <- pheno_data$Hist_siblings_diab | pheno_data$Hist_mother_diab | pheno_data$Hist_father_diab

clinical_model_df$fhdm <- pheno_data$Hist_family_diab


corticosteroids_path <- paste0("gpp_drugs_corticosteroids_AM_", APP_ID, "_GP-P_03_04_2024_12_43_35.tsv")
corticosteroids <- read.csv(corticosteroids_path,sep = '\t')
corticosteroids <- corticosteroids[, c("eid","issue_date")]
corticosteroids$issue_date <- as.Date(corticosteroids$issue_date, format = "%d/%m/%Y")

corticosteroids <- na.omit(corticosteroids)

earliest_dates <- corticosteroids %>%
  group_by(eid) %>%
  summarise(earliest_date = min(issue_date, na.rm = TRUE))

latest_dates <- corticosteroids %>%
  group_by(eid) %>%
  summarise(latest_date = max(issue_date, na.rm = TRUE))

sum(latest_dates$eid %in% pheno_data$id)

corticosteroids_dates <- merge(earliest_dates,latest_dates, by="eid") 
head(corticosteroids_dates)
str(corticosteroids_dates)

colnames(corticosteroids_dates) <- c("id", "use_corticosteroids_earliest_date", "use_corticosteroids_latest_date")
pheno_data <- merge(pheno_data,  corticosteroids_dates, by="id", all.x=T)

pheno_data$Hist_corticosteroids <- ifelse(is.na(pheno_data$use_corticosteroids_earliest_date), FALSE,
                                          ifelse(pheno_data$Assessment_date > pheno_data$use_corticosteroids_earliest_date,
                                                ifelse(pheno_data$Assessment_date - pheno_data$use_corticosteroids_latest_date > 365, FALSE,TRUE),
                                                FALSE))

table(pheno_data$Hist_corticosteroids, useNA = "ifany")

clinical_model_df$ster <- pheno_data$Hist_corticosteroids


statin_path <- file.path(UKB_RAW_DIR, "gp_prescription_19_12_2022_15_16_47.tsv")
statin <- read.csv(statin_path,sep = '\t')
statin <- statin[, c("eid","issue_date")]
head(statin)
statin$issue_date <- as.Date(statin$issue_date, format = "%d/%m/%Y")
statin <- na.omit(statin)

earliest_dates <- statin %>%
  group_by(eid) %>%
  summarise(earliest_date = min(issue_date, na.rm = TRUE))

latest_dates <- statin %>%
  group_by(eid) %>%
  summarise(latest_date = max(issue_date, na.rm = TRUE))

sum(latest_dates$eid %in% pheno_data$id)

statin_dates <- merge(earliest_dates,latest_dates, by="eid") 
head(statin_dates)

colnames(statin_dates) <- c("id", "use_statin_earliest_date", "use_statin_latest_date")
pheno_data <- merge(pheno_data,  statin_dates, by="id", all.x=T)

pheno_data$Hist_statin <- ifelse(is.na(pheno_data$use_statin_earliest_date), FALSE,
                                 ifelse(pheno_data$Assessment_date > pheno_data$use_statin_earliest_date,
                                        ifelse(pheno_data$Assessment_date - pheno_data$use_statin_latest_date > 365, FALSE,TRUE),
                                        FALSE))

table(pheno_data$Hist_statin, useNA = "ifany")
clinical_model_df$stat <- pheno_data$Hist_statin

antipsycotic_path <- paste0("gpp_drugs_antipsycotic_AM_", APP_ID, "_GP-P_03_04_2024_15_08_19.tsv")
antipsycotic <- read.csv(antipsycotic_path,sep = '\t')
antipsycotic <- antipsycotic[, c("eid","issue_date")]
head(antipsycotic)
antipsycotic$issue_date <- as.Date(antipsycotic$issue_date, format = "%d/%m/%Y")

antipsycotic <- na.omit(antipsycotic)

earliest_dates <- antipsycotic %>%
  group_by(eid) %>%
  summarise(earliest_date = min(issue_date, na.rm = TRUE))

latest_dates <- antipsycotic %>%
  group_by(eid) %>%
  summarise(latest_date = max(issue_date, na.rm = TRUE))

sum(latest_dates$eid %in% pheno_data$id)

antipsycotic_dates <- merge(earliest_dates,latest_dates, by="eid") 
head(antipsycotic_dates)

colnames(antipsycotic_dates) <- c("id", "use_antipsycotic_earliest_date", "use_antipsycotic_latest_date")
pheno_data <- merge(pheno_data,  antipsycotic_dates, by="id", all.x=T)

pheno_data$Hist_antipsycotic <- ifelse(is.na(pheno_data$use_antipsycotic_earliest_date), FALSE,
                                       ifelse(pheno_data$Assessment_date > pheno_data$use_antipsycotic_earliest_date,
                                              ifelse(pheno_data$Assessment_date - pheno_data$use_antipsycotic_latest_date > 365, FALSE,TRUE),
                                              FALSE))
table(pheno_data$Hist_antipsycotic, useNA = "ifany")

clinical_model_df$apsy <- pheno_data$Hist_antipsycotic

pheno_data$Obesity_diag_prev <- ifelse(is.na(pheno_data$Obesity_diag_prev), "Non Obese", pheno_data$Obesity_diag_prev)
pheno_data$Obesity_diag_prev <- recode(as.factor(pheno_data$Obesity_diag_prev),
                                 "TRUE"  = "Prevalent",
                                 "FALSE"  = "Incident")
table(pheno_data$Obesity_diag_prev, useNA = "ifany")

pheno_data$Hypertension_diag_prev <- ifelse(is.na(pheno_data$Hypertension_diag_prev), "Non Hypertensive", pheno_data$Hypertension_diag_prev)
pheno_data$Hypertension_diag_prev <- recode(as.factor(pheno_data$Hypertension_diag_prev),
                                      "TRUE"  = "Prevalent",
                                      "FALSE"  = "Incident")
table(pheno_data$Hypertension_diag_prev, useNA = "ifany")




follow_up_path <- paste0("md_additional_info_AM_", APP_ID, "_MD_15_07_2024_15_06_17.tsv")
follow_up_frame <- read.csv(follow_up_path,sep = '\t')
follow_up_frame <- follow_up_frame[,c(1:2)]
colnames(follow_up_frame) <- c("id", "Date_lost_to_follow.up")

pheno_data <- merge(pheno_data,  follow_up_frame, by="id", all.x=T)


# (in the original run this file was in the parent of UKB_WORK_DIR: set UKB_CPR_FILE)
CPR_path <- Sys.getenv("UKB_CPR_FILE", paste0("md_additional_info_AM_", APP_ID, "_MD_18_09_2025_09_40_46.tsv"))
CPR_df <-read.csv(CPR_path,sep = '\t')
CPR_df <- CPR_df[, c("f.eid","C.reactive_protein.0.0",
                     "Medication_for_cholesterol._blood_pressure._diabetes._or_take_exogenous_hormones.0.0",
                     "Medication_for_cholesterol._blood_pressure._diabetes._or_take_exogenous_hormones.0.1",
                     "Medication_for_cholesterol._blood_pressure._diabetes._or_take_exogenous_hormones.0.2",
                     "Medication_for_cholesterol._blood_pressure._diabetes._or_take_exogenous_hormones.0.3",
                     "Medication_for_cholesterol._blood_pressure_or_diabetes.0.0",
                     "Medication_for_cholesterol._blood_pressure_or_diabetes.0.1",
                     "Medication_for_cholesterol._blood_pressure_or_diabetes.0.2",
                     "Fasting_time.0.0",
                     "Time_blood_sample_collected.0.0",
                     "Time_blood_sample_collected.0.1",
                     "Time_blood_sample_collected.0.2",
                     "Time_blood_sample_collected.0.3",
                     "Time_blood_sample_collected.0.4",
                     "Time_blood_sample_collected.0.5",
                     "Time_blood_sample_collected.0.6"
                     )]
# shorthand for columns
med1 <- CPR_df$Medication_for_cholesterol._blood_pressure._diabetes._or_take_exogenous_hormones.0.0
med2 <- CPR_df$Medication_for_cholesterol._blood_pressure_or_diabetes.0.0

# define missing codes
missing_codes <- c(-7, -1, -3)

CPR_df$Cholesterol_lowering_medication <- case_when(
  med1 == 1 | med2 == 1 ~ TRUE,
  (is.na(med1) & is.na(med2)) ~ NA,
  (is.na(med1) & med2 %in% missing_codes) ~ NA,
  (is.na(med2) & med1 %in% missing_codes) ~ NA,
  TRUE ~ FALSE
)

# check distribution
table(CPR_df$Cholesterol_lowering_medication, useNA = "ifany")

time_cols <- grep("^Time_blood_sample_collected\\.", names(CPR_df), value = TRUE)

CPR_df <- CPR_df %>%
  mutate(across(all_of(time_cols), ~{
    x <- trimws(.)                                # strip spaces
    x[x == ""] <- NA_character_                   # treat empty as NA
    dt <- as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
    # Set out-of-range years to NA
    yr <- as.integer(format(dt, "%Y"))
    dt[!is.na(dt) & (yr < 2006 | yr > 2010)] <- NA
    dt
  })) %>%
  rowwise() %>%
  mutate(
    blood_time_min = {
      v <- c_across(all_of(time_cols))
      m <- suppressWarnings(min(v, na.rm = TRUE))
      if (is.infinite(m)) NA_real_ else as.numeric(m)
    },
    blood_time_min = as.POSIXct(blood_time_min, origin = "1970-01-01", tz = "UTC")
  ) %>%
  ungroup()

CPR_df <- CPR_df %>%
  mutate(
    Month_blood_sample_collected = factor(format(blood_time_min, "%B"),
                         levels = month.name, ordered = TRUE)
  )

CPR_df <- CPR_df %>%
  mutate(
    minutes_since_midnight = ifelse(
      is.na(blood_time_min), NA_integer_,
      as.integer(format(blood_time_min, "%H")) * 60 +
        as.integer(format(blood_time_min, "%M"))
    ),
    Time_of_day_blood_sample_collected = case_when(
      is.na(minutes_since_midnight)              ~ NA_character_,
      minutes_since_midnight <  8*60             ~ "early_morning",
      minutes_since_midnight < 12*60             ~ "morning",
      minutes_since_midnight < 16*60             ~ "afternoon",
      minutes_since_midnight < 19*60             ~ "late_afternoon",
      TRUE                                       ~ "night"
    ),
    Time_of_day_blood_sample_collected = factor(Time_of_day_blood_sample_collected,
                         levels = c("early_morning","morning","afternoon",
                                    "late_afternoon","night"),
                         ordered = TRUE)
  ) %>%
  dplyr::select(-minutes_since_midnight)
  


CPR_df <- CPR_df[, c("f.eid","C.reactive_protein.0.0","Cholesterol_lowering_medication",
                     "Fasting_time.0.0", "Month_blood_sample_collected", "Time_of_day_blood_sample_collected")]

colnames(CPR_df) <- c("id","C_reactive_protein","Cholesterol_lowering_medication", "Fasting_time",
                      "Month_blood_sample_collected", "Time_of_day_blood_sample_collected")
str(CPR_df)

pheno_data <- merge(pheno_data,  CPR_df, by="id", all.x=T)

write.csv(pheno_data, "full_pheno_data_20250911.csv", row.names=TRUE)
write.csv(clinical_model_df, "QDiab2018_data_20250911.csv", row.names=TRUE)






