### script to merge the clinical data from the bbs main and sub studies, so it is ready to merge with the metabolon metabolite data ###

####################
###### SET UP ######
####################

## read in parameter file
# parameter file
source("parameter_files/03.parameters.R")

## set working directory
#setwd(working_dir)

## capture output
sink(file=paste0(log_output_dir, "03.merge.bbs.clinical.data.txt"),split=TRUE)

## load libraries
library(tidyverse)
library(data.table)

###############################################################################
############################### READ IN DATA ##################################
###############################################################################

# read in the processed data from main and sub studies
bbs_main <- read.csv(paste0(data_output_dir, "01_clinical_data_bbsmain_wide.csv"))
bbs_substudy <- read.csv(paste0(data_output_dir, "02_clinical_data_bbssubstudy_wide.csv"))

# read in manifest (linker file)
manifest <- read.csv(paste0(inventory_dir, "manifest_linker_file.csv"))
dim(manifest)
# remove alspac samples from manifest
manifest <- manifest[-(which(manifest$site=="alspac")),]
dim(manifest)

###############################################################################
################################ CHOOSE VARS ##################################
###############################################################################

# standardise the common column names
setnames(bbs_substudy,
         old = c("studyId","baseline_wgt","op_wgt","diab_onset_m","diab_onset_y",
                 "diab_hba1c_lev","diab_typeone","diab_typetwo",
                 "blood_date_base","blood_date_m12","wgt_m12",
                 "diab_hba1c_12m_lev","bmikgm2_36months",
                 "timeinfreezer_36months","age_at_consent_years","ppt_gender"),
         new = c("studyid","wgt_Baseline","wgt_op_day","diab_date_m",
                 "diab_date_y","hba1c_Baseline","type1_diab_Baseline",
                 "type2_diab_Baseline","visitdate_Baseline","visitdate_end",
                 "wgt_end","hba1c_end","bmikgm2_end","timeinfreezer_end","age",
                 "sex"))

setnames(bbs_main,
         old = c("wgt_36months","visit_date_36months","bmikgm2_36months",
                 "hba1c_36months","timeinfreezer_36months","type1_diab0",
                 "type2_diab0", "visit_date_Baseline","age_at_rand"),
         new = c("wgt_end","visitdate_end","bmikgm2_end","hba1c_end",
                 "timeinfreezer_end","type1_diab_Baseline",
                 "type2_diab_Baseline", "visitdate_Baseline","age"))

# make a version that keeps all the variables (with NA for the substudy vars if its not recorded there)
colnames(bbs_main)
colnames(bbs_substudy)
bbs_main <- bbs_main[,c("studyid","opdate","sex","age","ethnicity","hgt_m",
                        "type1_diab_Baseline","type2_diab_Baseline",
                        "diab_date_m","diab_date_y","wgt_Baseline","wgt_op_day",
                        "wgt_end","wgt_change_to_opday","visitdate_Baseline",
                        "visitdate_end","bmikgm2_Baseline","bmikgm2_end",
                        "diabetes_month_year","diabetes_duration",
                        "baseline_to_op","op_to_end","baseline_to_end",
                        "date_thawed","timeinfreezer_Baseline",
                        "timeinfreezer_end","hba1c_Baseline","hba1c_end","employ",
                        "income0","smoke","hdlc_Baseline","hdlc_36months",
                        "ldlc_Baseline","ldlc_36months","fastgluc_Baseline",
                        "fastgluc_36months","waistaverage_Baseline",
                        "waistaverage_36months","sbpaverage_Baseline",
                        "sbpaverage_36months","dbpaverage_Baseline",
                        "dbpaverage_36months","chol_Baseline","chol_36months")]
bbs_substudy <- bbs_substudy[,c("studyid","opdate","sex","age","ethnicity","hgt_m",
                                "type1_diab_Baseline","type2_diab_Baseline",
                                "diab_date_m","diab_date_y","wgt_Baseline","wgt_op_day",
                                "wgt_end","wgt_change_to_opday","visitdate_Baseline",
                                "visitdate_end","bmikgm2_Baseline","bmikgm2_end",
                                "diabetes_month_year","diabetes_duration",
                                "baseline_to_op","op_to_end","baseline_to_end",
                                "date_thawed","timeinfreezer_Baseline",
                                "timeinfreezer_end","hba1c_Baseline","hba1c_end")]

# add NA columns in the substudy df for the vars that are only collected in main
bbs_substudy$employ <- NA
bbs_substudy$income0 <- NA
bbs_substudy$smoke <- NA
bbs_substudy$hdlc_Baseline <- NA
bbs_substudy$hdlc_36months <- NA
bbs_substudy$ldlc_Baseline <- NA
bbs_substudy$ldlc_36months <- NA
bbs_substudy$fastgluc_Baseline <- NA
bbs_substudy$fastgluc_36months <- NA
bbs_substudy$waistaverage_Baseline <- NA
bbs_substudy$waistaverage_36months <- NA
bbs_substudy$sbpaverage_Baseline <- NA
bbs_substudy$sbpaverage_36months <- NA
bbs_substudy$dbpaverage_Baseline <- NA
bbs_substudy$dbpaverage_36months <- NA
bbs_substudy$chol_Baseline <- NA
bbs_substudy$chol_36months <- NA


# other version - only keep the vars they have in common
colnames(bbs_main)
colnames(bbs_substudy)
bbs_main_common <- bbs_main[,c("studyid","opdate","sex","age","ethnicity","hgt_m",
                               "type1_diab_Baseline","type2_diab_Baseline",
                               "diab_date_m","diab_date_y","wgt_Baseline","wgt_op_day",
                               "wgt_end","wgt_change_to_opday","visitdate_Baseline",
                               "visitdate_end","bmikgm2_Baseline","bmikgm2_end",
                               "diabetes_month_year","diabetes_duration",
                               "baseline_to_op","op_to_end","baseline_to_end",
                               "date_thawed","timeinfreezer_Baseline",
                               "timeinfreezer_end","hba1c_Baseline","hba1c_end")]
bbs_substudy_common <- bbs_substudy[,c("studyid","opdate","wgt_Baseline","wgt_op_day",
                                       "diab_date_m","diab_date_y","hba1c_Baseline",
                                       "type1_diab_Baseline","type2_diab_Baseline",
                                       "ethnicity","visitdate_Baseline","visitdate_end",
                                       "wgt_end","hba1c_end","hgt_m","bmikgm2_Baseline",
                                       "bmikgm2_end","diabetes_month_year",
                                       "diabetes_duration","wgt_change_to_opday",
                                       "baseline_to_op","op_to_end","baseline_to_end",
                                       "date_thawed","timeinfreezer_Baseline",
                                       "timeinfreezer_end","age","sex")]


# order the cols for rbind
colnames(bbs_main_common)
colnames(bbs_substudy_common)
bbs_substudy_common <- bbs_substudy_common[,c("studyid","opdate","sex","age",
                                "ethnicity","hgt_m","type1_diab_Baseline","type2_diab_Baseline",
                                "diab_date_m","diab_date_y","wgt_Baseline","wgt_op_day","wgt_end","wgt_change_to_opday",
                                "visitdate_Baseline","visitdate_end","bmikgm2_Baseline","bmikgm2_end",
                                "diabetes_month_year","diabetes_duration","baseline_to_op","op_to_end","baseline_to_end",
                                "date_thawed","timeinfreezer_Baseline","timeinfreezer_end","hba1c_Baseline","hba1c_end")]


###############################################################################
############################### MERGE STUDIES #################################
###############################################################################

## full version
# add a column indicating study
bbs_main$study <- "main"
bbs_substudy$study <- "substudy"

bbs_all_cols <- rbind(bbs_main,bbs_substudy)

# replace all '36months' with 'end'
colnames(bbs_all_cols) <- gsub("36months", "end", colnames(bbs_all_cols))

# now make a long version - as some of the vars apply to a specific timepoint/sample
long_vars <- c("wgt_Baseline","wgt_end","visitdate_Baseline","visitdate_end",
               "bmikgm2_Baseline","bmikgm2_end","timeinfreezer_Baseline",
               "timeinfreezer_end","hba1c_Baseline","hba1c_end","hdlc_Baseline",
               "hdlc_end","ldlc_Baseline","ldlc_end","fastgluc_Baseline",
               "fastgluc_end","waistaverage_Baseline","waistaverage_end",
               "sbpaverage_Baseline","sbpaverage_end","dbpaverage_Baseline",
               "dbpaverage_end","chol_Baseline","chol_end")

bbs_all_long <- pivot_longer(bbs_all_cols,
                             cols=all_of(long_vars),
                             names_to=c(".value","timepoint"),
                             names_sep = "_")

## common version
# add a column indicating study
bbs_main_common$study <- "main"
bbs_substudy_common$study <- "substudy"

bbs_all_common <- rbind(bbs_main_common,bbs_substudy_common)

# save out
write.csv(bbs_all_common, file = paste0(data_output_dir, "03_bbs_all_wide.csv"), row.names = F)

# now make a long version - as some of the vars apply to a specific timepoint/sample
long_vars <- c("wgt_Baseline", "wgt_end", "visitdate_Baseline", "visitdate_end", "bmikgm2_Baseline", "bmikgm2_end",
               "timeinfreezer_Baseline", "timeinfreezer_end", "hba1c_Baseline", "hba1c_end")

bbs_all_long_common <- pivot_longer(bbs_all_common,
                             cols=all_of(long_vars),
                             names_to=c(".value","timepoint"),
                             names_sep = "_")

###############################################################################
########################### DERIVE WEIGHT LOSS ################################
###############################################################################

# in the wide all cols version
bbs_all_cols$wgt_change_to_36m <- bbs_all_cols$wgt_end - bbs_all_cols$wgt_Baseline
# as a percentage of baseline
bbs_all_cols$wgt_change_to_36m_percentage <- (bbs_all_cols$wgt_change_to_36m/bbs_all_cols$wgt_Baseline)*100

# save out
write.csv(bbs_all_cols, file = paste0(data_output_dir, "03_bbs_all_wide_all_cols.csv"), row.names = F)

###############################################################################
############################### SELECT SAMPLES ################################
###############################################################################

## select the data for the ppts that we have metabolite data for

manifest$timepoint <- gsub("36 months", "end", manifest$timepoint)
manifest$timepoint <- gsub("bbs_metabolomics_substudy_oneyear", "end", manifest$timepoint)
manifest$timepoint <- gsub("bbs_metabolomics_substudy_baseline", "Baseline", manifest$timepoint)
manifest$id_timepoint <- paste0(manifest$studyId,"_",manifest$timepoint)

bbs_all_long$id_timepoint <- paste0(bbs_all_long$studyid,"_",bbs_all_long$timepoint)

# merge with manifest - clinical data will repeat where there is more than one sample for a ppt
merge <- merge(manifest,bbs_all_long)

sample_clinical_data <- merge

# save out
write.csv(sample_clinical_data, file=paste0(data_output_dir,"03_sample_clinical_data_all_cols.csv"), row.names = F)


bbs_all_long_common$id_timepoint <- paste0(bbs_all_long_common$studyid,"_",bbs_all_long_common$timepoint)

# merge with manifest - clinical data will repeat where there is more than one sample for a ppt
merge <- merge(manifest,bbs_all_long_common)

sample_clinical_data_common <- merge

# save out
write.csv(sample_clinical_data_common, file=paste0(data_output_dir,"03_sample_clinical_data.csv"), row.names = F)


# capture session info
print("Session information:")
sessionInfo()

# save output
sink()

