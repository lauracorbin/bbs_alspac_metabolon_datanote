### script to merge the clinical data from the bbs main and sub studies, so it is ready to merge with the nmr metabolite data ###

####################
###### SET UP ######
####################

## read in parameter file
# parameter file
source("parameter_files/03a.parameters.R")

## capture output
sink(file=paste0(log_output_dir, "03a.merge.bbs.nmr.clinical.data.txt"),split=TRUE)

## load libraries
library(tidyverse)
library(data.table)

###############################################################################
############################### READ IN DATA ##################################
###############################################################################

# read in the processed data from main and sub studies
bbs_main <- read.csv(paste0(data_output_dir, "01_clinical_data_bbsmain_wide.csv"))
bbs_substudy <- read.csv(paste0(data_output_dir, "02_clinical_data_bbssubstudy_wide.csv"))

# read in manifest
manifest <- read.csv(paste0(inventory_dir, "2023Q2_nmr_annotated_manifest_2023-08-14.csv"))

###############################################################################
################################ CHOOSE VARS ##################################
###############################################################################

# only keep the vars they have in common, and standardise the column names
colnames(bbs_main)
colnames(bbs_substudy)
bbs_main <- bbs_main[,c("studyid","opdate","sex","age_at_rand","ethnicity","hgt_m","type1_diab0","type2_diab0",
                        "diab_date_m","diab_date_y","wgt_Baseline","wgt_op_day","wgt_36months","wgt_change_to_opday",
                        "visit_date_Baseline","visit_date_36months","bmikgm2_Baseline","bmikgm2_36months",
                        "diabetes_month_year","diabetes_duration","baseline_to_op","op_to_end","baseline_to_end","date_thawed",
                        "timeinfreezer_Baseline","timeinfreezer_36months","hba1c_Baseline","hba1c_36months")]
bbs_substudy <- bbs_substudy[,c("studyId","opdate","baseline_wgt","op_wgt","diab_onset_m","diab_onset_y","diab_hba1c_lev",
                                "diab_typeone","diab_typetwo","ethnicity","blood_date_base","blood_date_m12","wgt_m12",
                                "diab_hba1c_12m_lev","hgt_m","bmikgm2_Baseline","bmikgm2_36months","diabetes_month_year","diabetes_duration",
                                "wgt_change_to_opday","baseline_to_op","op_to_end","baseline_to_end","date_thawed","timeinfreezer_Baseline",
                                "timeinfreezer_36months","age_at_consent_years","ppt_gender")]

setnames(bbs_substudy,
         old = c("studyId","baseline_wgt","op_wgt","diab_onset_m","diab_onset_y","diab_hba1c_lev",
                               "diab_typeone","diab_typetwo","blood_date_base","blood_date_m12","wgt_m12",
                               "diab_hba1c_12m_lev","bmikgm2_36months","timeinfreezer_36months","age_at_consent_years","ppt_gender"),
         new = c("studyid","wgt_Baseline","wgt_op_day","diab_date_m","diab_date_y","hba1c_Baseline","type1_diab_Baseline","type2_diab_Baseline",
                 "visitdate_Baseline","visitdate_end","wgt_end","hba1c_end","bmikgm2_end","timeinfreezer_end","age","sex"))

setnames(bbs_main,
         old = c("wgt_36months","visit_date_36months","bmikgm2_36months","hba1c_36months","timeinfreezer_36months",
                 "type1_diab0","type2_diab0", "visit_date_Baseline","age_at_rand"),
         new = c("wgt_end","visitdate_end","bmikgm2_end","hba1c_end","timeinfreezer_end","type1_diab_Baseline",
                 "type2_diab_Baseline", "visitdate_Baseline","age"))


# order the cols for rbind
colnames(bbs_main)
colnames(bbs_substudy)
bbs_substudy <- bbs_substudy[,c("studyid","opdate","sex","age",
                                "ethnicity","hgt_m","type1_diab_Baseline","type2_diab_Baseline",
                                "diab_date_m","diab_date_y","wgt_Baseline","wgt_op_day","wgt_end","wgt_change_to_opday",
                                "visitdate_Baseline","visitdate_end","bmikgm2_Baseline","bmikgm2_end",
                                "diabetes_month_year","diabetes_duration","baseline_to_op","op_to_end","baseline_to_end",
                                "date_thawed","timeinfreezer_Baseline","timeinfreezer_end","hba1c_Baseline","hba1c_end")]


###############################################################################
############################### MERGE STUDIES #################################
###############################################################################

# add a column indicating study
bbs_main$study <- "main"
bbs_substudy$study <- "substudy"

bbs_all <- rbind(bbs_main,bbs_substudy) # do this once i have age and sex for substudy

# save out
write.csv(bbs_all, file = paste0(data_output_dir, "03_bbs_all_wide.csv"), row.names = F)

# now make a long version - as some of the vars apply to a specific timepoint/sample
long_vars <- c("wgt_Baseline", "wgt_end", "visitdate_Baseline", "visitdate_end", "bmikgm2_Baseline", "bmikgm2_end",
               "timeinfreezer_Baseline", "timeinfreezer_end", "hba1c_Baseline", "hba1c_end")

bbs_all_long <- pivot_longer(bbs_all,
                             cols=all_of(long_vars),
                             names_to=c(".value","timepoint"),
                             names_sep = "_")

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
write.csv(sample_clinical_data, file=paste0(data_output_dir,"03a_sample_clinical_data_all_nmr.csv"), row.names = F)


# capture session info
print("Session information:")
sessionInfo()

# save output
sink()

