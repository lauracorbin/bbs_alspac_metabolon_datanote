### script to read in and prepare participant clinical phenotype data for the BBS sub-study samples ###

####################
###### SET UP ######
####################

## read in parameter file
# parameter file
source("parameter_files/02.parameters.R")

## capture output
sink(file=paste0(log_output_dir, "02.clinical.data.prep.bbs.substudy.txt"),split=TRUE)

## load libraries
library(tidyverse)

###############################################################################
############################### READ IN DATA ##################################
###############################################################################

codebook <- read.table(paste0(bbs_substudy_data_dir, "Export_CodeBook.csv"), sep = "|", header = TRUE, na.strings = ".")
baseline <- read.table(paste0(bbs_substudy_data_dir, "Export_MET2.csv"), sep = "|", header = TRUE, na.strings = ".")
sample_baseline <- read.table(paste0(bbs_substudy_data_dir, "Export_MET3.csv"), sep = "|", header = TRUE, na.strings = ".")
sample_end <- read.table(paste0(bbs_substudy_data_dir, "Export_MET4.csv"), sep = "|", header = TRUE, na.strings = ".")
ethnicity <- read.table(paste0(bbs_substudy_data_dir, "Export_MET1_restricted_2023-04-03.csv"), sep = "|", header = TRUE, na.strings = ".")
withdraw <- read.table(paste0(bbs_substudy_data_dir, "Export_W.csv"), sep = "|", header = TRUE, na.strings = ".")
notes <- read.table(paste0(bbs_substudy_data_dir, "Export_N.csv"), sep = "|", header = TRUE, na.strings = ".")

###############################################################################
################################# VIEW DATA ###################################
###############################################################################

print("baseline (MET2) df dimensions:")
dim(baseline)
head(baseline)

print("sample_baseline (MET3) df dimensions:")
dim(sample_baseline)
head(sample_baseline)

print("sample_end (MET4) df dimensions:")
dim(sample_end)
head(sample_end)

print("ethnicity (MET1) df dimensions:")
dim(ethnicity)
head(ethnicity)

###############################################################################
############################### PROCESS DATA ##################################
###############################################################################

# select which columns to keep from each dataframe
baseline <- baseline[, c("studyId", "opdate", "baseline_wgt", "baseline_hgt", "op_wgt", "op_type", "diab_onset_m", "diab_onset_y",
                     "diab_hba1c_lev", "diab_hba1c_date", "diab_typeone", "diab_typetwo")]
ethnicity <- ethnicity[,c("studyId", "ethnicity", "age_at_consent_years", "ppt_gender")]
sample_baseline <- sample_baseline[,c("studyId", "blood_date_base", "blood_time_base", "last_food_base")]
sample_end <- sample_end[,c("studyId", "blood_date_m12", "blood_time_m12", "last_food_m12", "wgt_m12", "diab_typetwo_12m", "diab_hba1c_12m_lev",
                            "diab_hba1c_12m_date")]

# merge to make wide dataset
mydata_prep <- merge(baseline, ethnicity, all = T)
mydata_prep2 <- merge(mydata_prep, sample_baseline, all = T)
mydata <- merge(mydata_prep2, sample_end, all = T)


###############################################################################
############################# DERIVE VARIABLES ################################
###############################################################################

### DERIVE BMI ###
mydata$hgt_m <- mydata$baseline_hgt/100
mydata$bmikgm2_Baseline <- mydata$baseline_wgt / ((mydata$hgt_m)*(mydata$hgt_m)) 
mydata$bmikgm2_36months <- mydata$wgt_m12 / ((mydata$hgt_m )*(mydata$hgt_m)) 

### DERIVE TIME SINCE DIABETES DIAGNOSIS ###
## combine month and year variables. 
## For ease, the start of each month is used as there is no day variable
# there are a couple of years/months that are non-sensical, and some ppts that have 'no indication' of diabetes but then have a start date - remove these
mydata[which(mydata$type2_diab0==0),c("diab_onset_m","diab_onset_y")] <- NA
mydata[which(mydata$diab_onset_y<1950),c("diab_onset_m","diab_onset_y")] <- NA
mydata$diabetes_month_year <- na.omit(paste0("1/", mydata$diab_onset_m, "/",  mydata$diab_onset_y))
w <- which(mydata$diabetes_month_year == "1/NA/NA")
mydata$diabetes_month_year[w] <- NA
mydata$diabetes_month_year <- as.Date(mydata$diabetes_month_year, format="%d/%m/%Y")
mydata$blood_date_base <- as.Date(mydata$blood_date_base, format="%d/%m/%Y")

mydata$diabetes_duration <- difftime(mydata$blood_date_base,mydata$diabetes_month_year)
mydata$diabetes_duration <- as.numeric(mydata$diabetes_duration)
mydata$diabetes_duration <- round(mydata$diabetes_duration*0.032855, digits = 0)

### DERIVE difference between weight @ baseline and weight @ op day ###
mydata$wgt_change_to_opday <- mydata$op_wgt - mydata$baseline_wgt

# Time from baseline sampling to surgery - baseline visit date to operation date
mydata$opdate <- as.Date(mydata$opdate, format = "%d/%m/%Y")
mydata$baseline_to_op <- difftime(mydata$opdate, mydata$blood_date_base, units = "days")
mydata$baseline_to_op <- as.numeric(mydata$baseline_to_op)
mydata$baseline_to_op <- round(mydata$baseline_to_op*0.032855, digits = 1)
min(mydata$baseline_to_op, na.rm = T) #0
max(mydata$baseline_to_op, na.rm = T) #27.7

# Time from operation to end sampling - op date to visit 6
mydata$blood_date_m12 <- as.Date(mydata$blood_date_m12, format="%d/%m/%Y")
mydata$op_to_end <- difftime(mydata$blood_date_m12, mydata$opdate, units = "days")
mydata$op_to_end <- as.numeric(mydata$op_to_end)
mydata$op_to_end <- round(mydata$op_to_end*0.032855, digits = 1)
min(mydata$op_to_end, na.rm = T) #8.4
max(mydata$op_to_end, na.rm = T) #35.5

# Time from baseline sampling to end sampling - visit 0 to visit 6
mydata$baseline_to_end <- difftime(mydata$blood_date_m12, mydata$blood_date_base)
mydata$baseline_to_end <- as.numeric(mydata$baseline_to_end)
mydata$baseline_to_end <- round(mydata$baseline_to_end*0.032855, digits = 1)
min(mydata$baseline_to_end, na.rm = T) #12.6
max(mydata$baseline_to_end, na.rm = T) #39.4

# Sample time in freezer
mydata$date_thawed <- as.Date("2023-01-10")
mydata$timeinfreezer_Baseline <- difftime(mydata$date_thawed, mydata$blood_date_base, units = "days")
mydata$timeinfreezer_36months <- difftime(mydata$date_thawed, mydata$blood_date_m12, units = "days")
mydata$timeinfreezer_Baseline <- as.numeric(mydata$timeinfreezer_Baseline)
mydata$timeinfreezer_36months <- as.numeric(mydata$timeinfreezer_36months)
mydata$timeinfreezer_Baseline <- round(mydata$timeinfreezer_Baseline*0.032855, digits = 1)
mydata$timeinfreezer_36months <- round(mydata$timeinfreezer_36months*0.032855, digits = 1)
min(mydata$timeinfreezer_Baseline, na.rm = T) #8.4
max(mydata$timeinfreezer_Baseline, na.rm = T) #39.3
min(mydata$timeinfreezer_36months, na.rm = T) #-1.2
max(mydata$timeinfreezer_36months, na.rm = T) #23.5


## QC on age - some of them don't make sense due to errors in recording DOB
# set any age <18 to missing
hist(as.numeric(mydata$age_at_consent_years))
mydata$age_at_consent_years <- as.numeric(mydata$age_at_consent_years)
w <- which(mydata$age_at_consent_years<18)
mydata[w,"age_at_consent_years"] <- NA
hist(mydata$age_at_consent_years)

###############################################################################
################################# SAVE OUT ####################################
###############################################################################

write.csv(mydata, file = paste0(data_output_dir, "02_clinical_data_bbssubstudy_wide.csv"), row.names = F)


# capture session info
print("Session information:")
sessionInfo()

# save output
sink()
