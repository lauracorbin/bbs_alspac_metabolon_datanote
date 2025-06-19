### script to read in and prepare participant clinical phenotype data for the BBS main study samples ###

####################
###### SET UP ######
####################

## read in parameter file
# parameter file
source("parameter_files/01.parameters.R")

## set working directory

## capture output
sink(file=paste0(log_output_dir, "01.clinical.data.prep.bbs.main.txt"),split=TRUE)

## load libraries
library(tidyverse)

###############################################################################
############################### READ IN DATA ##################################
###############################################################################

baseline <- read.table(paste0(bbs_main_data_dir, "metabolomics_sample_baseline_22JAN2024.csv"), sep = "|", header = TRUE, na.strings = ".")
followup <- read.table(paste0(bbs_main_data_dir, "metabolomics_sample_followup_22JAN2024.csv"), sep = "|", header = TRUE, na.strings = ".")
bloods <- read.table(paste0(bbs_main_data_dir, "metabolomics_sample_bloods_22JAN2024.csv"), sep = "|", header = TRUE, na.strings = ".")
meds <- read.table(paste0(bbs_main_data_dir, "metabolomics_sample_meds_22JAN2024.csv"), sep = "|", header = TRUE, na.strings = ".")

###############################################################################
################################# VIEW DATA ###################################
###############################################################################

print("baseline df dimensions:")
dim(baseline)
head(baseline)

print("followup df dimensions:")
dim(followup)
head(followup)

print("bloods df dimensions:")
dim(bloods)
head(bloods)

print("meds df dimensions:")
dim(meds)
head(meds)

###############################################################################
############################### PROCESS DATA ##################################
###############################################################################

### MEDS ###
## create a table of medication  categories - yes or no for each cat for each patient

# create the skeleton df
med_cats <- data.frame(matrix(NA, nrow = length(unique(meds$StudyID)), ncol = 22))
colnames(med_cats) <- c(1:11, 1:11)
colnames(med_cats)[1:11] <- paste("Baseline_", colnames(med_cats[,c(1:11)]), sep = "")
colnames(med_cats)[12:22] <- paste("36months_", colnames(med_cats[,c(12:22)]), sep = "")
rownames(med_cats) <- unique(meds$StudyID)

meds$Timepoint <- "Baseline"
w <- which(meds$visitid == 6)
meds$Timepoint[w] <- "36months"

# make wide so one column per timepoint
meds_wide <- pivot_wider(
  meds,
  id_cols = StudyID,
  names_from = c(Timepoint, drug_cat),
  names_prefix = "",
  names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "unique",
  values_from = c(drug_name),
  values_fill = NA,
  values_fn = NULL
)

# if there is 'NULL' in a cell for that category at that timepoint,
# then put 'No' in the corresponding cell of the med_cats df, or 'Yes' if not 'NULL'
cols <- c(paste0("Baseline_", c(1:5,7:11,NA)), paste0("36months_", c(1:5,7:11,NA)))
for (i in cols) {
  print(i)
  med_cats[,i] <- ifelse(meds_wide[,i]=="NULL", "No", "Yes")
}

# no category 6 drugs so add No to these columns
med_cats$Baseline_6 <- "No"
med_cats$`36months_6` <- "No"

## make long version of med cats
med_cats <- tibble::rownames_to_column(med_cats, var = "StudyID")
med_cats_long <- pivot_longer(
  med_cats,
  !StudyID,
  names_to=c("timepoint", ".value"),
  names_sep = "_"
)

#make a column of studyid + timepoint then delete originals
med_cats_long$sample <- paste(med_cats_long$StudyID, med_cats_long$timepoint, sep = "_")
med_cats_long <- subset(med_cats_long, select=-c(StudyID, timepoint))
colnames(med_cats_long)[1:12] <- paste0("drug_cat", colnames(med_cats_long)[1:12], sep = "")

### BLOODS ###

## Restrict bloods to timepoints of interest
bloods$Timepoint <- "Baseline"
w <- which(bloods$VisitID == 6)
bloods$Timepoint[w] <- "36months"
w <- which(is.na(bloods$VisitID))
bloods$Timepoint[w] <- NA

# make a column of studyid + timepoint for long merge (and remove original cols)
bloods$sample <- paste(bloods$StudyID, bloods$Timepoint, sep = "_")

# check for duplicates and remove
dim(unique(bloods))
length(unique(bloods$sample))
duplicates <- bloods[duplicated(bloods$sample),]
bloods <- unique(bloods)

duplicates <- bloods[duplicated(bloods$sample),] # still 9 samples duplicated but must have diff data, need to choose which to keep

dups <- which((bloods$StudyID=="REDACTED"|bloods$StudyID=="REDACTED"|bloods$StudyID=="REDACTED"|bloods$StudyID=="REDACTED"
              |bloods$StudyID=="REDACTED"|bloods$StudyID=="REDACTED"|bloods$StudyID=="REDACTED"|bloods$StudyID=="REDACTED"
              |bloods$StudyID=="REDACTED") & bloods$Timepoint=="36months")

view_dups <- bloods[dups,]
# only thing that differs for these is the ELF data - keep the one with the date that makes sense with the blood date (i.e. is a few days after it)
bloods <- bloods %>%
  arrange(StudyID, VisitID, desc(as.Date(ELF_date))) %>%
  group_by(StudyID, VisitID) %>%
  slice(1) %>%
  ungroup()

# remove those that have NA for timepoint
w <- which(is.na(bloods$VisitID))
bloods <- bloods[-w,]
table(bloods$VisitID, useNA = "ifany")

# Change blood into wide format
bloodswide <- pivot_wider(
  bloods,
  id_cols = StudyID,
  names_from = Timepoint,
  names_prefix = "",
  names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "unique",
  values_from = c(hb, plate, wbc, rbc, hct, mch, mcv, mchc, lymph, neutro, sodium, potass,
                  urea, creat, alt, alp, albumin, bili, bili_sign, hba1c, serumiron,
                  ferrit, vitb12, folate, folate_sign, hydroxy, calc, parathyroid, crp,
                  mag, phos, totalprot, blood_date, blood_time, chol, hdlc, ldlc, triglyc,
                  fastgluc, fblood_date, fblood_time, blood_liver, blood_liver_date, blood_liver_time,
                  ELF_Fibrosis, blood_research, blood_research_date, blood_research_time),
  values_fill = NA,
  values_fn = NULL
)

### FOLLOW-UP ###
## Restrict follow up to timepoint 0 and 6 and 1 (opdate)
w <- which(followup$visitid == 0 | followup$visitid == 6 | followup$visitid == 1)
followup <- followup[w,]

table(followup$visitid, useNA = "ifany")

followup$Timepoint <- "Baseline"
w <- which(followup$visitid == 6)
followup$Timepoint[w] <- "36months"
w <- which(followup$visitid == 1)
followup$Timepoint[w] <- "op_day"

# make wide
followup_wide <- pivot_wider(
  followup,
  id_cols = studyid,
  names_from = Timepoint,
  names_prefix = "",
  names_sep = "_",
  names_glue = NULL,
  names_sort = FALSE,
  names_repair = "unique",
  values_from = c(visit_date, wgt, preop_wgt_date, waist_1, waist_2, waist_3,
                  waist_4, sbp_1, dbp_1, sbp_2, dbp_3, sbp_3, dbp_3),
  values_fill = NA,
  values_fn = NULL
)

# make long follow up for baseline and endpoint
w <- which(followup$visitid == 0 | followup$visitid == 6)
followup_long <- followup[w,]
table(followup_long$visitid, useNA = "ifany")

followup_long$Timepoint <- "Baseline"
w <- which(followup_long$visitid == 6)
followup_long$Timepoint[w] <- "36months"

# make a column of studyid + timepoint
followup_long$sample <- paste(followup_long$studyid, followup_long$Timepoint, sep = "_")
followup_long <- subset(followup_long, select=-c(studyid, Timepoint))

### MERGE ALL WIDE DATA ###
mydata_prep <- merge(baseline, followup_wide, by="studyid")
mydata_prep2 <- merge(mydata_prep, bloodswide, by.x="studyid", by.y="StudyID")

## keeping meds seperate for now as not all of them have meds data so don't want to lose any data when merging
mydata <- mydata_prep2
dim(mydata)

###############################################################################
############################# DERIVE VARIABLES ################################
###############################################################################

### DERIVE BMI ###
mydata$hgt_m <- mydata$hgt/100
mydata$bmikgm2_Baseline <- mydata$wgt_Baseline / ((mydata$hgt_m)*(mydata$hgt_m)) 
mydata$bmikgm2_36months <- mydata$wgt_36months / ((mydata$hgt_m )*(mydata$hgt_m)) 

### DERIVE WAIST CIRCUMFERENCE ###
# make continuous then calculate average from repeats
mydata$waist_3_Baseline <- as.numeric(mydata$waist_3_Baseline)
mydata$waist_4_Baseline <- as.numeric(mydata$waist_4_Baseline)
mydata$waistaverage_Baseline <- rowMeans(mydata[, c("waist_1_Baseline", 
                                                     "waist_2_Baseline", "waist_3_Baseline", "waist_4_Baseline")], na.rm=TRUE )

mydata$waist_3_36months <- as.numeric(mydata$waist_3_36months)
mydata$waist_4_36months <- as.numeric(mydata$waist_4_36months)
mydata$waistaverage_36months <- rowMeans(mydata[, c("waist_1_36months", 
                                                      "waist_2_36months", "waist_3_36months", "waist_4_36months")], na.rm=TRUE )

### DERIVE BLOOD PRESSURE ###
# make bp vars continous
voi <- mydata[,c("sbp_1_Baseline", "sbp_1_36months", "dbp_1_Baseline", "dbp_1_36months",
                 "sbp_2_Baseline", "sbp_2_36months", "dbp_3_Baseline", "dbp_3_36months",
                 "sbp_3_Baseline", "sbp_3_36months")]

for (i in 1:10) {
  voi[,i] <- as.numeric(voi[,i])
}

## Derive average blood pressures from repeats
mydata$sbpaverage_Baseline <- rowMeans(voi[, c("sbp_1_Baseline", 
                                                "sbp_2_Baseline", "sbp_2_Baseline")], na.rm=TRUE )

mydata$sbpaverage_36months <- rowMeans(voi[, c("sbp_1_36months", 
                                                 "sbp_2_36months", "sbp_2_36months")], na.rm=TRUE ) ## 1 NA

mydata$dbpaverage_Baseline <- rowMeans(voi[, c("dbp_1_Baseline", 
                                                "dbp_3_Baseline")], na.rm=TRUE )

mydata$dbpaverage_36months <- rowMeans(voi[, c("dbp_1_36months", 
                                                 "dbp_3_36months")], na.rm=TRUE ) ## 1 NA

### DERIVE TIME SINCE DIABETES DIAGNOSIS ###
## combine month and year variables. 
## For ease, the start of each month is used as there is no day variable
# there are a couple of years/months that are non-sensical, and some ppts that have 'no indication' of diabetes but then have a start date - remove these
mydata[which(mydata$type2_diab0==0),c("diab_date_m","diab_date_y")] <- NA
mydata[which(mydata$diab_date_y<1950),c("diab_date_m","diab_date_y")] <- NA
mydata$diabetes_month_year <- na.omit(paste0("1/", mydata$diab_date_m, "/",  mydata$diab_date_y))
w <- which(mydata$diabetes_month_year == "1/NA/NA")
mydata$diabetes_month_year[w] <- NA
mydata$diabetes_month_year <- as.Date(mydata$diabetes_month_year, format="%d/%m/%Y")
mydata$visit_date_Baseline <- as.Date(mydata$visit_date_Baseline, format="%d/%m/%Y")

mydata$diabetes_duration <- difftime(mydata$visit_date_Baseline,mydata$diabetes_month_year)
mydata$diabetes_duration <- as.numeric(mydata$diabetes_duration)
mydata$diabetes_duration <- round(mydata$diabetes_duration*0.032855, digits = 0)

### DERIVE difference between weight @ baseline and weight @ op day ###
mydata$wgt_change_to_opday <- mydata$wgt_op_day - mydata$wgt_Baseline

### DERIVE OTHER TIME VARS ###
# is randomisation date the same as baseline visit date?
mydata$randdate <- as.Date(mydata$randdate, format="%d/%m/%Y")
length(which(mydata$randdate == mydata$visit_date_Baseline))
length(which(mydata$randdate != mydata$visit_date_Baseline))
w <- which(mydata$randdate != mydata$visit_date_Baseline)
diff <- mydata[w,c("randdate","visit_date_Baseline")]
dim(diff)

table(mydata$visit_date_Baseline, useNA = "ifany") # 27 NAs
table(mydata$randdate, useNA = "ifany") # 25 NAs

# is baseline blood date the same as baseline visit date?
mydata$blood_date_Baseline <- unlist(mydata$blood_date_Baseline)
mydata$blood_date_Baseline <- as.Date(mydata$blood_date_Baseline, format = "%d/%m/%Y")
length(which(mydata$visit_date_Baseline == mydata$blood_date_Baseline))
length(which(mydata$visit_date_Baseline != mydata$blood_date_Baseline))

table(mydata$blood_date_Baseline, useNA = "ifany") # 192 NAs

# Time from baseline sampling to surgery - baseline visit date to operation date
mydata$opdate <- as.Date(mydata$opdate, format = "%d/%m/%Y")
mydata$baseline_to_op <- difftime(mydata$opdate, mydata$visit_date_Baseline)
mydata$baseline_to_op <- as.numeric(mydata$baseline_to_op)
mydata$baseline_to_op <- round(mydata$baseline_to_op*0.032855, digits = 1)
min(mydata$baseline_to_op, na.rm = T) #0.2
max(mydata$baseline_to_op, na.rm = T) #42.8

# Time from operation to end sampling - op date to visit 6
mydata$visit_date_36months <- as.Date(mydata$visit_date_36months, format="%d/%m/%Y")
mydata$op_to_end <- difftime(mydata$visit_date_36months, mydata$opdate)
mydata$op_to_end <- as.numeric(mydata$op_to_end)
mydata$op_to_end <- round(mydata$op_to_end*0.032855, digits = 1)
min(mydata$op_to_end, na.rm = T) #-8.3
max(mydata$op_to_end, na.rm = T) #50.8

# Time from baseline sampling to end sampling - visit 0 to visit 6
mydata$baseline_to_end <- difftime(mydata$visit_date_36months, mydata$visit_date_Baseline)
mydata$baseline_to_end <- as.numeric(mydata$baseline_to_end)
mydata$baseline_to_end <- round(mydata$baseline_to_end*0.032855, digits = 1)
min(mydata$baseline_to_end, na.rm = T) #28.8
max(mydata$baseline_to_end, na.rm = T) #56.5

# Sample time in freezer
mydata$date_thawed <- as.Date("2023-01-10")
mydata$blood_research_date_Baseline <- unlist(mydata$blood_research_date_Baseline)
mydata$blood_research_date_Baseline <- as.Date(mydata$blood_research_date_Baseline, format = "%d/%m/%Y")
mydata$timeinfreezer_Baseline <- difftime(mydata$date_thawed, mydata$blood_research_date_Baseline)
mydata$fblood_date_36months <- gsub("NULL", NA, mydata$fblood_date_36months)
mydata$fblood_date_36months <- unlist(mydata$fblood_date_36months)
mydata$fblood_date_36months <- as.Date(mydata$fblood_date_36months, format = "%d/%m/%Y")
mydata$timeinfreezer_36months <- difftime(mydata$date_thawed, mydata$fblood_date_36months)
mydata$timeinfreezer_Baseline <- as.numeric(mydata$timeinfreezer_Baseline)
mydata$timeinfreezer_36months <- as.numeric(mydata$timeinfreezer_36months)
mydata$timeinfreezer_Baseline <- round(mydata$timeinfreezer_Baseline*0.032855, digits = 1)
mydata$timeinfreezer_36months <- round(mydata$timeinfreezer_36months*0.032855, digits = 1)
min(mydata$timeinfreezer_Baseline, na.rm = T) #38.6
max(mydata$timeinfreezer_Baseline, na.rm = T) #119.6
min(mydata$timeinfreezer_36months, na.rm = T) #3.7
max(mydata$timeinfreezer_36months, na.rm = T) #94.9

# Time difference between weight measures and sampling - visit date to blood research date (should mostly be the same)
mydata$visitdatetobloods_Baseline <- difftime(mydata$blood_research_date_Baseline, mydata$visit_date_Baseline, units = "days")
mydata$visitdatetobloods_36months <- difftime(mydata$fblood_date_36months, mydata$visit_date_36months, units = "days")
mydata$visitdatetobloods_Baseline <- as.numeric(mydata$visitdatetobloods_Baseline)
mydata$visitdatetobloods_36months <- as.numeric(mydata$visitdatetobloods_36months)
mydata$visitdatetobloods_Baseline <- round(mydata$visitdatetobloods_Baseline*0.032855, digits = 1)
mydata$visitdatetobloods_36months <- round(mydata$visitdatetobloods_36months*0.032855, digits = 1)

# collect derived variables for long format
derived_vars <- mydata[,c("studyid", "bmikgm2_Baseline", "bmikgm2_36months", "waistaverage_Baseline",
                          "waistaverage_36months", "sbpaverage_Baseline", "sbpaverage_36months",
                          "dbpaverage_Baseline", "dbpaverage_36months", "timeinfreezer_Baseline",
                          "timeinfreezer_36months", "visitdatetobloods_Baseline",
                          "visitdatetobloods_36months")]
# then pivot_longer
derived_vars_long <- pivot_longer(
  derived_vars,
  !studyid,
  names_to=c(".value", "timepoint"),
  names_sep = "_"
)

# and merge in other vars that have one value (not baseline and 36 months)
other_vars <- mydata[,c("studyid", "diabetes_duration", "wgt_change_to_opday", "baseline_to_op", "op_to_end", "baseline_to_end")]
all_derived_vars <- merge(derived_vars_long, other_vars, by = "studyid")
all_derived_vars$sample <- paste(all_derived_vars$studyid, all_derived_vars$timepoint, sep = "_")
all_derived_vars <- subset(all_derived_vars, select=-c(studyid, timepoint))


###############################################################################
################################### MERGE #####################################
###############################################################################

### merge all long data ###
mydata_prep <- merge(baseline, bloods, by.x="studyid", by.y="StudyID")
mydata_prep2 <- merge(mydata_prep, followup_long, by="sample")
mydata_prep3 <- merge(mydata_prep2, all_derived_vars, by="sample")
mydata_long <- merge(mydata_prep3, med_cats_long, by="sample")


###############################################################################
########################## REMOVE WITHDRAWALS #################################
###############################################################################

woc <- read.table(file = paste0(withdrawals, "WoC_20230608"))

mydata <- mydata[-(which(mydata$studyid%in%woc$V1)),]

mydata_long <- mydata_long[-(which(mydata_long$studyid%in%woc$V1)),]

###############################################################################
################################# SAVE OUT ####################################
###############################################################################

### SAVE MERGED DF ### (wide)
write.csv(mydata, file = paste0(data_output_dir, "01_clinical_data_bbsmain_wide.csv"), row.names = F)

#save long data
# version without the meds (as lose samples here)
write.csv(mydata_prep3, file = paste0(data_output_dir, "01_clinical_data_bbsmain_long.csv"), row.names = F)
write.csv(mydata_long, file = paste0(data_output_dir, "01_clinical_data_bbsmain_long_with_meds.csv"), row.names = F)


###############################################################################
############################# EXPLORE DISTRIBUTIONS ###########################
###############################################################################

### HISTOGRAMS - CONTINOUS VARS ###
## Subset all continuous variables and remove those that are just NAs
data_hist <- as.data.frame(lapply(mydata,as.numeric))
data_hist <- data_hist[,colSums(is.na(data_hist))<nrow(data_hist)]

## Remove all the categorical variables from histograms
data_hist <-
  data_hist[, -which(names(data_hist) %in%
                       c("type2_diab0", "BMI", "sex", "ethnicity", "employ",
                         "income0", "income_36m", "type1_diab0", "T2D.status",
                          "diab_date_m", "angina0", "prevmi0", "prevcabg0",
                         "prevstroke0", "prevclaud0", "nyha0",
                          "smoke", "visit_date_Baseline", "visit_data_36months", "folate_sign_36_months",
                         "visit_date_36months", "diab_date_y", "fblood_date_36months", "blood_liver_36months",
                         "blood_research_36months",
                         "diabetes_month_year", "randdate", "opdate", "visit_date_36_months", "blood_date_Baseline",
                         "fblood_date_36_months", "blood_liver_Baseline", "blood_liver_36_months",
                         "blood_research_Baseline", "blood_research_36_months", "blood_research_date_Baseline",
                         "date_thawed"))]

pdf(file = paste0(results_fig_dir, "01_clinic_vars_bbsmain_histograms.pdf"), height = 10, width = 7)
par(mfrow=c(3,2))
for (i in 1:108) {
  DataDescribed = psych::describe(data_hist[,i]) 
  meanvar<-DataDescribed$mean
  medianvar<-DataDescribed$median
  minvar<-DataDescribed$min
  maxvar<-DataDescribed$max
  kurtosisvar<-DataDescribed$kurtosis
  skewnessvar<-DataDescribed$skew
  N<- nrow(data_hist) - sum(is.na(data_hist[,i]))
  missingness <- (sum(is.na(data_hist[,i]))/nrow(data_hist))*100
  
  a<-density(data_hist[,i], na.rm=T)
  thresholdx<-(maxvar+(maxvar/100))
  thresholdy<-min(a$y)+(max(a$y)/4)
  
  hist(data_hist[,i], col="red",main=(names(data_hist)[i]),prob=TRUE,xlab="peak area") 
  lines(density(data_hist[,i], na.rm = TRUE),col="blue", lwd=2)
  text(thresholdx,thresholdy, cex=0.6, 
       paste("N=", N, "\npercent missing=", 
             signif(missingness, 3), "\nmin=", 
             signif(minvar, 3), " \nmax=",
             signif(maxvar, 3), 
             "\nmean=", 
             signif(meanvar, 3), " \nmedian=", signif(medianvar, 3), 
             "\nkurt=", 
             signif(kurtosisvar, 3), " \nskew=", 
             signif(skewnessvar, 3), sep = ''), pos = 3,xpd = NA)
}
dev.off() 


# capture session info
print("Session information:")
sessionInfo()

# save output
sink()


