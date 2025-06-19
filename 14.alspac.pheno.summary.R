### script to summarise alspac demographic, etc phenos for data note ###

####################
###### SET UP ######
####################

library(haven)
library(tidyverse)

source("parameter_files/14.parameters.R")

## capture output
sink(file=paste0(log_output_dir, "14.alspac.pheno.summary.log"),split=TRUE)

## record today's date
today = Sys.Date()
today = gsub("-","_",today)


###############################################################################
################################# READ IN DATA ################################
###############################################################################

# read in data extracted from R drive using script 11
dat_released <- read_dta(paste0(data_rel_input_dir,"11.alspac_pheno.dta"))
dim(dat_released)

# read in biosamples related vars provided by Kate N (since not released yet)
# provided as a stata file
dat_bespoke <- read_dta(paste0(data_bespoke_dir,"clinic_LauraC_18082023.dta")) 
dim(dat_bespoke) 

# read in the linked metabolite data (so can see which samples we have metabolon data for)
linked_dat <- read.csv(paste0(linked_data_dir,"Corbin_metabolomics_built_20250612.csv")) 
dim(linked_dat)

###############################################################################
################################# EXTRACT DATA ################################
###############################################################################

# make an alnqlet variable
dat_released$alnqlet <- paste0(dat_released$aln,dat_released$qlet)
dat_bespoke$alnqlet <- paste0(dat_bespoke$aln,dat_bespoke$qlet)
linked_dat$alnqlet <- paste0(linked_dat$aln,linked_dat$qlet)

# merge pheno
dat <- merge(dat_released, dat_bespoke, by="alnqlet", all.x=T)
dim(dat)

# reduce dat to only those we have metabolon data for
dat_red <- dat[which(dat$alnqlet%in%linked_dat$alnqlet),]
dim(dat_red)

# which participants dont have release data?
linked_dat[which(!linked_dat$alnqlet %in% dat$alnqlet),1:10] # should be none 

#######################################################################
########################### RUN SUMMARIES #############################
#######################################################################

# summaries for Table 2 in datanote
print("sex:")
table(dat_red$sex, useNA="always")
print("bmi")
summary(dat_red$bmi_kgm2)
sd(dat_red$bmi_kgm2, na.rm=T)
print("weight")
summary(dat_red$weight_kg)
sd(dat_red$weight_kg, na.rm=T)
print("age")
summary(dat_red$age_yrs)
sd(dat_red$age_yrs, na.rm=T)
print("ethnicity")
table(dat_red$ethnicity_cat, useNA="always") 
  ## 416 English/ Welsh/ Scottish/ Northern Irish/ British 
  ## 75 not completed
  ## (416/(519-75)=0.9369369)
print("insulin user at age 30")
table(dat_red$insulin_user, useNA="always") 

###############################################################################
################################## SAVE OUT ###################################
###############################################################################

# save out to use in 
write.csv(dat_red, file=paste0(data_output_dir, "14.alspac_pheno_processed.csv"), row.names = F)

####################################################################################
####################################################################################
## finish
####################################################################################
####################################################################################

# capture session info
print("Session information:")
sessionInfo()

# save output
sink()

# remove data
rm(list = ls())


