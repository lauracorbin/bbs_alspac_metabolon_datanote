### script to remove the withdrawal of consents from the post-metaboprepped joint ALSPAC and BBS filtered data files ###

####################
###### SET UP ######
####################

# load libraries
library(dplyr)
library(readxl)
library(purrr)
library(tidyverse)
library(magrittr)
library(ggplot2)

## read in parameters
source("parameter_files/06.parameters.R")

## capture output
sink(file=paste0(log_output_dir, "06.remove_woc_from_joint_data.log"),split=TRUE)

## record today's date
today = Sys.Date()
today = gsub("-","_",today)


###############################################################################
############################### READ IN DATA ##################################
###############################################################################

# read in filtered data from metaboprep to be separated by study
metabolitedata <- read.csv(paste0(data_dir, "bbs_plus_alspac_no_outliers_2025_04_04_Filtered_metabolite_data.txt"), sep = "\t", header = T)
sampledata <- read.csv(paste0(data_dir, "bbs_plus_alspac_no_outliers_2025_04_04_Filtered_sample_data.txt"), sep = "\t", header = T)


# read in manifest linker file
manifest <- read.csv(file=paste0(manifest_dir,"manifest_linker_file.csv"))


###############################################################################
########################### REMOVE WITHDRAWALS ################################
###############################################################################

# remove wtihdrawals from manifest so they will be removed for the metabolite data on merging
woc <- read.table(file = paste0(withdrawals, "WoC_20230608"))

manifest <- manifest[-(which(manifest$studyId%in%woc$V1)),]

# sample data
sampledata <- sampledata[which(sampledata$CLIENT_SAMPLE_ID%in%manifest$Client.Sample.ID.),]

# metabolite data
metabolitedata <- metabolitedata[which(rownames(metabolitedata)%in%sampledata$PARENT_SAMPLE_NAME),]


###############################################################################
################################# SAVE OUT ####################################
###############################################################################

write.table(metabolitedata, file=paste0(data_output_dir, "bbs_alspac/bbs_plus_alspac_metabolite_data.txt"), sep = "\t", row.names = T, col.names = T)
write.table(sampledata, file=paste0(data_output_dir, "bbs_alspac/bbs_plus_alspac_sample_data.txt"), sep = "\t", row.names = F, col.names = T)



# capture session info
print("Session information:")
sessionInfo()

# save output
sink()

# remove data
rm(list = ls())

