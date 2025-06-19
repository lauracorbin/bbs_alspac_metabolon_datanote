### code to read in metaboprepped data (bbs+alspac) and export separate study versions

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
source("parameter_files/07.parameters.R")

## capture output
sink(file=paste0(log_output_dir, "07.split_data_by_study.log"),split=TRUE)

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


###############################################################################
############################## SPLIT BY STUDY #################################
###############################################################################

### extract sample IDs 
# for bbs
bbs_manifest <- manifest[manifest$source_study == "bbs_mainstudy" | manifest$source_study == "bbs_substudy",]
dim(bbs_manifest)
bbs_ids <- bbs_manifest$`Client.Sample.ID.`
length(bbs_ids)
length(unique(bbs_ids))

# for alspac
alspac_manifest <- manifest[manifest$source_study == "alspac",]
dim(alspac_manifest)
alspac_ids <- alspac_manifest$`Client.Sample.ID.`
length(alspac_ids)
length(unique(alspac_ids))


### make new sample meta data
# for bbs
bbs_sample_metadata <- sampledata[sampledata$CLIENT_SAMPLE_ID %in% bbs_ids,]
dim(bbs_sample_metadata) # 30 lost during filtering plus 142 site G outliers

# for alspac
alspac_sample_metadata <- sampledata[sampledata$CLIENT_SAMPLE_ID %in% alspac_ids,]
dim(alspac_sample_metadata) # 3 lost during filtering


### make new metabolite data
# for bbs
bbs_metabolite_data <- metabolitedata[rownames(metabolitedata) %in% bbs_sample_metadata$PARENT_SAMPLE_NAME,]
dim(bbs_metabolite_data)

# for alspac
alspac_metabolite_data <- metabolitedata[rownames(metabolitedata) %in% alspac_sample_metadata$PARENT_SAMPLE_NAME,]
dim(alspac_metabolite_data)


###############################################################################
################################# SAVE OUT ####################################
###############################################################################

write.table(bbs_metabolite_data, file=paste0(data_output_dir, "bbs/bbs_only_metabolite_data.txt"), sep = "\t", row.names = T, col.names = T)
write.table(alspac_metabolite_data, file=paste0(data_output_dir, "alspac/alspac_only_metabolite_data.txt"), sep = "\t", row.names = T, col.names = T)

write.table(bbs_sample_metadata, file=paste0(data_output_dir, "bbs/bbs_only_sample_data.txt"), sep = "\t", row.names = F, col.names = T)
write.table(alspac_sample_metadata, file=paste0(data_output_dir, "alspac/alspac_only_sample_data.txt"), sep = "\t", row.names = F, col.names = T)



# capture session info
print("Session information:")
sessionInfo()

# save output
sink()

# remove data
rm(list = ls())

