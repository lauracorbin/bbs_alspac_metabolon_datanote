### code to read in Metabolon raw data and export separate study versions
## run from script dir

####################
###### SET UP ######
####################

# load libraries
library(dplyr)
library(readxl)
library(purrr)
library(tidyverse)
# for writing out excel files
library(openxlsx)
library(magrittr)
library(ggplot2)

## read in parameters
source("parameter_files/03.parameters.R")

## capture output
sink(file=paste0(log_output_dir, "03.process.metabolon.log"),split=TRUE)

## record today's date
today = Sys.Date()
today = gsub("-","_",today)


####################################################################################
####################################################################################
## read in Metabolon data from excel file
####################################################################################
####################################################################################

## read in excel file 
inputfile = paste0(data_dir, file2process)
## read xls sheet names
sheetnames = readxl::excel_sheets( inputfile )

## read in all sheets
input_excel <- inputfile %>%
                    excel_sheets() %>%
                    purrr::set_names() %>%
                    map(read_excel, path = inputfile)

# read in manifest linker file
manifest <- read_delim(file=paste0(input_dir,"manifest_linker_file.csv"),col_names = T)
spec(manifest)

# read in report data from round 1 of metaboprep
load(paste0(metaboprep_rel, "ReportData.Rdata"))

####################################################################################
####################################################################################
## process files to split into studies
####################################################################################
####################################################################################

### extract sample IDs 
# for bbs
bbs_manifest <- manifest[manifest$source_study == "bbs_mainstudy" | manifest$source_study == "bbs_substudy",]
dim(bbs_manifest)
bbs_ids <- bbs_manifest$`Client Sample ID*`
length(bbs_ids)
length(unique(bbs_ids))

# for alspac
alspac_manifest <- manifest[manifest$source_study == "alspac",]
dim(alspac_manifest)
alspac_ids <- alspac_manifest$`Client Sample ID*`
length(alspac_ids)
length(unique(alspac_ids))

# combined list for testing
all_ids <- c(bbs_ids,alspac_ids)
length(all_ids)
length(unique(all_ids))

# identify samples to exclude based on PC values
## define tibble
pcs = as_tibble(qc_data$sample_data)

## Extract outliers
siteG_exc = as.data.frame(pcs %>% filter( PC1 > 0 & PC2 < c(-5) ) %>% pull(CLIENT_SAMPLE_ID))
dim(siteG_exc)
names(siteG_exc)[1] <- "CLIENT_SAMPLE_ID"

# joint dataset minus siteG exclusions
filtered_manifest <- manifest[!manifest$`Client Sample ID*` %in% siteG_exc$CLIENT_SAMPLE_ID,]
dim(filtered_manifest)
filtered_ids <- filtered_manifest$`Client Sample ID*`
length(filtered_ids)
length(unique(filtered_ids))

####################################################################################
####################################################################################

### make new sample meta data
# for bbs
bbs_sample_metadata <- input_excel$`Sample Meta Data`[input_excel$`Sample Meta Data`$CLIENT_SAMPLE_ID %in% bbs_ids,]
dim(bbs_sample_metadata)

# for alspac
alspac_sample_metadata <- input_excel$`Sample Meta Data`[input_excel$`Sample Meta Data`$CLIENT_SAMPLE_ID %in% alspac_ids,]
dim(alspac_sample_metadata)

# for filtered
filtered_sample_metadata <- input_excel$`Sample Meta Data`[input_excel$`Sample Meta Data`$CLIENT_SAMPLE_ID %in% filtered_ids,]
dim(filtered_sample_metadata)


### make new Peak Area Data
# for bbs
bbs_peak_area_data <- input_excel$`Peak Area Data`[input_excel$`Peak Area Data`$PARENT_SAMPLE_NAME %in% bbs_sample_metadata$PARENT_SAMPLE_NAME,]
dim(bbs_peak_area_data)

# for alspac
alspac_peak_area_data <- input_excel$`Peak Area Data`[input_excel$`Peak Area Data`$PARENT_SAMPLE_NAME %in% alspac_sample_metadata$PARENT_SAMPLE_NAME,]
dim(alspac_peak_area_data)

# for bbs filtered
filtered_peak_area_data <- input_excel$`Peak Area Data`[input_excel$`Peak Area Data`$PARENT_SAMPLE_NAME %in% filtered_sample_metadata$PARENT_SAMPLE_NAME,]
dim(filtered_peak_area_data)


### make new Batch-normalized Data
# for bbs
bbs_batch_norm_data <- input_excel$`Batch-normalized Data`[input_excel$`Batch-normalized Data`$PARENT_SAMPLE_NAME %in% bbs_sample_metadata$PARENT_SAMPLE_NAME,]
dim(bbs_batch_norm_data)

# for alspac
alspac_batch_norm_data <- input_excel$`Batch-normalized Data`[input_excel$`Batch-normalized Data`$PARENT_SAMPLE_NAME %in% alspac_sample_metadata$PARENT_SAMPLE_NAME,]
dim(alspac_batch_norm_data)

# for bbs filtered
filtered_batch_norm_data <- input_excel$`Batch-normalized Data`[input_excel$`Batch-normalized Data`$PARENT_SAMPLE_NAME %in% filtered_sample_metadata$PARENT_SAMPLE_NAME,]
dim(filtered_batch_norm_data)


### make new Batch-norm Imputed Data
# for bbs
bbs_batch_norm_imp_data <- input_excel$`Batch-norm Imputed Data`[input_excel$`Batch-norm Imputed Data`$PARENT_SAMPLE_NAME %in% bbs_sample_metadata$PARENT_SAMPLE_NAME,]
dim(bbs_batch_norm_imp_data)

# for alspac
alspac_batch_norm_imp_data <- input_excel$`Batch-norm Imputed Data`[input_excel$`Batch-norm Imputed Data`$PARENT_SAMPLE_NAME %in% alspac_sample_metadata$PARENT_SAMPLE_NAME,]
dim(alspac_batch_norm_imp_data)

# for bbs filtered
filtered_batch_norm_imp_data <- input_excel$`Batch-norm Imputed Data`[input_excel$`Batch-norm Imputed Data`$PARENT_SAMPLE_NAME %in% filtered_sample_metadata$PARENT_SAMPLE_NAME,]
dim(filtered_batch_norm_imp_data)


### make new Log Transformed Data
# for bbs
bbs_log_data <- input_excel$`Log Transformed Data`[input_excel$`Log Transformed Data`$PARENT_SAMPLE_NAME %in% bbs_sample_metadata$PARENT_SAMPLE_NAME,]
dim(bbs_log_data)

# for alspac
alspac_log_data <- input_excel$`Log Transformed Data`[input_excel$`Log Transformed Data`$PARENT_SAMPLE_NAME %in% alspac_sample_metadata$PARENT_SAMPLE_NAME,]
dim(alspac_log_data)

# for bbs filtered
filtered_log_data <- input_excel$`Log Transformed Data`[input_excel$`Log Transformed Data`$PARENT_SAMPLE_NAME %in% filtered_sample_metadata$PARENT_SAMPLE_NAME,]
dim(filtered_log_data)


### make new Vol_extracted-norm Data
# for bbs
bbs_vol_norm_data <- input_excel$`Vol_extracted-norm Data`[input_excel$`Vol_extracted-norm Data`$PARENT_SAMPLE_NAME %in% bbs_sample_metadata$PARENT_SAMPLE_NAME,]
dim(bbs_vol_norm_data)

# for alspac
alspac_vol_norm_data <- input_excel$`Vol_extracted-norm Data`[input_excel$`Vol_extracted-norm Data`$PARENT_SAMPLE_NAME %in% alspac_sample_metadata$PARENT_SAMPLE_NAME,]
dim(alspac_vol_norm_data)

# for bbs filtered
filtered_vol_norm_data <- input_excel$`Vol_extracted-norm Data`[input_excel$`Vol_extracted-norm Data`$PARENT_SAMPLE_NAME %in% filtered_sample_metadata$PARENT_SAMPLE_NAME,]
dim(filtered_vol_norm_data)

####################################################################################
####################################################################################
## build new excel files in Metabolon format - one per study and write out
####################################################################################
####################################################################################

# all by-band-sleeve
bbs <- list("Data Key & Explanation" = input_excel$`Data Key & Explanation`,"Chemical Annotation" = input_excel$`Chemical Annotation`,
            "Sample Meta Data"=bbs_sample_metadata, "Peak Area Data" = bbs_peak_area_data, "Batch-normalized Data" = bbs_batch_norm_data,
            "Batch-norm Imputed Data" = bbs_batch_norm_imp_data, "Log Transformed Data" = bbs_log_data, "Vol_extracted-norm Data" = bbs_vol_norm_data)

write.xlsx(bbs, file=paste0(bbs_output_raw,today,"_bbs_raw_metabolon_data.xlsx"))

# all alspac
alspac <- list("Data Key & Explanation" = input_excel$`Data Key & Explanation`,"Chemical Annotation" = input_excel$`Chemical Annotation`,
            "Sample Meta Data"=alspac_sample_metadata, "Peak Area Data" = alspac_peak_area_data, "Batch-normalized Data" = alspac_batch_norm_data,
            "Batch-norm Imputed Data" = alspac_batch_norm_imp_data, "Log Transformed Data" = alspac_log_data, "Vol_extracted-norm Data" = alspac_vol_norm_data)

write.xlsx(alspac, file=paste0(alspac_output_raw,today,"_alspac_raw_metabolon_data.xlsx"))
# also write out ID list for linking
write.table(alspac_ids,file=paste0(alspac_output_raw,today,"_alspac_metabolon_labids_B4132.csv"),row.names=F,quote = F,sep=",",col.names = F)

# combined bbs and alspac after erroneous samples removed
filtered <- list("Data Key & Explanation" = input_excel$`Data Key & Explanation`,"Chemical Annotation" = input_excel$`Chemical Annotation`,
            "Sample Meta Data"=filtered_sample_metadata, "Peak Area Data" = filtered_peak_area_data, "Batch-normalized Data" = filtered_batch_norm_data,
            "Batch-norm Imputed Data" = filtered_batch_norm_imp_data, "Log Transformed Data" = filtered_log_data, "Vol_extracted-norm Data" = filtered_vol_norm_data)

write.xlsx(filtered, file=paste0(filtered_output_raw,today,"_filtered_raw_metabolon_data.xlsx"))

# write out exclusion list
write.table(siteG_exc, file=paste0(data_output_dir, "intermediate/siteG_sample_outlier_anno.txt"), row.names = F, quote=F, sep="\t")


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

