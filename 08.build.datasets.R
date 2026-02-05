### script to run final run in metaboprep ###
### remove n<5 

####################
###### SET UP ######
####################

# load libraries
library('metaboprep')
library('readxl')
library('kableExtra')
library('openxlsx')

## read in parameters
source("parameter_files/08.parameters.R")

## capture output
sink(file=paste0(log_output_dir, "08.build.datasets.log"),split=TRUE)

## record today's date
today = Sys.Date()
today = gsub("-","_",today)

###############################################################################
################################# READ IN DATA ################################
###############################################################################

# feature metadata
features <- as.data.frame(read.csv(file=paste0(data_dir,"/bbs_alspac/bbs_plus_alspac_no_outliers_2025_04_04_Filtered_feature_data.txt"),h=T, sep="\t"))
names(features)[1] <- "feature_id"

# combined dataset
samples <- as.data.frame(read.table(file=paste0(data_dir,"/bbs_alspac/bbs_plus_alspac_sample_data.txt"),h=T))
names(samples)[1] <- "sample_id"
input_dat <- read.table(file=paste0(data_dir,"/bbs_alspac/bbs_plus_alspac_metabolite_data.txt"),h=T)
data <- as.matrix(input_dat)

# ALSPAC dataset
alspac_samples <- as.data.frame(read.table(file=paste0(data_dir,"/alspac/alspac_only_sample_data.txt"),h=T))
names(alspac_samples)[1] <- "sample_id"
alspac_input_dat <- read.table(file=paste0(data_dir,"/alspac/alspac_only_metabolite_data.txt"),h=T)
alspac_data <- as.matrix(alspac_input_dat)

# BBS dataset
bbs_samples <- as.data.frame(read.table(file=paste0(data_dir,"/bbs/bbs_only_sample_data.txt"),h=T))
names(bbs_samples)[1] <- "sample_id"
bbs_input_dat <- read.table(file=paste0(data_dir,"/bbs/bbs_only_metabolite_data.txt"),h=T)
bbs_data <- as.matrix(bbs_input_dat)

###############################################################################
######################## IDENTIFY FEATURES TO DROP (N<5) ######################
###############################################################################

## process alspac
# create metaboprep object
alspac_mydata <- Metaboprep(data=alspac_data, 
                     features = features, 
                     samples  = alspac_samples)

# run feature summary
alspac_feature_sum <- feature_summary(metaboprep     = alspac_mydata, 
                                source_layer    = "input", 
                                outlier_udist   = 1.0,
                                tree_cut_height = 0.5,
                                output          = "data.frame")

# identify features with n<5 
alspac_features_to_remove <- alspac_feature_sum[alspac_feature_sum$n < 5, c("feature_id")]
length(alspac_features_to_remove)

###############################################################################

## process bbs
# create metaboprep object
bbs_mydata <- Metaboprep(data=bbs_data, 
                            features = features, 
                            samples  = bbs_samples)

# run feature summary
bbs_feature_sum <- feature_summary(metaboprep     = bbs_mydata, 
                                      source_layer    = "input", 
                                      outlier_udist   = 1.0,
                                      tree_cut_height = 0.5,
                                      output          = "data.frame")

# identify features with n<5 
bbs_features_to_remove <- bbs_feature_sum[bbs_feature_sum$n < 10, c("feature_id")]
length(bbs_features_to_remove)

###############################################################################

# combine lists
features_to_remove <- c(alspac_features_to_remove, bbs_features_to_remove)
length(features_to_remove)
features_to_remove <- unique(features_to_remove)
length(features_to_remove)

###############################################################################
################################ REMOVE FEATURES ##############################
###############################################################################

# create new metaboprep objects with reduced set of features

# all
cols_to_keep <- which(!names(input_dat) %in% features_to_remove)
data2 <- as.matrix(input_dat[,cols_to_keep])
features2 <- features[which(!features$feature_id %in% features_to_remove),]
  
all_mydata <- Metaboprep(data=data2, 
                     features = features2, 
                     samples  = samples)

# alspac
cols_to_keep <- which(!names(alspac_input_dat) %in% features_to_remove)
alspac_data2 <- as.matrix(alspac_input_dat[,cols_to_keep])

alspac_mydata2 <- Metaboprep(data=alspac_data2, 
                         features = features2, 
                         samples  = alspac_samples)

# bbs
cols_to_keep <- which(!names(bbs_input_dat) %in% features_to_remove)
bbs_data2 <- as.matrix(bbs_input_dat[,cols_to_keep])

bbs_mydata2 <- Metaboprep(data=bbs_data2, 
                             features = features2, 
                             samples  = bbs_samples)

###############################################################################
##################################### RUN QC ##################################
###############################################################################

## Identify the Xenobiotics to exclude from the QC steps
xenos <- all_mydata@features[which(all_mydata@features$SUPER_PATHWAY == "Xenobiotics"), "feature_id" ]

## how many xenobiotics identified
length(xenos)

## QC 

all_mydata <- all_mydata |>
  quality_control(source_layer        = "input", 
                  sample_missingness  = NA, 
                  feature_missingness = NA, 
                  total_peak_area_sd  = NA, 
                  outlier_udist       = NA, 
                  outlier_treatment   = "leave_be", 
                  winsorize_quantile  = 1.0, 
                  tree_cut_height     = 0.5, 
                  pc_outlier_sd       = NA,
                  feature_selection   = "max_var_exp", 
                  features_exclude_but_keep = xenos 
  )

alspac_mydata2 <- alspac_mydata2 |>
  quality_control(source_layer        = "input", 
                  sample_missingness  = NA, 
                  feature_missingness = NA, 
                  total_peak_area_sd  = NA, 
                  outlier_udist       = NA, 
                  outlier_treatment   = "leave_be", 
                  winsorize_quantile  = 1.0, 
                  tree_cut_height     = 0.5, 
                  pc_outlier_sd       = NA,
                  feature_selection   = "max_var_exp", 
                  features_exclude_but_keep = xenos 
  )

bbs_mydata2 <- bbs_mydata2 |>
  quality_control(source_layer        = "input", 
                  sample_missingness  = NA, 
                  feature_missingness = NA, 
                  total_peak_area_sd  = NA, 
                  outlier_udist       = NA, 
                  outlier_treatment   = "leave_be", 
                  winsorize_quantile  = 1.0, 
                  tree_cut_height     = 0.5, 
                  pc_outlier_sd       = NA,
                  feature_selection   = "max_var_exp", 
                  features_exclude_but_keep = xenos 
  )

###############################################################################
########################### GENERATE REPORTS ##################################
###############################################################################

## Generate the metaboprep report

# render reports
generate_report(all_mydata,
                project         = "bbs_alspac",
                output_dir      = paste0(data_output_dir, "bbs_alspac/metaboprep_output/"),
                output_filename = "bbs_alspac_metaboprep_qc_report.html",
                format          = "html",
                template        = "qc_report")

generate_report(bbs_mydata2,
                project         = "bbs",
                output_dir      = paste0(data_output_dir, "bbs/metaboprep_output/"),
                output_filename = "bbs_metaboprep_qc_report.html",
                format          = "html",
                template        = "qc_report")

generate_report(alspac_mydata2,
                project         = "alspac",
                output_dir      = paste0(data_output_dir, "alspac/metaboprep_output/"),
                output_filename = "alspac_metaboprep_qc_report.html",
                format          = "html",
                template        = "qc_report")


###############################################################################
######################### BUILD DATA DICTIONARIES #############################
###############################################################################

## file 1: data file
data_dictionary <- all_mydata@features[,c("feature_id", "CHEMICAL_NAME")]
data_dictionary$column_index <- seq(2,nrow(data_dictionary)+1,1)
names(data_dictionary)[1] <- "variable_name"
names(data_dictionary)[2] <- "variable_label"
data_dictionary <- data_dictionary[,c("column_index", "variable_name", "variable_label")]

# add row with sample ID
new_row <- c(1,"sample_id","Unique sample identifier (assigned by Metabolon)")
data_dictionary <- rbind(new_row,data_dictionary)

###############################################################################

## file 2: sample metadata
## merge sample metadata and sample summary 
# all samples
all_sample_summary <- as.data.frame(all_mydata@sample_summary[,,"input"])
all_sample_summary$sample_id <- rownames(all_sample_summary)
all_sampleinfo_merged <- merge(all_sample_summary,as.data.frame(all_mydata@samples), by="sample_id")
# drop pca vars
all_sampleinfo_merged <- all_sampleinfo_merged[,c(1:5,23:37)]

# alspac samples
alspac_sample_summary <- as.data.frame(alspac_mydata2@sample_summary[,,"input"])
alspac_sample_summary$sample_id <- rownames(alspac_sample_summary)
alspac_sampleinfo_merged <- merge(alspac_sample_summary,as.data.frame(alspac_mydata2@samples), by="sample_id")
# drop pca vars
alspac_sampleinfo_merged <- alspac_sampleinfo_merged[,c(1:5,23:37)]

# bbs samples
bbs_sample_summary <- as.data.frame(bbs_mydata2@sample_summary[,,"input"])
bbs_sample_summary$sample_id <- rownames(bbs_sample_summary)
bbs_sampleinfo_merged <- merge(bbs_sample_summary,as.data.frame(bbs_mydata2@samples), by="sample_id")
# drop pca vars
bbs_sampleinfo_merged <- bbs_sampleinfo_merged[,c(1:5,23:37)]

## build dictionaries
variable_labels <- c("Unique sample identifier (assigned by Metabolon)","Proportion of missing (NA) data (from metaboprep)", "total peak area (from metaboprep)",
                    "Total peak area based on complete features only (from metaboprep)", "Number of outlying values (+/- 50 IQR from median) (from metaboprep)",
                    "Batch identifiers for the neg platform", "Batch identifiers for the polar platform", "Batch identifiers for the pos.early platform",
                    "Batch identifiers for the pos.late platform", "Shipment box number", "Sample type", "Sample_barcode", "Run order","Plate identifier",
                    "Assigned plate number (36 per plate)", "Volume of material sent to Metabolon", "Units of volume", "Location in shipment box",
                    "Volume of material extracted by Metabolon for analysis", "Within plate order")

# all 
all_sample_dictionary <- as.data.frame(names(all_sampleinfo_merged))
names(all_sample_dictionary)[1] <- "variable_name"
all_sample_dictionary$column_index <- seq(1,nrow(all_sample_dictionary),1)
all_sample_dictionary$variable_label <- variable_labels
all_sample_dictionary <- all_sample_dictionary[,c("column_index", "variable_name", "variable_label")]

# alspac samples 
alspac_sample_dictionary <- as.data.frame(names(alspac_sampleinfo_merged))
names(alspac_sample_dictionary)[1] <- "variable_name"
alspac_sample_dictionary$column_index <- seq(1,nrow(alspac_sample_dictionary),1)
alspac_sample_dictionary$variable_label <- variable_labels
alspac_sample_dictionary <- alspac_sample_dictionary[,c("column_index", "variable_name", "variable_label")]

# bbs samples 
bbs_sample_dictionary <- as.data.frame(names(bbs_sampleinfo_merged))
names(bbs_sample_dictionary)[1] <- "variable_name"
bbs_sample_dictionary$column_index <- seq(1,nrow(bbs_sample_dictionary),1)
bbs_sample_dictionary$variable_label <- variable_labels
bbs_sample_dictionary <- bbs_sample_dictionary[,c("column_index", "variable_name", "variable_label")]

###############################################################################

## file 3: feature metadata
## merge feature metadata and feature summary 

# all samples
feature_summary <- as.data.frame(t(all_mydata@feature_summary[,,"input"]))
feature_summary$feature_id <- rownames(feature_summary)
feature_summary_merged <- merge(as.data.frame(all_mydata@features),feature_summary, by="feature_id")
# drop exclusion columns - not relevant for this summary run
feature_summary_merged <- feature_summary_merged[,c(1,15,2:14,16:18, 21:40)]

# alspac samples
alspac_feature_summary <- as.data.frame(t(alspac_mydata2@feature_summary[,,"input"]))
alspac_feature_summary$feature_id <- rownames(alspac_feature_summary)
alspac_feature_summary_merged <- merge(as.data.frame(all_mydata@features),alspac_feature_summary, by="feature_id")
# drop exclusion columns - not relevant for this summary run
alspac_feature_summary_merged <- alspac_feature_summary_merged[,c(1,15,2:14,16:18, 21:40)]

# bbs samples
bbs_feature_summary <- as.data.frame(t(bbs_mydata2@feature_summary[,,"input"]))
bbs_feature_summary$feature_id <- rownames(bbs_feature_summary)
bbs_feature_summary_merged <- merge(as.data.frame(all_mydata@features),bbs_feature_summary, by="feature_id")
# drop exclusion columns - not relevant for this summary run
bbs_feature_summary_merged <- bbs_feature_summary_merged[,c(1,15,2:14,16:18, 21:40)]

## build dictionary
feature_dictionary <-as.data.frame(names(feature_summary_merged))
names(feature_dictionary)[1] <- "variable_name"
feature_dictionary$column_index <- seq(1,nrow(feature_dictionary),1)
feature_dictionary$variable_label <- c("An internally-generated Metabolon compound identifier (concatenated with compid_); identifies metabolites in data file (Metabolon supplied)", "The name of the identified biochemical (Metabolon supplied)",
                                       "The general biochemical class assigned to the identified biochemical (Metabolon supplied)", "A more specific biochemical class assigned to the identified biochemical (Metabolon supplied)",
                                       "Identifier assigned by the NCBI and searchable in the PubChem database (Metabolon supplied)","Method used to generate data (Metabolon supplied)",
                                       "An alphanumeric compound identifier and link to compound information in KEGG","An alphanumeric compound identifier and link to compound information in HMDB",
                                       "Unique identifier for each biochemical (Metabolon supplied)", "Library reference (Metabolon supplied)",
                                       "Chromotography library reference (Metabolon supplied)","Metabolon-generated number intended for sorting compounds by SUPER PATHWAY and SUB PATHWAY",
                                       "Indicates whether a compound is named or unnamed/unknown (Metabolon supplied)", "IUPAC textual chemical identifier derived from INChI (Metabolon supplied)",
                                       "Simplified molecular-input line-entry system (SMILES) line notation string (Metabolon supplied)","A shorter name of the identified biochemical, convenient for plotting (Metabolon supplied)",
                                       "A unique numerical identifier assigned by the Chemical Abstracts Service (CAS) (Metabolon supplied)", "A numerical identifier as maintained in the ChemSpider database (Metabolon supplied)",
                                       "Proportion of missing (NA) data (metaboprep output)","Number of outlying values (+/- 50 IQR from median) (metaboprep output)",
                                       "Sample size (non-NA) (metaboprep output)", "Mean (metaboprep output)", "Standard deviation (metaboprep output)", 
                                       "Median (metaboprep output)","Minimimum value (metaboprep output)", "Maximum value (metaboprep output)", "Range (max - min) (metaboprep output)",
                                       "Skewness (metaboprep output)", "Kurtosis (metaboprep output)", "Standard error (metaboprep output)", "Count of NA (metaboprep output)", 
                                       "Variance (metaboprep output)", "Variance to mean ratio (metaboprep output)", "Coefficient of variation (metaboprep output)", 
                                       "Shaprio's W-statistic of normality on raw distribution (metaboprep output)", "Shaprio's W-statistic of normality on log10 distribution (metaboprep output)",
                                       "Cluster assignment (metaboprep output)", "Feature is a representative independent feature (1=yes, 0=no) (metaboprep output)")
feature_dictionary <- feature_dictionary[,c("column_index", "variable_name", "variable_label")]


###############################################################################
################################## SAVE OUT ###################################
###############################################################################

## save out data & dictionary

# combined dataset
metabdata <- as.data.frame(all_mydata@data[,,"input"])
metabdata$sample_id <- rownames(metabdata)
metabdata <- metabdata[,c(ncol(metabdata),1:(ncol(metabdata)-1))]
all_tabs <- list("MetaboliteData" = metabdata,"MetaboliteDataDictionary" = data_dictionary,
            "SampleMetadata" = alspac_sampleinfo_merged, "SampleMetadataDictionary" = all_sample_dictionary, 
            "FeatureMetadata" = feature_summary_merged, "FeatureMetadataDictionary" = feature_dictionary)

# write as excel
write.xlsx(all_tabs, file=paste0(data_output_dir,"bbs_alspac/metaboprep_output/",today,"_bbs_alspac_dataset.xlsx"))

# write as text files
filenames <- names(all_tabs)
for (i in 1:length(all_tabs)){
  outname <- paste0(data_output_dir,"bbs_alspac/metaboprep_output/",today, "_", filenames[i], ".txt")
  write.table(all_tabs[[i]], outname, quote= F, row.names = F, sep="\t")
}

# combined dataset - PUBLIC
all_tabs <- list("MetaboliteData" = NULL,"MetaboliteDataDictionary" = data_dictionary,
                 "SampleMetadata" = NULL, "SampleMetadataDictionary" = all_sample_dictionary, 
                 "FeatureMetadata" = feature_summary_merged, "FeatureMetadataDictionary" = feature_dictionary)

# write as excel
write.xlsx(all_tabs, file=paste0(data_output_dir,"bbs_alspac/metaboprep_output/",today,"_bbs_alspac_dataset_public.xlsx"))

# bbs only
metabdata <- as.data.frame(bbs_mydata2@data[,,"input"])
metabdata$sample_id <- rownames(metabdata)
metabdata <- metabdata[,c(ncol(metabdata),1:(ncol(metabdata)-1))]
bbs_tabs <- list("MetaboliteData" = metabdata,"MetaboliteDataDictionary" = data_dictionary,
                 "SampleMetadata" = bbs_sampleinfo_merged, "SampleMetadataDictionary" = bbs_sample_dictionary, 
                 "FeatureMetadata" = bbs_feature_summary_merged, "FeatureMetadataDictionary" = feature_dictionary)

# write as excel
write.xlsx(bbs_tabs, file=paste0(data_output_dir,"bbs/metaboprep_output/",today,"_bbs_dataset.xlsx"))

# write as text files
filenames <- names(bbs_tabs)
for (i in 1:length(bbs_tabs)){
  outname <- paste0(data_output_dir,"bbs/metaboprep_output/",today, "_", filenames[i], ".txt")
  write.table(bbs_tabs[[i]], outname, quote= F, row.names = F, sep="\t")
}

# bbs only - PUBLIC
bbs_tabs <- list("MetaboliteData" = NULL,"MetaboliteDataDictionary" = data_dictionary,
                 "SampleMetadata" = NULL, "SampleMetadataDictionary" = bbs_sample_dictionary, 
                 "FeatureMetadata" = bbs_feature_summary_merged, "FeatureMetadataDictionary" = feature_dictionary)

# write as excel
write.xlsx(bbs_tabs, file=paste0(data_output_dir,"bbs/metaboprep_output/",today,"_bbs_dataset_public.xlsx"))

# alspac only
metabdata <- as.data.frame(alspac_mydata2@data[,,"input"])
metabdata$sample_id <- rownames(metabdata)
metabdata <- metabdata[,c(ncol(metabdata),1:(ncol(metabdata)-1))]
alspac_tabs <- list("MetaboliteData" = metabdata,"MetaboliteDataDictionary" = data_dictionary,
                 "SampleMetadata" = alspac_sampleinfo_merged, "SampleMetadataDictionary" = alspac_sample_dictionary, 
                 "FeatureMetadata" = alspac_feature_summary_merged, "FeatureMetadataDictionary" = feature_dictionary)

# write as excel
write.xlsx(alspac_tabs, file=paste0(data_output_dir,"alspac/metaboprep_output/",today,"_alspac_dataset.xlsx"))

# write as text files
filenames <- names(alspac_tabs)
for (i in 1:length(alspac_tabs)){
  outname <- paste0(data_output_dir,"alspac/metaboprep_output/",today, "_", filenames[i], ".txt")
  write.table(alspac_tabs[[i]], outname, quote= F, row.names = F, sep="\t")
}

# alspac only - PUBLIC
alspac_tabs <- list("MetaboliteData" = NULL,"MetaboliteDataDictionary" = data_dictionary,
                    "SampleMetadata" = NULL, "SampleMetadataDictionary" = alspac_sample_dictionary, 
                    "FeatureMetadata" = alspac_feature_summary_merged, "FeatureMetadataDictionary" = feature_dictionary)

# write as excel
write.xlsx(alspac_tabs, file=paste0(data_output_dir,"alspac/metaboprep_output/",today,"_alspac_dataset_public.xlsx"))

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
#rm(list = ls())
