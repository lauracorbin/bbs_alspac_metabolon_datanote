### script to analyse the @30 glucose and cholesterol data (compare to Metabolon) for data note ###
### combine with the plots that have the BBS metab/clin comparisons

####################
###### SET UP ######
####################

# load libraries
library(tidyverse)
library(data.table)
library(blandr)

## read in parameters
source("parameter_files/15.parameters.R")

## capture output
sink(file=paste0(log_output_dir, "15.clin.vs.metabolon.log"),split=TRUE)

## record today's date
today = Sys.Date()
today = gsub("-","_",today)

###############################################################################
################################# READ IN DATA ################################
###############################################################################

# read in @30 cleaned data and metabolon data
data_30 <- read.csv(paste0(data_input_dir, "14.alspac_pheno_processed.csv"))
alspac <- read.csv(paste0(linked_data_dir, "Corbin_metabolomics_built_20250612.csv")) 
alspac$alnqlet <- paste0(alspac$aln,alspac$qlet)

# read in bbs data
load(paste0(bbs_dir, "intermediate/01_dataset_list.RData"))
bbs_metabolon <- dataset_list$raw

# check sample list of build dataset i am using matches latest
new_bbs_metabolon <- read.table(paste0(bbs_dir_updated),h=T,sep="\t")
length(which(new_bbs_metabolon$sample_id %in% bbs_metabolon$PARENT_SAMPLE_NAME))
length(which(!new_bbs_metabolon$sample_id %in% bbs_metabolon$PARENT_SAMPLE_NAME))

# read in clinical bloods
bbs_clinic <- read.csv(paste0(clinical_data_input_dir, "01_clinical_data_bbsmain_long.csv"))
# change "36 months" to "end"
bbs_clinic$sample <- gsub("36months", "end", bbs_clinic$sample)

###############################################################################
###################### APPLY WITHDRAWALS FOR ALSPAC ###########################
###############################################################################
# list of withdrawn aln extracted from ALSPAC-Data/Syntax/Withdrawal of consent/child_completed_WoC.do on 7th April 2025 (32 on list, none in our dataset)

# clinic data
dim(data_30)
data_30 <- data_30[!data_30$aln.x %in% withdrawals,]
dim(data_30)

# metabolite data - restrict to those with alspac pheno
dim(alspac)
alspac <- alspac[which(alspac$alnqlet %in% data_30$alnqlet),]
dim(alspac)

###############################################################################
########################## CORRELATIONS AND PLOTS #############################
###############################################################################


# merge bbs metabolon and clinic (only the cols with the bloods we want)
bbs_metabolon_clinic <- merge(bbs_metabolon, bbs_clinic[,c("sample","fastgluc","chol")], by.x = "id_timepoint", by.y = "sample")
# add a 'study' col
bbs_metabolon_clinic$study <- "BBS"

# merge alspac metabolon and clinic
alspac_metab_clinic <- merge(alspac, data_30[,c("alnqlet","Cholesterol","Glucose")], by = "alnqlet")
# rename cols to match bbs
setnames(alspac_metab_clinic,
          old = c("Cholesterol","Glucose"),
         new = c("chol","fastgluc"))
# add a 'study' col
alspac_metab_clinic$study <- "ALSPAC"

common_columns <- intersect(names(bbs_metabolon_clinic), names(alspac_metab_clinic))

mydata <- rbind(bbs_metabolon_clinic[common_columns], alspac_metab_clinic[common_columns])

# filter outliers as there are some values in the clinical chemistry that must be reporting errors
replace_outliers_with_na <- function(column) {
  mean_val <- mean(column, na.rm=T)
  sd_val <- sd(column, na.rm=T)
  lower_bound <- mean_val - 5 * sd_val
  upper_bound <- mean_val + 5 * sd_val
  column[column < lower_bound | column > upper_bound] <- NA
  return(column)
}

column_names <- c("fastgluc","chol")

# Loop over columns and apply the function to replace outliers with NA
for (col in column_names) {
  mydata[[col]] <- replace_outliers_with_na(mydata[[col]])
}

plot_list = list(glucose=NA, cholesterol=NA)
plot_list_ba = list(glucose=NA, cholesterol=NA)

cor_res_alspac <- cor.test(mydata[which(mydata$study=="ALSPAC"),"fastgluc"], mydata[which(mydata$study=="ALSPAC"),"compid_48152"], method = "pearson", use = "complete.obs")
cor_res_bbs <- cor.test(mydata[which(mydata$study=="BBS"),"fastgluc"], mydata[which(mydata$study=="BBS"),"compid_48152"], method = "pearson", use = "complete.obs")
plot_list[["glucose"]] <- ggplot(mydata, aes(x=fastgluc, y=compid_48152)) +
  geom_point(aes(x=fastgluc, y=compid_48152, color=study), shape=1) +
  ylim(0,4) +
  geom_smooth(data = mydata, aes(x = fastgluc, y = compid_48152, color = study),
              method = 'lm', formula = y ~ x, se = FALSE, linewidth=0.75) +
  xlab("clinical (mmol/L)") + ylab("Metabolon (scaled)") + ggtitle("A: Glucose") +
  theme(panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_line(color = "lightgray", linewidth = 0.2)) +
  annotate("text",x = 5, y = 3.6,
            label = paste0("ALSPAC r = ", round(cor_res_alspac$estimate,2), " (",round(cor_res_alspac$conf.int[1],2),", ",round(cor_res_alspac$conf.int[2],2),")"),
            size = 3.5) +
  annotate("text", x = 5, y = 3.2,
            label = paste0("BBS r = ", round(cor_res_bbs$estimate,2), " (",round(cor_res_bbs$conf.int[1],2),", ",round(cor_res_bbs$conf.int[2],2),")"),
            size = 3.5)

blandr_out <- blandr.statistics(scale(mydata[which(mydata$study=="ALSPAC"),"fastgluc"]), scale(mydata[which(mydata$study=="ALSPAC"),"compid_48152"]), sig.level = 0.95, LoA.mode = 1)
blandr_out
plot_list_ba[["glucose"]] <- blandr.plot.ggplot(blandr_out, method1name = "clinical", method2name = "MS", 
                                                   plotTitle= "A: Glucose")

cor_res_alspac <- cor.test(mydata[which(mydata$study=="ALSPAC"),"chol"], mydata[which(mydata$study=="ALSPAC"),"compid_63"], method = "pearson", use = "complete.obs")
cor_res_bbs <- cor.test(mydata[which(mydata$study=="BBS"),"chol"], mydata[which(mydata$study=="BBS"),"compid_63"], method = "pearson", use = "complete.obs")
plot_list[["cholesterol"]] <- ggplot(mydata, aes(x=chol, y=compid_63)) +
  geom_point(aes(x=chol, y=compid_63, color=study), shape=1) +
  ylim(0,4) +
  geom_smooth(data = mydata, aes(x = chol, y = compid_63, color = study),
              method = 'lm', formula = y ~ x, se = FALSE, linewidth=0.75) +
  xlab("clinical (mmol/L)") + ylab("Metabolon (scaled)") + ggtitle("B: Cholesterol") +
  theme(panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_line(color = "lightgray", linewidth = 0.2)) +
  annotate("text", x = 3.5, y = 3.6,
            label = paste0("ALSPAC r = ", round(cor_res_alspac$estimate,2), " (",round(cor_res_alspac$conf.int[1],2),", ",round(cor_res_alspac$conf.int[2],2),")"),
            size = 3.5) +
  annotate("text", x = 3.5, y = 3.2,
            label = paste0("BBS r = ", round(cor_res_bbs$estimate,2), " (",round(cor_res_bbs$conf.int[1],2),", ",round(cor_res_bbs$conf.int[2],2),")"),
            size = 3.5)

blandr_out <- blandr.statistics(scale(mydata[which(mydata$study=="ALSPAC"),"chol"]), scale(mydata[which(mydata$study=="ALSPAC"),"compid_63"]), sig.level = 0.95, LoA.mode = 1)
blandr_out
plot_list_ba[["cholesterol"]] <- blandr.plot.ggplot(blandr_out, method1name = "clinical", method2name = "MS", 
                                                plotTitle= "B: Cholesterol")
# save out plots and table
filename = paste0(fig_dir,"FigS1ab_assay_comparison.pdf")
outputfig <- ggpubr::ggarrange(plotlist=plot_list,
                               ncol=2, nrow=2,
                               common.legend = T,
                               legend="top",
                               heights = 2,
                               widths = 2)
ggpubr::ggexport(outputfig, filename=filename)

filename = paste0(fig_dir,"Fig3ab_assay_comparison.jpeg")
outputfig <- ggpubr::ggarrange(plotlist=plot_list_ba,
                               ncol=2, nrow=2,
                               common.legend = T,
                               legend="top",
                               heights = 2,
                               widths = 2)
ggpubr::ggexport(outputfig, filename=filename,dpi=600)

# get N numbers for each comparison
na <- colSums(is.na(mydata[,c("fastgluc","compid_48152","chol","compid_63")]))
na <- as.data.frame(na)

gluc_n <- nrow(mydata)-na["fastgluc",]-na["compid_48152",]
gluc_n
chol_n <- nrow(mydata)-na["chol",]-na["compid_63",]
chol_n


# capture session info
print("Session information:")
sessionInfo()

# save output
sink()

# remove data
rm(list = ls())


