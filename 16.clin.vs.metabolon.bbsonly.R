### this script is to compare clinical assay results to metabolon data
### created by Maddy Smith 

####################
###### SET UP ######
####################

## read in parameter file
source("parameter_files/16.parameters.R")

## capture output
sink(file=paste0(log_output_dir, "16.clin.vs.metabolon.bbsonly.log"),split=TRUE)

## load libraries
library(readxl)
library(tidyr)
library(dplyr)
library(lme4)
library(nlme)
library(ggplot2)
library(tidyverse)
library(gtools)
library(qqman)
library(lmtest)
library(stats)
library(blandr)

##############################################################################
################################ READ IN DATA ################################
##############################################################################

# read in dataset list
load(paste0(bbs_dir, "intermediate/01_dataset_list.RData"))

##############################################################################
###################### COMPARE CLINIC TO METABOLON ###########################
##############################################################################

metabolon <- dataset_list$raw

# read in clinical bloods
clinic <- read.csv(paste0(clinical_data_input_dir, "01_clinical_data_bbsmain_long.csv"))
# change "36 months" to "end"
clinic$sample <- gsub("36months", "end", clinic$sample)

# merge into dataset (only the cols with the bloods we want)
metabolon_clinic <- merge(metabolon, clinic[,c("sample","urea","creat","bili")], by.x = "id_timepoint", by.y = "sample")

# which samples aren't in the clinical data
no_clin <- metabolon[which(!metabolon$id_timepoint %in% clinic$sample),]

# filter outliers as there are some values in the clinical chemistry that must be reporting errors
replace_outliers_with_na <- function(column) {
  mean_val <- mean(column, na.rm=T)
  sd_val <- sd(column, na.rm=T)
  lower_bound <- mean_val - 5 * sd_val
  upper_bound <- mean_val + 5 * sd_val
  column[column < lower_bound | column > upper_bound] <- NA
  return(column)
}

column_names <- c("urea", "creat", "bili")

# Loop over columns and apply the function to replace outliers with NA
for (col in column_names) {
  metabolon_clinic[[col]] <- replace_outliers_with_na(metabolon_clinic[[col]])
}

## now recreate scatter plots
plot_list = list(urea=NA, creatinine=NA, bilirubin=NA)
plot_list_ba = list(urea=NA, creatinine=NA, bilirubin=NA)

## correlations
correlations <- as.data.frame(matrix(nrow=5, ncol = 2), row.names = names(plot_list))
colnames(correlations) <- c("pearson",  "spearman")
correlations["urea","pearson"] <- cor(metabolon_clinic$urea, metabolon_clinic$compid_1670, method = "pearson", use = "complete.obs")
correlations["urea","spearman"] <- cor(metabolon_clinic$urea, metabolon_clinic$compid_1670, method = "spearman", use = "complete.obs")
correlations["creatinine","pearson"] <- cor(metabolon_clinic$creat, metabolon_clinic$compid_513, method = "pearson", use = "complete.obs")
correlations["creatinine","spearman"] <- cor(metabolon_clinic$creat, metabolon_clinic$compid_513, method = "spearman", use = "complete.obs")
correlations["bilirubin","pearson"] <- cor(metabolon_clinic$bili, metabolon_clinic$compid_43807, method = "pearson", use = "complete.obs")
correlations["bilirubin","spearman"] <- cor(metabolon_clinic$bili, metabolon_clinic$compid_43807, method = "spearman", use = "complete.obs")

cor_res <- cor.test(metabolon_clinic$urea, metabolon_clinic$compid_1670, method = "pearson", use = "complete.obs")
plot_list[["urea"]] <- ggplot(metabolon_clinic, aes(x=urea, y=compid_1670)) +
  geom_point(aes(x=urea, y=compid_1670)) +
  geom_smooth(method = 'lm', formula = y~x) +
  xlab("clinical (mmol/L)") + ylab("Metabolon (scaled)") + ggtitle("C: Urea") +
  theme(panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_line(color = "lightgray", linewidth = 0.2))+
  annotate("text", x = 5, y = 3.5, label = paste0("r = ", round(correlations["urea","pearson"],2), " (",round(cor_res$conf.int[1],2),", ",round(cor_res$conf.int[2],2),")"), size = 3.5)

blandr_out <- blandr.statistics(scale(metabolon_clinic$urea), scale(metabolon_clinic$compid_1670), sig.level = 0.95, LoA.mode = 1)
blandr_out
plot_list_ba[["urea"]] <- blandr.plot.ggplot(blandr_out, method1name = "clinical", method2name = "MS", 
                                                plotTitle= "C: Urea")

cor_res <- cor.test(metabolon_clinic$creat, metabolon_clinic$compid_513, method = "pearson", use = "complete.obs")
plot_list[["creatinine"]] <- ggplot(metabolon_clinic, aes(x=creat, y=compid_513)) +
  geom_point(aes(x=creat, y=compid_513)) +
  geom_smooth(method = 'lm', formula = y~x) +
  xlab("clinical (µmol/L)") + ylab("Metabolon (scaled)") + ggtitle("D: Creatinine") +
  theme(panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_line(color = "lightgray", linewidth = 0.2))+
  annotate("text", x = 50, y = 2, label = paste0("r = ", round(correlations["creatinine","pearson"],2), " (",round(cor_res$conf.int[1],2),", ",round(cor_res$conf.int[2],2),")"), size = 3.5)

blandr_out <- blandr.statistics(scale(metabolon_clinic$urea), scale(metabolon_clinic$compid_513), sig.level = 0.95, LoA.mode = 1)
blandr_out
plot_list_ba[["creatinine"]] <- blandr.plot.ggplot(blandr_out, method1name = "clinical", method2name = "MS", 
                                             plotTitle= "D: Creatinine")

cor_res <- cor.test(metabolon_clinic$bili, metabolon_clinic$compid_43807, method = "pearson", use = "complete.obs")
plot_list[["bilirubin"]] <- ggplot(metabolon_clinic, aes(x=bili, y=compid_43807)) +
  geom_point(aes(x=bili, y=compid_43807)) +
  geom_smooth(method = 'lm', formula = y~x) +
  xlab("clinical (µmol/L)") + ylab("Metabolon (scaled)") + ggtitle("E: Bilirubin") +
  theme(panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_line(color = "lightgray", linewidth = 0.2))+
  annotate("text", x = 10, y = 4, label = paste0("r = ", round(correlations["bilirubin","pearson"],2), " (",round(cor_res$conf.int[1],2),", ",round(cor_res$conf.int[2],2),")"), size = 3.5)

blandr_out <- blandr.statistics(scale(metabolon_clinic$bili), scale(metabolon_clinic$compid_43807), sig.level = 0.95, LoA.mode = 1)
blandr_out
plot_list_ba[["bilirubin"]] <- blandr.plot.ggplot(blandr_out, method1name = "clinical", method2name = "MS", 
                                                   plotTitle= "E: Bilirubin")

# save out plots and table
filename = paste0(fig_dir,"FigS1cde_assay_comparison.pdf")
outputfig <- ggpubr::ggarrange(plotlist=plot_list,
                               ncol=2, nrow=2,
                               common.legend = T,
                               legend="top",
                               heights = 2,
                               widths = 2)
ggpubr::ggexport(outputfig, filename=filename)

filename = paste0(fig_dir,"Fig3cde_assay_comparison.jpeg")
outputfig <- ggpubr::ggarrange(plotlist=plot_list_ba,
                               ncol=2, nrow=2,
                               common.legend = T,
                               legend="top",
                               heights = 2,
                               widths = 2)
ggpubr::ggexport(outputfig, filename=filename,dpi=600)

## correlations
correlations <- as.data.frame(matrix(nrow=3, ncol = 2), row.names = names(plot_list))
colnames(correlations) <- c("pearson",  "spearman")
correlations["urea","pearson"] <- cor(metabolon_clinic$urea, metabolon_clinic$compid_1670, method = "pearson", use = "complete.obs")
correlations["urea","spearman"] <- cor(metabolon_clinic$urea, metabolon_clinic$compid_1670, method = "spearman", use = "complete.obs")
correlations["creatinine","pearson"] <- cor(metabolon_clinic$creat, metabolon_clinic$compid_513, method = "pearson", use = "complete.obs")
correlations["creatinine","spearman"] <- cor(metabolon_clinic$creat, metabolon_clinic$compid_513, method = "spearman", use = "complete.obs")
correlations["bilirubin","pearson"] <- cor(metabolon_clinic$bili, metabolon_clinic$compid_43807, method = "pearson", use = "complete.obs")
correlations["bilirubin","spearman"] <- cor(metabolon_clinic$bili, metabolon_clinic$compid_43807, method = "spearman", use = "complete.obs")

write.csv(correlations, file = paste0(table_dir, "clinical_metabolon_correlations.csv"), row.names = T)

# get N numbers for each comparison
na <- colSums(is.na(metabolon_clinic[,c("urea","compid_1670","creat","compid_513","bili","compid_43807")]))
na <- as.data.frame(na)

urea_n <- nrow(metabolon_clinic)-na["urea",]-na["compid_1670",]
urea_n
creat_n <- nrow(metabolon_clinic)-na["creat",]-na["compid_513",]
creat_n
bili_n <- nrow(metabolon_clinic)-na["bili",]-na["compid_43807",]
bili_n

# capture session info
print("Session information:")
sessionInfo()

# save output
sink()

# remove data
rm(list = ls())
