# bbs\_plus\_alspac_datanote

This repository contains code used to process the By-Band-Sleeve and ALSPAC Metabolon dataset and accompanies the datanote that describes the dataset creation (Untargeted metabolomics data in the By-Band-Sleeve trial and ALSPAC: integrating clinical trial and population cohort data').

## Code list

1) _01.read.in.metabolon.new.R_ - contains and runs a function to read in new format Metabolon Excel files. Adapted from original Metabolon read in code contained in R package _metaboprep_ by Maddy Smith. Also, reinstates missing data points in volume normalised version of data.  

2) _02.run.metaboprep.step1.sh_ - contains code to run the _metaboprep_ pipeline on the original dataset (n=2128 samples, m=1331 features), both the peak area data (_parameter\_file\_peak\_area\_full\_dataset.txt_) and the volume normalised data with missing data reinstated (_parameter\_file\_vol\_norm\_full\_dataset.txt_). Calls the _run\_metaboprep\_pipeline.R_ script.

3) _03.process.metabolon.R_ - contains code to separate out the three datasets (BBS only, ALSPAC only, BBS & ALSPAC) and remove the outlying cluster of samples from a single site identified in (2).   

4) _04.read.in.metabolon.new.R_ - contains and runs a function to read in new format Metabolon Excel files (as (1)).   

5) _05.run.metaboprep.step3.sh_ - contains code to run the _metaboprep_ pipeline on the combined volume normalised dataset (with missing data reinstated) with outlying cluster removed (n=1986 samples, m=1331 features) (_parameter\_file\_vol\_norm\_full\_dataset\_no\_outliers.txt_). Calls the _run\_metaboprep\_pipeline.R_ script.  

6) _06.remove.woc.from.joint.data.R_ - contains code to remove data belonging to BBS participants who have removed consent from the combined dataset.  

7) _07.split.data.by.study.R_ - contains code to read in output from (5) and (6) and generate BBS only and ALSPAC only datasets.  

8) _08.build.datasets.R_ - contains code to run _metaboprep_ pipeline on each of the three datasets: BBS, ALSPAC, BBS & ALSPAC. Also produces datasets and associated dictionaries for sharing in OSF. Output: html report from _metaboprep_ (v2) and data dictionary for each of the three datasets.

9) _09.DataSetSplittingNormalization.Rmd_ - script written by David Hughes to explore the impact of normalizing the whole dataset versus normalizing within study (BBS main, BBS sub-study, ALSPAC). 

10) _10.Explore.PCA.Structure.Rmd_ - code written to explore source of substructure in original dataset. 

11) _11.extract.alspac.pheno.do_ - code to extract alspac variables from age 30 clinic (when samples were collected). 

12) _12.modify.metabolite.data.format.sh_ - code to change format of metabolite data file output from metaboprep so that it will read into stata with the import delim function.  

13) _13.prepare.alspacf30.release.do_ - code to produce a proposed data return for alspac to be incorporated into the main catalogue.  

14) _14.alspac.pheno.summary.R_ - code to run phenotypic summaries for Table in datanote.

15) _15.clin.vs.metabolon.R_ - contains code to compare Metabolon data to clinical assays for glucose and cholesterol.  

16) _16.clin.vs.metabolon.bbsonly.R_ - contains code to compare Metabolon data to clinical assays for creatine, bilirubin & urea. 

17) _17.compare.siteB.Rmd_ - contains code to compare pilot data to main analysis data.  
  
clinical\_data\_prep:

1) 01.clinical.data.prep.bbsmain.R - takes csv extracts that contain main trial data and were supplied by Bristol Trials Centre and reformats into flat text as well as deriving additional variables.

2) 02.clinical.data.prep.bbssubstudy.R - takes csv extracts that contain substudy data and were supplied by Bristol Trials Centre and reformats into flat text as well as deriving additional variables.

3) merge.bbs.clinical.data.R - merges clinical data with metabolon manifest to restrict to those with metabolite data

4) 03a.merge.bbs.nmr.clinical.data.R - merges clinical data with nmr manifest to restrict to those with metabolite data

5) 04.explore.bbs.clinical.data.Rmd - summarises clinical data for datanote.

parameter\_file\_* - parameter files for running _metaboprep_.


