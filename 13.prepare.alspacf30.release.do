*** TASK: Create ALSPAC at 30 Metabolon data in format for ALSPAC release

clear all
set more off
version 19.0
set varabbrev off

cd "/REDACTED/working/data/metabolon/step4/alspac/data_return"

****************************************************************************************************************

** set up a log of activity
capture log close
log using "/REDACTED/working/results/metabolon/logs/13.prepare.alspacf30.release.log", replace

****************************************************************************************************************

** READ IN DATA AND CONVERT TO STATA FORMAT


*** sample metadata & summary data combined

* read in and convert labels
import delimited "../metaboprep_output/2026_05_05_SampleMetadataDictionary.txt", delim("\t") clear
replace variable_name = lower(variable_name)
replace variable_name = "poslate" if variable_name == "pos.late"
replace variable_name = "posearly" if variable_name == "pos.early"


tempname fh

local N = c(N)

// create a new do-file
file open `fh' using my_labels.do , write replace

forvalues i = 1/`N' {
          file write `fh' "label variable `= variable_name[`i']' "
          file write `fh' `""`= variable_label[`i']'""' _newline
}

file close `fh'


import delimited "../metaboprep_output/2026_05_05_SampleMetadata.txt", varnames(1) delim("\t") clear
do my_labels.do
rm my_labels.do

save "temp_SampleMetadata.dta", replace


*** metabolite data

* read in and convert labels
import delimited "../metaboprep_output/2026_05_05_MetaboliteDataDictionary.txt", delim("\t") clear

tempname fh

local N = c(N)

// create a new do-file
file open `fh' using my_labels.do , write replace

forvalues i = 1/`N' {
          file write `fh' "label variable `= variable_name[`i']' "
          file write `fh' `""`= variable_label[`i']'""' _newline
}

file close `fh'

import delimited "../metaboprep_output/2026_05_05_MetaboliteData.txt", varnames(1) delim("\t") clear

do my_labels.do
rm my_labels.do

save "temp_MetaboliteData.dta", replace

merge 1:1 sample_id using "temp_SampleMetadata.dta", nogenerate

save "B4132_metabolite_data_20260505.dta", replace

****************************************************************************************************************

** full (pre-QC sample list)
import delimited "../../../step2/alspac/excel/2025_04_04_alspac_metabolon_labids_B4132.csv", delim(",") clear
rename v1 client_sample_id
save "complete_sample_list.dta", replace

merge 1:1 client_sample_id using "B4132_metabolite_data_20260505.dta", nogenerate

****************************************************************************************************************

save "B4132_metabolite_data_20260505.dta", replace

****************************************************************************************************************
* delete temp files
rm temp_MetaboliteData.dta
rm temp_SampleMetadata.dta
rm complete_sample_list.dta
****************************************************************************************************************

*---------------------------------------------------------------
clear all
*---------------------------------------------------------------

log close
