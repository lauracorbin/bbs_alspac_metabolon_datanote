*** TASK: Create ALSPAC at 30 Metabolon data in format for ALSPAC release

clear all
set more off
version 18.0

cd "/REDACETD/working/data/metabolon/step4/alspac/data_return"

****************************************************************************************************************

** set up a log of activity
capture log close
log using "/REDACETD/working/results/metabolon/logs/13.prepare.alspacf30.release.log", replace

****************************************************************************************************************

** READ IN DATA AND CONVERT TO STATA FORMAT

* metabolite data (modified header to read into stata)
** Create a macro with columns 2 to 1254
local numlist
forvalues i = 2/1254 {
    local numlist `numlist' `i'
}

// Use the macro in the import command
import delimited "temp_metab_data2.txt", varnames(1) numericcols(`numlist') clear

save "temp_metab_data.dta", replace

* sample metadata & summary data combined
import delimited "../metaboprep_release_2025_04_04/sumstats/filtered_dataset/alspac_only_2025_04_04_sample_anno_sumstats.txt", varnames(1) clear
drop client_identifier source_study

save "temp_sample_summarydata.dta", replace

* feature metadata and summary data combined
import delimited "../metaboprep_release_2025_04_04/sumstats/filtered_dataset/alspac_only_2025_04_04_feature_anno_sumstats.txt", varnames(1) clear
* update labels based on Metabolon key
label variable feature_names "An internally-generated Metabolon compound identifier (concatenated with compid_)"
label variable chem_id "Unique identifier for each biochemical : Use to merge with metabolite data"
label variable super_pathway "The general biochemical class assigned to the identified biochemical"
label variable sub_pathway "A more specific biochemical class assigned to the identified biochemical"
label variable pubchem "Identifier assigned by the NCBI and searchable in the PubChem database"
label variable platform "Method used to generate data"
label variable kegg "An alphanumeric compound identifier and link to compound information in KEGG"
label variable hmdb "An alphanumeric compound identifier and link to compound information HMDB"
label variable pathway_sortorder "Internally-generated number intended for sorting compounds by SUPER PATHWAY and SUB PATHWAY"
label variable type "Indicates whether a compound is named or unnamed/unknown."
label variable inchikey "IUPAC textual chemical identifier derived from INChI"
label variable smiles "Simplified molecular-input line-entry system (SMILES) line notation string"
label variable chemical_name "The name of the identified biochemical"
label variable plot_name "A shorter name of the identified biochemical, convenient for plotting"
label variable cas "A unique numerical identifier assigned by the Chemical Abstracts Service (CAS)"
label variable chemspider "A numerical identifier as maintained in the ChemSpider database"

label variable feature_missingness "Proportion of missing (NA) data"
label variable outlier_count "Number of outlying values (+/- 50 IQR from median)"
label variable n "Sample size (non-NA)" 
label variable mean "Mean"
label variable sd "Standard deviation"
label variable median "Median"
label variable min "Minimimum value"
label variable max "Maximum value" 
label variable range "Range (max - min)"
label variable skew "Skewness"
label variable kurtosis "Kurtosis"
label variable se "Syandard error"
label variable missing "Count of NA"
label variable var "Variance"
label variable disp_index "Variance to mean ratio"
label variable coef_variance "Coefficient of variation"
label variable w "shaprio's W-statistic of normality on raw distribution"
label variable log10_w "shaprio's W-statistic of normality on log10 distribution"
label variable k "Cluster assignment"
label variable independent_features_binary "Feature is a representative independent feature (1=yes)"

save "B4132_metabolite_metadata_20250604.dta", replace

** full (pre-QC sample list)

import delimited "../../../step2/alspac/excel/2025_04_04_alspac_metabolon_labids_B4132.csv", delim(",") clear
rename v1 client_sample_id
save "complete_sample_list.dta", replace

****************************************************************************************************************

* MERGE METABOLITE DATA WITH SAMPLE METADATA

use "temp_sample_summarydata.dta", clear

merge 1:1 parent_sample_name using "temp_metab_data.dta", nogen

* lab
label variable neg "Batch identifiers for the neg platform."
label variable polar "Batch identifiers for the polar platform."
label variable posearly "Batch identifiers for the posearly platform."
label variable poslate "Batch identifiers for the poslate platform."
label variable box_number "Shipment box number"
label variable client_matrix "Sample type"
label variable client_sample_id "Individual ID"
label variable client_sample_number "Run order"
label variable group_name "Assigned plate number (36 per plate)"
label variable sample_amount "Volume of material sent to Metabolon"
label variable sample_amount_units "Units of volume"
label variable sample_box_location "Shipment box location"
label variable vol_extracted "Volume of material extracted by Metabolon for analysis"
label variable within_set_sample_number "Within plate order"
label variable sample_missingness "Proportion of missing (NA) data"
label variable sample_missingness_w_exclusions "Proportion of missing (NA) data after removing xenobiotics"
label variable tpa_total "total peak area"
label variable tpa_completefeature "total peak area based on complete features only"
label variable outlier_count "Number of outlying values (+/- 50 IQR from median)"

* drop the metabolon assigned ID (not needed now sample metadata merged)
drop parent_sample_name
* bring lab ID to the first column - this is what needs to be replaced with aln/qlet during linkage
order client_sample_id

* merge with sample list
merge 1:1 client_sample_id using "complete_sample_list.dta", nogen


save "B4132_metabolite_data_20250604.dta", replace

****************************************************************************************************************
* delete temp files
rm "temp_metab_data.dta"
rm "temp_metab_data2.txt"
rm "temp_sample_summarydata.dta"
rm "complete_sample_list.dta"
****************************************************************************************************************

*---------------------------------------------------------------
clear all
*---------------------------------------------------------------

log close
