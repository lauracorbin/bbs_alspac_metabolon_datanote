*** Syntax template for direct users preparing datasets using child and parent based datasets.


* Updated 14th March 


****************************************************************************************************************************************************************************************************************************
* This template is based on that used by the data buddy team and they include a number of variables by default.
* To ensure the file works we suggest you keep those in and just add any relevant variables that you need for your project.
* To add data other than that included by default you will need to add the relvant files and pathnames in each of the match commands below.
* There is a separate command for mothers questionnaires, mothers clinics, partner, mothers providing data on the child and data provided by the child themselves.
* Each has different withdrawal of consent issues so they must be considered separately.
* You will need to replace 'YOUR PATHNAME' in each section with your working directory pathname.

*****************************************************************************************************************************************************************************************************************************.

*****************************************************************************************************************************************************************************************************************************.

* SETUP

clear all
set more off
version 18.0

** set up a log of activity
capture log close
log using "/REDACTED/working/results/metabolon/logs/11.extract.alspac.pheno.log", replace
****************************************************************************************************************



* G0 Mother (pregnancy) based files - include here all files related to the pregnancy and/or mother

* If no mother variables are required, KEEP this section and remove the instruction below to run it..

clear
set maxvar 32767	
use "/REDACTED/current/other/cohort profile/G0/mother/mz_6a.dta", clear
sort aln
gen in_mz=1
merge 1:1 aln using "/REDACTED/current/quest/G0/mother/a_3e.dta", nogen
merge 1:1 aln using "/REDACTED/current/quest/G0/mother/b_4f.dta", nogen
merge 1:1 aln using "/REDACTED/current/quest/G0/mother/c_8a.dta", nogen
keep aln mz010a preg_in_alsp preg_in_core preg_enrol_status mz005l mz005m mz013 mz014 mz028b ///
a006 a525 ///
b032 b650 b663 - b667 ///
c645a c755 c765 ///
bestgest


*keep only those pregnancies enrolled in ALSPAC
keep if preg_enrol_status < 3 | aln == REDACTED

* Dealing with withdrawal of consent: For this to work additional variables required have to be inserted before bestgest, so replace the *** line above with additional variables. 
* If none are required remember to delete the *** line.
* An additional do file is called in to set those withdrawing consent to missing so that this is always up to date whenever you run this do file. Note that mother based WoCs are set to .a

order aln mz010a, first
order bestgest, last

do "/REDACTED/Syntax/Withdrawal of consent/mother_WoC.do"

* Check withdrawal of consent frequencies=30 and baseline number is 15447
tab1 mz010a, mis

save "/REDACTED/working/data/metabolon/alspac_phenotype_data/mother.dta", replace



*****************************************************************************************************************************************************************************************************************************.
/* G0 PARTNER - ***UNBLOCK SECTION WHEN REQUIRED***
* G0 Partner files - include here all files related to the G0 partner/father


use "\\ads.bris.ac.uk\filestore\SSCM ALSPAC\Data\current\other\cohort profile\G0\partner\pz_1a.dta"
sort aln
/*merge in any additional datasets if required*/

keep aln partner_in_alspac partner_data partner_enrolled partner_in_core pz_mult pz_multid partner_changed partner_changed_when ///
/* add your variable list here; remove this line if you are not adding any extra*/
partner_age second_partner_age


*keep in just those enrolled partners
keep if partner_in_alspac >= 1


* Removing withdrawl of consent cases *** FOR LARGE DATASETS THIS CAN TAKE A FEW MINUTES
* An additional do file is called in to set those withdrawing consent to missing so that this is always up to date whenever you run this do file. Note that partner based WoCs are set to .c

order aln partner_in_alspac, first
order partner_age second_partner_age, last

do "\\ads.bris.ac.uk\filestore\SSCM ALSPAC\Data\Syntax\Withdrawal of consent\partner_WoC.do"

* Check there is a total of n=12113 G0 partners.
* Check withdrawal of consent frequencies partner=5  - currently none of these are in the formal cohort
tab1 partner_in_alspac, mis

save "YOUR PATHNAME\partner.dta", replace 
*/





*****************************************************************************************************************************************************************************************************************************.
* G1 Child BASED files - in this section the following file types need to be placed:
* Mother completed Qs about YP
* Obstetrics file OA

* ALWAYS KEEP THIS SECTION EVEN IF ONLY CHILD COMPLETED REQUESTED, although you will need to remove the *****

use "/REDACTED/current/other/cohort profile/G1/cp_3a.dta", clear
sort aln qlet
gen in_kz=1
/* remove this line and end of commenting on the next line to add in any other file names if relevant files are needed in this section
merge 1:1 aln qlet using "", nogen   */


keep aln qlet kz021 kz030 kz011b ///
in_core in_alsp in_phase2 in_phase3 in_phase4 tripquad


* Dealing with withdrawal of consent: For this to work additional variables required have to be inserted before in_core, so replace the ***** line with additional variables.
* If none are required remember to delete the ***** line.
* An additional do file is called in to set those withdrawing consent to missing so that this is always up to date whenever you run this do file. Note that child based WoCs are set to .b


order aln qlet kz021, first
order in_alsp tripquad, last

do "/REDACTED/Syntax/Withdrawal of consent/child_based_WoC.do"

* Check withdrawal of consent frequencies child based=32 (two mums of twins have withdrawn consent)
tab1 kz011b, mis

drop kz021
save "/REDACTED/working/data/metabolon/alspac_phenotype_data/childB.dta", replace

*****************************************************************************************************************************************************************************************************************************.
* G1 Child COMPLETED files - in this section the following file types need to be placed:
* YP completed Qs
* Puberty Qs
* Child clinic data
* Child biosamples data
* School Qs
* Obstetrics file OC
* G1 IMD for years goes in this section e.g. jan1999imd2010_crimeq5_YP
* Child longitudinal data

* NOTE: Always keep this section even if no child completed variables requested. Remove the ***** line.
* NOTE: Keep kz021 tripquad just to make the withdrawal of consent work - these are dropped for this file as the ones in the child BASED file are the important ones and should take priority

use "/REDACTED/current/other/cohort profile/G1/cp_3a.dta", clear
sort aln qlet
gen in_kz=1
merge 1:1 aln qlet using "/REDACTED/current/quest/G1/self/YPH_4a.dta", nogen
merge 1:1 aln qlet using "/REDACTED/current/clinic/G1/F30_G1_1a.dta", nogen

keep aln qlet kz021 kz030 ///
YPH2010 FLSA1020 FLAR0010 FLMS1053 FLMS1060 ///
tripquad

* Dealing with withdrawal of consent: For this to work additional variables required have to be inserted before tripquad, so replace the ***** line with additional variables.
* An additional do file is called in to set those withdrawing consent to missing so that this is always up to date whenever you run this do file.  Note that mother based WoCs are set to .b

order aln qlet kz021, first
order tripquad, last

do "/REDACTED/Syntax/Withdrawal of consent/child_completed_WoC.do"

* Check withdrawal of consent frequencies child completed=32
tab1 kz021, mis

drop tripquad
save "/REDACTED/working/data/metabolon/alspac_phenotype_data/childC.dta", replace

*****************************************************************************************************************************************************************************************************************************.
** Matching all data together and saving out the final file*.
* NOTE: any linkage data should be added here*.

use "/REDACTED/working/data/metabolon/alspac_phenotype_data/childB.dta", clear
merge 1:1 aln qlet using "/REDACTED/working/data/metabolon/alspac_phenotype_data/childC.dta", nogen
merge m:1 aln using "/REDACTED/working/data/metabolon/alspac_phenotype_data/mother.dta", nogen
* IF partner data is required please unstar the following line
/* merge m:1 aln using "YOUR PATHWAY\partner.dta", nogen */


* Remove non-alspac children.
drop if in_alsp!=1 & aln != 33893

* Remove trips and quads.
drop if tripquad==1

drop in_alsp tripquad
save "/REDACTED/working/data/metabolon/alspac_phenotype_data/alspac_pheno.dta", replace

*****************************************************************************************************************************************************************************************************************************.
* QC checks*
use "/REDACTED/working/data/metabolon/alspac_phenotype_data/alspac_pheno.dta", clear

* Check that there are 15645 records.
count

*****************************************************************************************************************************************************************************************************************************.
* rename vars
rename kz021 sex
rename YPH2010 ethnicity
rename FLSA1020 insulin_user
rename FLAR0010 age_yrs
rename FLMS1053 weight_kg
rename FLMS1060 bmi_kgm2

* preprocess so ready for reading into R
replace bmi_kgm2 = . if bmi_kgm2 < 0
replace weight_kg = . if weight_kg < 0
replace age_yrs = . if age_yrs <0
replace insulin_user = . if insulin_user <0
decode ethnicity, generate(ethnicity_cat)
tab ethnicity ethnicity_cat

* export as csv
save "/REDACTED/working/data/metabolon/alspac_phenotype_data/11.alspac_pheno.dta", replace
***************************************************************************************************************************************************************************************************************************.

* remove files not needed
rm "/REDACTED/working/data/metabolon/alspac_phenotype_data/childC.dta"
rm "/REDACTED/working/data/metabolon/alspac_phenotype_data/childB.dta"
rm "/REDACTED/working/data/metabolon/alspac_phenotype_data/mother.dta"

***************************************************************************************************************************************************************************************************************************.


*---------------------------------------------------------------
clear all
*---------------------------------------------------------------

log close


