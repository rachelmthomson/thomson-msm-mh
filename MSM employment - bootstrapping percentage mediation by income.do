clear all    // drops any data currently in the working directory, 
macro drop _all // drops an macros currently in the memory

global filepath "" // enter filepath for data here
display "$filepath"

log using "", replace // enter desired log name here

*NOTES*
*******

/* 
Code for replication of IPTWs used in teffects command obtained from Statalist forum
https://www.statalist.org/forums/forum/general-stata-discussion/general/1419349-obtaining-weighted-means-and-sds-using-tebalance-summarize-after-teffects-ipw
Appears to be correct based on comparisons and generates plausible value for weighted total prevalence
*/

/* To calculate PAF:
1. First need PAR i.e. prevalence in whole population minus prevalence in unexposed
2. Can get unexposed prevalence from ipw output
3. Need to run simple proportion command to get prevalence in whole population
4. Once have PAR, divide by total prevalence to get PAF
*/

*Setting up and general descriptives*
*************************************

use "", clear // enter data file name/location
set more off 

mi xtset pidp wavenum
mi tsset, clear

preserve // saving dataset in this form, already mi xtset - should then be usable for stratified analyses when restored

*Creating global macros containing confounders for analysis (need to be global rather than local so can run within programs)

global tinvarlist "i.sex i.education i.bame" // time invariant confounders

global tvardirectlist "i.wavenum dvage age2 i.employlag1 nchild_dvlag1 i.partneredlag1 i.benstatuslag1 i.home_ownerlag1 i.gor_dvlag1 i.ghqcaselag1 sf12pcs_dvlag1 sf12mcs_dvlag1 i.poverty2lag1 i.poverty2 equivhhincminuscosts equivhhincminuscostslag1" // time-varying confounders for direct effect

global tvartotallist "i.wavenum dvage age2 i.employlag1 nchild_dvlag1 i.partneredlag1 i.benstatuslag1 i.home_ownerlag1 i.gor_dvlag1 i.ghqcaselag1 sf12pcs_dvlag1 sf12mcs_dvlag1 i.poverty2lag1 equivhhincminuscostslag1" // time-varying confounders for total effect

global tvardirectsenslist "i.wavenum dvage age2 nchild_dvlag1 i.partneredlag1 i.benstatuslag1 i.home_ownerlag1 i.gor_dvlag1 i.ghqcaselag1 sf12pcs_dvlag1 sf12mcs_dvlag1 i.poverty2lag1 i.poverty2 equivhhincminuscosts equivhhincminuscostslag1" // time-varying confounders for direct effect excluding lagged employment

global tvartotalsenslist "i.wavenum dvage age2 nchild_dvlag1 i.partneredlag1 i.benstatuslag1 i.home_ownerlag1 i.gor_dvlag1 i.ghqcaselag1 sf12pcs_dvlag1 sf12mcs_dvlag1 i.poverty2lag1 equivhhincminuscostslag1" // time-varying confounders for total effect excluding lagged employment

*Creating programs for each analysis

program define analyse, rclass // program for main analysis

display "Direct effect"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 1/20 {
	qui logit employ $tinvarlist $tvardirectlist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	qui predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	qui replace ipw0 = 0.employ/(1-ps) if _mi_m == `i' & incCC == 1 
	qui replace ipw1 = 1.employ/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

forvalues i = 1/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

logistic ghqcase i.employ $tinvarlist $tvardirectlist if incCC == 1 & _mi_m > 0 [ pw = ipw ], cluster(pidp) 
margins, dydx(employ)

local directeffect = r(table)[1,2]
display "Direct effect:  " `directeffect'

drop ipw ipw0 ipw1

display "Total effect"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 1/20 {
	qui logit employ $tinvarlist $tvartotallist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	qui predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	qui replace ipw0 = 0.employ/(1-ps) if _mi_m == `i' & incCC == 1 
	qui replace ipw1 = 1.employ/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

forvalues i = 1/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

logistic ghqcase i.employ $tinvarlist $tvartotallist if incCC == 1 & _mi_m > 0 [ pw = ipw ], cluster(pidp) 
margins, dydx(employ)

local totaleffect = r(table)[1,2]
display "Total effect:  " `totaleffect'

return scalar ratio = `directeffect' / `totaleffect'
display `directeffect' / `totaleffect'

drop ipw ipw0 ipw1

end

program define analysesens, rclass // program for analysis excluding lagged employment for use with transition exposure variable

display "Direct effect"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 1/20 {
	qui logit employ $tinvarlist $tvardirectsenslist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	qui predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	qui replace ipw0 = 0.employ/(1-ps) if _mi_m == `i' & incCC == 1 
	qui replace ipw1 = 1.employ/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

forvalues i = 1/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

logistic ghqcase i.employ $tinvarlist $tvardirectsenslist if incCC == 1 & _mi_m > 0 [ pw = ipw ], cluster(pidp) 
margins, dydx(employ)

local directsenseffect = r(table)[1,2]
display "Direct effect:  " `directsenseffect'

drop ipw ipw0 ipw1

display "Total effect"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 1/20 {
	qui logit employ $tinvarlist $tvartotalsenslist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	qui predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	qui replace ipw0 = 0.employ/(1-ps) if _mi_m == `i' & incCC == 1 
	qui replace ipw1 = 1.employ/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

forvalues i = 1/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

logistic ghqcase i.employ $tinvarlist $tvartotalsenslist if incCC == 1 & _mi_m > 0 [ pw = ipw ], cluster(pidp) 
margins, dydx(employ)

local totalsenseffect = r(table)[1,2]
display "Total effect:  " `totalsenseffect'

return scalar ratio = `directsenseffect' / `totalsenseffect'
display `directsenseffect' / `totalsenseffect'

drop ipw ipw0 ipw1

end

program define analysecc, rclass // program for analysis of complete cases only

display "Direct effect"

qui gen double ipw0 = .
qui gen double ipw1 = .

qui logit employ $tinvarlist $tvardirectlist if incCC == 1, nolog cluster(pidp) // propensity score model for exposure
qui predict double ps if incCC == 1 // propensity score prediction
qui replace ipw0 = 0.employ/(1-ps) if incCC == 1 
qui replace ipw1 = 1.employ/ps if incCC == 1 
drop ps

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

qui sum ipw
qui replace ipw = ipw/r(mean) // normalising weights

logistic ghqcase i.employ $tinvarlist $tvardirectlist if incCC == 1 [ pw = ipw ], cluster(pidp) 
margins, dydx(employ)

local directeffect = r(table)[1,2]
display "Direct effect:  " `directeffect'

drop ipw ipw0 ipw1

display "Total effect"

qui gen double ipw0 = .
qui gen double ipw1 = .

qui logit employ $tinvarlist $tvartotallist if incCC == 1, nolog cluster(pidp) // propensity score model for exposure
qui predict double ps if incCC == 1 // propensity score prediction
qui replace ipw0 = 0.employ/(1-ps) if incCC == 1 
qui replace ipw1 = 1.employ/ps if incCC == 1 
drop ps

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

qui sum ipw
qui replace ipw = ipw/r(mean) // normalising weights

logistic ghqcase i.employ $tinvarlist $tvartotallist if incCC == 1 [ pw = ipw ], cluster(pidp) 
margins, dydx(employ)

local totaleffect = r(table)[1,2]
display "Total effect:  " `totaleffect'

return scalar ratio = `directeffect' / `totaleffect'
display `directeffect' / `totaleffect'

drop ipw ipw0 ipw1

end

*****EVERYONE*****
******************

bootstrap ratio = r(ratio), size(45497) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analyse

log close

**********************************************
************ONLY THOSE IN POVERTY*************
**********************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping in pov - manual MSM.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if poverty2 != 1

bootstrap ratio = r(ratio), size(18135) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analyse

log close

**************************************************
************ONLY THOSE NOT IN POVERTY*************
**************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping not in pov - manual MSM.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if poverty2 != 0

bootstrap ratio = r(ratio), size(40023) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analyse

log close

*************************************************
*********************MEN ONLY********************
*************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping men - manual MSM.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if sex != 1

bootstrap ratio = r(ratio), size(20891) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analyse

log close

*****************************************************************
*************************WOMEN ONLY******************************
*****************************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping women - manual MSM.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if sex != 0

bootstrap ratio = r(ratio), size(24634) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analyse

log close

*****************************************************************
************************HIGH EDUCATION***************************
*****************************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping high ed - manual MSM.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 1

bootstrap ratio = r(ratio), size(19987) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analyse

log close

*****************************************************************
***********************MEDIUM EDUCATION**************************
*****************************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping med ed - manual MSM.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 2 

bootstrap ratio = r(ratio), size(19576) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analyse

log close

*****************************************************************
*************************LOW EDUCATION***************************
*****************************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping low ed - manual MSM.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 3

bootstrap ratio = r(ratio), size(12638) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analyse

log close

*****************************************************************
******************YOUNGER WORKING AGE ONLY***********************
*****************************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping young - manual MSM.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if agecat != 0

bootstrap ratio = r(ratio), size(20384) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analyse

log close

*****************************************************************
*******************OLDER WORKING AGE ONLY************************
*****************************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping old - manual MSM.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if agecat != 1

bootstrap ratio = r(ratio), size(29711) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analyse

log close

*****************************************************************
***********SENS. ANALYSIS - INTO UNEMPLOYMENT ONLY***************
*****************************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping job loss - manual MSM.smcl", replace

restore
preserve

*reconfiguring employment variable so that usual files can be used

drop employ
rename newjobloss employ
replace incCC = 0 if missing(employ)

bootstrap ratio = r(ratio), size(36610) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analysesens

log close

*****************************************************************
***********SENS. ANALYSIS - INTO WORK ONLY***********************
*****************************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping job gain - manual MSM.smcl", replace

restore
preserve

*reconfiguring employment variable so that usual files can be used

drop employ
rename gotajob employ
replace incCC = 0 if missing(employ)

bootstrap ratio = r(ratio), size(18353) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: analysesens

log close

*****************************************************************
********************COMPLETE CASE ANALYSIS***********************
*****************************************************************

log using "${filepath}\Log Files\MSM imputed employment bootstrapping CCs - manual MSM.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if _mi_m > 0

bootstrap ratio = r(ratio), reps(10) cluster(pidp) seed(6837) noisily: analysecc

log close
