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

global tvardirectlist "dvage age2 i.employlag1 nchild_dvlag1 i.partneredlag1 i.benstatuslag1 i.home_ownerlag1 i.gor_dvlag1 i.ghqcaselag1 sf12pcs_dvlag1 sf12mcs_dvlag1 i.poverty2lag1 i.poverty2 equivhhincminuscosts equivhhincminuscostslag1" // time-varying confounders for direct effect

global tvartotallist "dvage age2 i.employlag1 nchild_dvlag1 i.partneredlag1 i.benstatuslag1 i.home_ownerlag1 i.gor_dvlag1 i.ghqcaselag1 sf12pcs_dvlag1 sf12mcs_dvlag1 i.poverty2lag1 equivhhincminuscostslag1" // time-varying confounders for total effect

global tvardirectsenslist "dvage age2 nchild_dvlag1 i.partneredlag1 i.benstatuslag1 i.home_ownerlag1 i.gor_dvlag1 i.ghqcaselag1 sf12pcs_dvlag1 sf12mcs_dvlag1 i.poverty2lag1 i.poverty2 equivhhincminuscosts equivhhincminuscostslag1" // time-varying confounders for direct effect excluding lagged employment

global tvartotalsenslist "dvage age2 nchild_dvlag1 i.partneredlag1 i.benstatuslag1 i.home_ownerlag1 i.gor_dvlag1 i.ghqcaselag1 sf12pcs_dvlag1 sf12mcs_dvlag1 i.poverty2lag1 equivhhincminuscostslag1" // time-varying confounders for total effect excluding lagged employment

*Creating programs for each analysis

program define analysedirect // program for direct effect

display "*CAUSAL ESTIMATES*"
display "******************"

display "Generating weights"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit employ i.wavenum $tinvarlist $tvardirectlist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.employ/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.employ/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

display "Odds ratio"
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.employ $tinvarlist $tvardirectlist if incCC == 1 [ pw = ipw ], cluster(pidp) 

display "Average treatment effect"
mimrgns, dydx(employ) predict(default) esampvaryok 

display "Prevalence"
mimrgns employ, predict(default) esampvaryok

display "PAF"
local unexpprev = r(table)[1,1]
display "Prevalence in unexposed:  " `unexpprev'

prop ghqcase if incCC == 1 & _mi_m > 0 [pw = ipw]
local totprev = r(table)[1,2]
display `totprev'

local PAR = `totprev' - `unexpprev'
display `PAR'

display "PAF:"
display `PAR'/`totprev'

display "Confounder balance"

teffects ipw (ghqcase) (employ i.wavenum $tinvarlist $tvardirectlist, logit) if incCC == 1 & _mi_m > 0, vce (cluster pidp)
tebalance summarize // shows confounder balance

end

program define analysetotal // program for total effect

display "*CAUSAL ESTIMATES*"
display "******************"

display "Generating weights"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit employ i.wavenum $tinvarlist $tvartotallist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.employ/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.employ/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

display "Odds ratio"
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.employ $tinvarlist $tvartotallist if incCC == 1 [ pw = ipw ], cluster(pidp) 

display "Average treatment effect"
mimrgns, dydx(employ) predict(default) esampvaryok 

display "Prevalence"
mimrgns employ, predict(default) esampvaryok

display "PAF"
local unexpprev = r(table)[1,1]
display "Prevalence in unexposed:  " `unexpprev'

prop ghqcase if incCC == 1 & _mi_m > 0 [pw = ipw]
local totprev = r(table)[1,2]
display `totprev'

local PAR = `totprev' - `unexpprev'
display `PAR'

display "PAF:"
display `PAR'/`totprev'

display "Confounder balance"

teffects ipw (ghqcase) (employ i.wavenum $tinvarlist $tvartotallist, logit) if incCC == 1 & _mi_m > 0, vce (cluster pidp)
tebalance summarize // shows confounder balance

end

program define analysedirectsens // program for direct effect excluding lagged employment for use with transition exposure variable

display "*CAUSAL ESTIMATES*"
display "******************"

display "Generating weights"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit employ i.wavenum $tinvarlist $tvardirectsenslist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.employ/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.employ/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

display "Odds ratio"
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.employ $tinvarlist $tvardirectsenslist if incCC == 1 [ pw = ipw ], cluster(pidp) 

display "Average treatment effect"
mimrgns, dydx(employ) predict(default) esampvaryok 

display "Prevalence"
mimrgns employ, predict(default) esampvaryok

display "PAF"
local unexpprev = r(table)[1,1]
display "Prevalence in unexposed:  " `unexpprev'

prop ghqcase if incCC == 1 & _mi_m > 0 [pw = ipw]
local totprev = r(table)[1,2]
display `totprev'

local PAR = `totprev' - `unexpprev'
display `PAR'

display "PAF:"
display `PAR'/`totprev'

display "Confounder balance"

teffects ipw (ghqcase) (employ i.wavenum $tinvarlist $tvardirectsenslist, logit) if incCC == 1 & _mi_m > 0, vce (cluster pidp)
tebalance summarize // shows confounder balance

end

program define analysetotalsens // program for total effect excluding lagged employment for use with transition exposure variable

display "*CAUSAL ESTIMATES*"
display "******************"

display "Generating weights"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit employ i.wavenum $tinvarlist $tvartotalsenslist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.employ/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.employ/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

display "Odds ratio"
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.employ $tinvarlist $tvartotalsenslist if incCC == 1 [ pw = ipw ], cluster(pidp) 

display "Average treatment effect"
mimrgns, dydx(employ) predict(default) esampvaryok 

display "Prevalence"
mimrgns employ, predict(default) esampvaryok

display "PAF"
local unexpprev = r(table)[1,1]
display "Prevalence in unexposed:  " `unexpprev'

prop ghqcase if incCC == 1 & _mi_m > 0 [pw = ipw]
local totprev = r(table)[1,2]
display `totprev'

local PAR = `totprev' - `unexpprev'
display `PAR'

display "PAF:"
display `PAR'/`totprev'

display "Confounder balance"

teffects ipw (ghqcase) (employ i.wavenum $tinvarlist $tvartotalsenslist, logit) if incCC == 1 & _mi_m > 0, vce (cluster pidp)
tebalance summarize // shows confounder balance

end

********************************************************************************
*****************DIRECT EFFECT I.E. CONTROLLING FOR CONCURRENT INCOME***********
********************************************************************************

analysedirect

**********************************************
************ONLY THOSE IN POVERTY*************
**********************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if poverty2 != 1

analysedirect

**************************************************
************ONLY THOSE NOT IN POVERTY*************
**************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if poverty2 != 0

analysedirect

*************************************************
*********************MEN ONLY********************
*************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if sex != 1

analysedirect

*****************************************************************
*************************WOMEN ONLY******************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if sex != 0

analysedirect

*****************************************************************
************************HIGH EDUCATION***************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 1

analysedirect

*****************************************************************
***********************MEDIUM EDUCATION**************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 2

analysedirect

*****************************************************************
*************************LOW EDUCATION***************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 3

analysedirect

*****************************************************************
******************YOUNGER WORKING AGE ONLY***********************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if agecat != 0

analysedirect

*****************************************************************
*******************OLDER WORKING AGE ONLY************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if agecat != 1

analysedirect

*****************************************************************
***********SENS. ANALYSIS - INTO UNEMPLOYMENT ONLY***************
*****************************************************************

restore
preserve

*reconfiguring employment variable so that usual files can be used

drop employ
rename newjobloss employ
replace incCC = 0 if missing(employ)

analysedirectsens

*****************************************************************
***********SENS. ANALYSIS - INTO WORK ONLY***********************
*****************************************************************

restore
preserve

*reconfiguring employment variable so that usual files can be used

drop employ
rename gotajob employ
replace incCC = 0 if missing(employ)

analysedirectsens

*****************************************************************
********************COMPLETE CASE ANALYSIS***********************
*****************************************************************
restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if mnum > 0

analysedirect

************************************************************************************************
*extra sensitivity analysis using survey weights (only run ons those with longitudinal weights)*
************************************************************************************************

restore
preserve

replace incCC = 0 if lweight == 0 
replace incCC = 0 if missing(lweight)

*Causal estimate, logistic - no use of survey weights

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit employ i.wavenum $tinvarlist $tvardirectlist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.employ/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.employ/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

*Odds ratio
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.employ $tinvarlist $tvardirectlist if incCC == 1 [ pw = ipw ], cluster(pidp) 

*Average treatment effect
mimrgns, dydx(employ) predict(default) esampvaryok 

*Causal estimate, logistic - survey weights multiplied in

replace ipw = ipw * lweight

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

*Odds ratio
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.employ $tinvarlist $tvardirectlist if incCC == 1 [ pw = ipw ], cluster(pidp) 

*Average treatment effect
mimrgns, dydx(employ) predict(default) esampvaryok 

********************************************************************************
*****TOTAL EFFECT I.E. NOT CONTROLLING FOR CONCURRENT INCOME********
********************************************************************************

restore
preserve

analysetotal

**********************************************
************ONLY THOSE IN POVERTY*************
**********************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if poverty2 != 1

analysetotal

**************************************************
************ONLY THOSE NOT IN POVERTY*************
**************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if poverty2 != 0

analysetotal

*************************************************
*********************MEN ONLY********************
*************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if sex != 1

analysetotal

*****************************************************************
*************************WOMEN ONLY******************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if sex != 0

analysetotal

*****************************************************************
************************HIGH EDUCATION***************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 1

analysetotal

*****************************************************************
***********************MEDIUM EDUCATION**************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 2 

analysetotal

*****************************************************************
*************************LOW EDUCATION***************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 3

analysetotal

*****************************************************************
******************YOUNGER WORKING AGE ONLY***********************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if agecat != 0

analysetotal

*****************************************************************
*******************OLDER WORKING AGE ONLY************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if agecat != 1

analysetotal

*****************************************************************
***********SENS. ANALYSIS - INTO UNEMPLOYMENT ONLY***************
*****************************************************************

restore
preserve

*reconfiguring employment variable so that usual files can be used

drop employ
rename newjobloss employ
replace incCC = 0 if missing(employ)

analysetotalsens

*****************************************************************
***********SENS. ANALYSIS - INTO WORK ONLY***********************
*****************************************************************

restore
preserve

*reconfiguring employment variable so that usual files can be used

drop employ
rename gotajob employ
replace incCC = 0 if missing(employ)

analysetotalsens

*****************************************************************
********************COMPLETE CASE ANALYSIS***********************
*****************************************************************
restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if mnum > 0

analysetotal

************************************************************************************************
*extra sensitivity analysis using survey weights (only run ons those with longitudinal weights)*
************************************************************************************************

restore
preserve

replace incCC = 0 if lweight == 0 
replace incCC = 0 if missing(lweight)

*Causal estimate, logistic - no use of survey weights

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit employ i.wavenum $tinvarlist $tvartotallist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.employ/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.employ/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if employ==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if employ==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

*Odds ratio
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.employ $tinvarlist $tvartotallist if incCC == 1 [ pw = ipw ], cluster(pidp) 

*Average treatment effect
mimrgns, dydx(employ) predict(default) esampvaryok 

*Causal estimate, logistic - survey weights multiplied in

replace ipw = ipw * lweight

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

*Odds ratio
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.employ $tinvarlist $tvartotallist if incCC == 1 [ pw = ipw ], cluster(pidp) 

*Average treatment effect
mimrgns, dydx(employ) predict(default) esampvaryok 

log close
