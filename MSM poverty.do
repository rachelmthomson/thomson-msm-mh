clear all    // drops any data currently in the working directory, 
macro drop _all // drops an macros currently in the memory

global filepath "" // enter filepath for data here
display "$filepath"

log using "", replace // enter desired log name here

*Setting up and general descriptives*
*************************************

use "", clear // enter data file name/location
set more off 

mi xtset pidp wavenum
mi tsset, clear

preserve // saving dataset in this form, already mi xtset - should then be usable for stratified analyses when restored

*Creating global macros containing confounders for analysis (need to be global rather than local so can run within programs)

global tinvarlist "i.sex i.education i.bame" // time invariant confounders

global tvarlist "i.wavenum dvage age2 i.employ i.employlag1 nchild_dvlag1 i.partneredlag1 i.benstatuslag1 i.home_ownerlag1 i.gor_dvlag1 i.ghqcaselag1 sf12pcs_dvlag1 sf12mcs_dvlag1 i.poverty2lag1" // time-varying confounders

global tvarsenslist "i.wavenum dvage age2 i.employ i.employlag1 nchild_dvlag1 i.partneredlag1 i.benstatuslag1 i.home_ownerlag1 i.gor_dvlag1 i.ghqcaselag1 sf12pcs_dvlag1 sf12mcs_dvlag1" // time-varying confounders excluding lagged poverty

*Creating programs for each analysis

program define analyse // program for main analysis

display "*TRADITIONAL REGRESSION ESTIMATES*"
display "**********************************"

display "LOGISTIC, UNWEIGHTED:"
mi estimate, or esampvaryok: logistic ghqcase i.poverty2 if incCC == 1, cluster(pidp)
mimrgns, dydx(poverty2) predict(default) esampvaryok

display "FIXED-EFFECTS, UNWEIGHTED:"
mi xtset pidp wavenum
mi estimate, or esampvaryok: xtlogit ghqcase i.poverty2 if incCC == 1, fe
mimrgns, dydx(poverty2) predict(default) esampvaryok

display "FIXED-EFFECTS, UNWEIGHTED but with time-varying confounders:"
mi xtset pidp wavenum
mi estimate, or esampvaryok: xtlogit ghqcase i.poverty2 $tvarlist if incCC == 1, fe
mimrgns, dydx(poverty2) predict(default) esampvaryok

display "*CAUSAL ESTIMATES*"
display "******************"

display "Generating weights"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit poverty2 $tinvarlist $tvarlist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.poverty2/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.poverty2/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if poverty2==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if poverty2==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

display "Odds ratio"
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.poverty2 $tinvarlist $tvarlist if incCC == 1 [ pw = ipw ], cluster(pidp) 

display "Average treatment effect"
mimrgns, dydx(poverty2) predict(default) esampvaryok 

display "Prevalence"
mimrgns poverty2, predict(default) esampvaryok

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

teffects ipw (ghqcase) (poverty2 $tinvarlist $tvarlist, logit) if incCC == 1 & _mi_m > 0, vce (cluster pidp)
tebalance summarize // shows confounder balance

end

program define analysesens // program for analysis excluding lagged poverty variable for use with transition exposure variable

display "*TRADITIONAL REGRESSION ESTIMATES*"
display "**********************************"

display "FIXED-EFFECTS, UNWEIGHTED but with time-varying confounders:"
mi xtset pidp wavenum
mi estimate, or esampvaryok: xtlogit ghqcase i.poverty2 $tvarsenslist if incCC == 1, fe
mimrgns, dydx(poverty2) predict(default) esampvaryok

display "*CAUSAL ESTIMATES*"
display "******************"

display "Generating weights"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit poverty2 $tinvarlist $tvarsenslist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.poverty2/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.poverty2/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if poverty2==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if poverty2==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

display "Odds ratio"
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.poverty2 $tinvarlist $tvarsenslist if incCC == 1 [ pw = ipw ], cluster(pidp) 

display "Average treatment effect"
mimrgns, dydx(poverty2) predict(default) esampvaryok 

display "Prevalence"
mimrgns poverty2, predict(default) esampvaryok

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

teffects ipw (ghqcase) (poverty2 $tinvarlist $tvarsenslist, logit) if incCC == 1 & _mi_m > 0, vce (cluster pidp)
tebalance summarize // shows confounder balance - can't think of faster way to do this

end

*************************************************************************
*********************Main analysis on full sample************************
*************************************************************************

analyse

*extra sensitivity analysis using survey weights (only run ons those with longitudinal weights)*
************************************************************************************************

drop ps
drop ipw

replace incCC = 0 if lweight == 0 
replace incCC = 0 if missing(lweight)

*logistic, unweighted
mi estimate, or esampvaryok: logistic ghqcase i.poverty2 if incCC == 1, cluster(pidp)
mimrgns, dydx(poverty2) predict(default) esampvaryok

*logistic, longitudinal weights
mi estimate, or esampvaryok: logistic ghqcase i.poverty2 if incCC == 1 [pweight = lweight], cluster(pidp)
mimrgns, dydx(poverty2) predict(default) esampvaryok

*fixed-effects, unweighted
mi xtset pidp wavenum
mi estimate, or esampvaryok: xtlogit ghqcase poverty2 if incCC == 1, fe
mimrgns, dydx(poverty2) predict(default) esampvaryok

*fixed-effects, longitudinal weights - to compare with FE with IPTWs

mi svyset pidp [pw=lweight]
mi estimate, or esampvaryok: svy: clogit ghqcase poverty2 if incCC == 1, group(pidp)
mimrgns, dydx(poverty2) predict(default) esampvaryok

*fixed-effects, unweighted plus time-varying confounders
mi xtset pidp wavenum
mi estimate, or esampvaryok: xtlogit ghqcase i.poverty2 $tvarlist if incCC == 1, fe
mimrgns, dydx(poverty2) predict(default) esampvaryok

*fixed-effects, longitudinal weights plus time-varying confounders
mi svyset pidp [pw=lweight]
mi estimate, or esampvaryok: svy: clogit ghqcase i.poverty2 $tvarlist if incCC == 1, group(pidp)
mimrgns, dydx(poverty2) predict(default) esampvaryok

*Causal estimate, logistic - no use of survey weights

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit poverty2 $tinvarlist $tvarlist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.poverty2/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.poverty2/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if poverty2==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if poverty2==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

*Odds ratio
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.poverty2 $tinvarlist $tvarlist if incCC == 1 [ pw = ipw ], cluster(pidp) 

*Average treatment effect
mimrgns, dydx(poverty2) predict(default) esampvaryok 

*Causal estimate, logistic - survey weights multiplied in

replace ipw = ipw * lweight

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

*Odds ratio
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.poverty2 $tinvarlist $tvarlist if incCC == 1 [ pw = ipw ], cluster(pidp) 

*Average treatment effect
mimrgns, dydx(poverty2) predict(default) esampvaryok 

*************************************************
*********************MEN ONLY********************
*************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if sex != 1

analyse

*****************************************************************
*************************WOMEN ONLY******************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if sex != 0

analyse

*****************************************************************
************************HIGH EDUCATION***************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 1

analyse

*****************************************************************
***********************MEDIUM EDUCATION**************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 2 // this causes issue with N. Ireland in imputation 16 
// need to run these manually

*TRADITIONAL REGRESSION ESTIMATES*
**********************************

*LOGISTIC, UNWEIGHTED:
mi estimate, or esampvaryok: logistic ghqcase i.poverty2 if incCC == 1, cluster(pidp)
mimrgns, dydx(poverty2) predict(default) esampvaryok

*FIXED-EFFECTS, UNWEIGHTED:
mi xtset pidp wavenum
mi estimate, or esampvaryok: xtlogit ghqcase i.poverty2 if incCC == 1, fe
mimrgns, dydx(poverty2) predict(default) esampvaryok

drop if gor_dvlag1==12 // need to do this in advance of running this FE model because there is nobody with outcome variation in this group *except* in imputed dataset 16, which means omitted terms vary and mi estimate can't generate summary effects

*"FIXED-EFFECTS, UNWEIGHTED but with time-varying confounders:"
mi xtset pidp wavenum
mi estimate, or esampvaryok: xtlogit ghqcase i.poverty2 $tvarlist if incCC == 1, fe
mimrgns, dydx(poverty2) predict(default) esampvaryok

restore // getting back to original sample for causal estimates
preserve

replace incCC = 0 if education != 2 

display "Generating weights"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit poverty2 $tinvarlist $tvarlist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.poverty2/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.poverty2/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if poverty2==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if poverty2==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

display "Odds ratio"
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.poverty2 $tinvarlist $tvarlist if incCC == 1 [ pw = ipw ], cluster(pidp) 

display "Average treatment effect"
mimrgns, dydx(poverty2) predict(default) esampvaryok 

display "Prevalence"
mimrgns poverty2, predict(default) esampvaryok

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

teffects ipw (ghqcase) (poverty2 $tinvarlist $tvarlist, logit) if incCC == 1 & _mi_m > 0, vce (cluster pidp)
tebalance summarize // shows confounder balance

*****************************************************************
*************************LOW EDUCATION***************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 3

analyse

*****************************************************************
******************YOUNGER WORKING AGE ONLY***********************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if agecat != 0

analyse

*****************************************************************
*******************OLDER WORKING AGE ONLY************************
*****************************************************************

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if agecat != 1

analyse

*****************************************************************
***********SENS. ANALYSIS - INTO POVERTY ONLY********************
*****************************************************************

restore
preserve

*reconfiguring poverty variable so that usual files can be used

drop poverty2
rename intopov poverty2
replace incCC = 0 if missing(poverty2)

analysesens

*****************************************************************
***********SENS. ANALYSIS - OUT OF POVERTY ONLY******************
*****************************************************************

restore
preserve

*reconfiguring poverty variable so that usual files can be used

drop poverty2
rename povstatus poverty2
replace incCC = 0 if missing(poverty2) // this causes issue with N. Ireland in imputation 2 for FE analysis with time-varying confounders
// need to run these manually

*TRADITIONAL REGRESSION ESTIMATES*
**********************************

*LOGISTIC, UNWEIGHTED:
mi estimate, or esampvaryok: logistic ghqcase i.poverty2 if incCC == 1, cluster(pidp)
mimrgns, dydx(poverty2) predict(default) esampvaryok

*FIXED-EFFECTS, UNWEIGHTED:
mi xtset pidp wavenum
mi estimate, or esampvaryok: xtlogit ghqcase i.poverty2 if incCC == 1, fe
mimrgns, dydx(poverty2) predict(default) esampvaryok

drop if gor_dvlag1==12 // need to do this in advance of running this FE model because there is nobody with outcome variation in this group *except* in imputed dataset 16, which means omitted terms vary and mi estimate can't generate summary effects

*"FIXED-EFFECTS, UNWEIGHTED but with time-varying confounders:"
mi xtset pidp wavenum
mi estimate, or esampvaryok: xtlogit ghqcase i.poverty2 $tvarsenslist if incCC == 1, fe
mimrgns, dydx(poverty2) predict(default) esampvaryok

restore // getting back to original sample for causal estimates
preserve

drop poverty2
rename povstatus poverty2
replace incCC = 0 if missing(poverty2)

*CAUSAL ESTIMATES*
******************

display "Generating weights"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 0/20 {
	qui logit poverty2 $tinvarlist $tvarsenslist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.poverty2/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.poverty2/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if poverty2==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if poverty2==0 // Sampling weights for the non - treated group

forvalues i = 0/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

display "Odds ratio"
mi estimate, esampvaryok cmdok or post: logistic ghqcase i.poverty2 $tinvarlist $tvarsenslist if incCC == 1 [ pw = ipw ], cluster(pidp) 

display "Average treatment effect"
mimrgns, dydx(poverty2) predict(default) esampvaryok 

display "Prevalence"
mimrgns poverty2, predict(default) esampvaryok

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

teffects ipw (ghqcase) (poverty2 $tinvarlist $tvarsenslist, logit) if incCC == 1 & _mi_m > 0, vce (cluster pidp)
tebalance summarize // shows confounder balance

*****************************************************************
********************COMPLETE CASE ANALYSIS***********************
*****************************************************************

restore

*dropping those we don't want to include before generating weights

replace incCC = 0 if mnum > 0
preserve

analyse

log close
