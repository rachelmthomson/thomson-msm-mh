clear all    // drops any data currently in the working directory, 
macro drop _all // drops an macros currently in the memory

global filepath "" // enter filepath for data here
display "$filepath"

log using "All.smcl", replace // enter desired log name here

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

program define PAFanalyse, rclass // program for main analysis

display "Estimate"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 1/20 {
	qui logit poverty2 $tinvarlist $tvarlist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.poverty2/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.poverty2/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if poverty2==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if poverty2==0 // Sampling weights for the non - treated group

forvalues i = 1/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

logistic ghqcase i.poverty2 $tinvarlist $tvarlist if incCC == 1 & _mi_m > 0 [ pw = ipw ], cluster(pidp) 

display "PAF calculation"

margins poverty2

local unexpprev = r(table)[1,1]
display "Prevalence in unexposed:  " `unexpprev'

prop ghqcase if incCC == 1 & _mi_m > 0 [pw = ipw]
local totprev = r(table)[1,2]
display `totprev'

local PAR = `totprev' - `unexpprev'
display `PAR'

return scalar ratio = `PAR'/`totprev'
display `PAR'/`totprev'

drop ipw ipw0 ipw1

end

program define PAFanalysesens, rclass // program for analysis excluding lagged poverty variable for use with transition exposure variable

display "Estimate"

qui gen double ipw0 = .
qui gen double ipw1 = .

forvalues i = 1/20 {
	qui logit poverty2 $tinvarlist $tvarsenslist if incCC == 1 & _mi_m == `i', nolog cluster(pidp) // propensity score model for exposure
	predict double ps if _mi_m == `i' & incCC == 1 // propensity score prediction
	replace ipw0 = 0.poverty2/(1-ps) if _mi_m == `i' & incCC == 1 
	replace ipw1 = 1.poverty2/ps if _mi_m == `i' & incCC == 1 
	drop ps
}

qui gen double ipw=. 
qui replace ipw = ipw1 if poverty2==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if poverty2==0 // Sampling weights for the non - treated group

forvalues i = 1/20 {
	qui sum ipw if _mi_m == `i'
	qui replace ipw = ipw/r(mean) if _mi_m == `i'
} // normalising weights

logistic ghqcase i.poverty2 $tinvarlist $tvarsenslist if incCC == 1 & _mi_m > 0 [ pw = ipw ], cluster(pidp) 

display "PAF calculation"

margins poverty2

local unexpprev = r(table)[1,1]
display "Prevalence in unexposed:  " `unexpprev'

prop ghqcase if incCC == 1 & _mi_m > 0 [pw = ipw]
local totprev = r(table)[1,2]
display `totprev'

local PAR = `totprev' - `unexpprev'
display `PAR'

return scalar ratio = `PAR'/`totprev'
display `PAR'/`totprev'

drop ipw ipw0 ipw1

end

program define PAFanalysecc, rclass // program for complete case analysis

display "Estimate"

qui gen double ipw0 = .
qui gen double ipw1 = .

qui logit poverty2 $tinvarlist $tvarlist if incCC == 1, nolog cluster(pidp) // propensity score model for exposure
predict double ps if incCC == 1 // propensity score prediction
replace ipw0 = 0.poverty2/(1-ps) if incCC == 1 
replace ipw1 = 1.poverty2/ps if incCC == 1 
drop ps
	
qui gen double ipw=. 
qui replace ipw = ipw1 if poverty2==1 // Sampling weights for the treated group
qui replace ipw = ipw0 if poverty2==0 // Sampling weights for the non - treated group

qui sum ipw
qui replace ipw = ipw/r(mean) // normalising weights

logistic ghqcase i.poverty2 $tinvarlist $tvarlist if incCC == 1 [ pw = ipw ], cluster(pidp) 

display "PAF calculation"

margins poverty2

local unexpprev = r(table)[1,1]
display "Prevalence in unexposed:  " `unexpprev'

prop ghqcase if incCC == 1 [pw = ipw]
local totprev = r(table)[1,2]
display `totprev'

local PAR = `totprev' - `unexpprev'
display `PAR'

return scalar ratio = `PAR'/`totprev'
display `PAR'/`totprev'

drop ipw ipw0 ipw1

end

*****EVERYONE*****
******************

bootstrap ratio = r(ratio), size(45497) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: PAFanalyse

log close

*************************************************
*********************MEN ONLY********************
*************************************************

log using "${filepath}\Men.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if sex != 1

bootstrap ratio = r(ratio), size(20891) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: PAFanalyse

log close

*****************************************************************
*************************WOMEN ONLY******************************
*****************************************************************

log using "${filepath}\Women.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if sex != 0

bootstrap ratio = r(ratio), size(24634) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: PAFanalyse

log close

*****************************************************************
************************HIGH EDUCATION***************************
*****************************************************************

log using "${filepath}\High ed.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 1

bootstrap ratio = r(ratio), size(19987) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: PAFanalyse

log close

*****************************************************************
***********************MEDIUM EDUCATION**************************
*****************************************************************

log using "${filepath}\Med ed.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 2 

bootstrap ratio = r(ratio), size(19576) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: PAFanalyse

log close

*****************************************************************
*************************LOW EDUCATION***************************
*****************************************************************

log using "${filepath}\Low ed.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if education != 3

bootstrap ratio = r(ratio), size(12636) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: PAFanalyse

log close

*****************************************************************
******************YOUNGER WORKING AGE ONLY***********************
*****************************************************************

log using "${filepath}\Young.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if agecat != 0

bootstrap ratio = r(ratio), size(20384) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: PAFanalyse

log close

*****************************************************************
*******************OLDER WORKING AGE ONLY************************
*****************************************************************

log using "${filepath}\Old.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if agecat != 1

bootstrap ratio = r(ratio), size(29711) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: PAFanalyse

log close

*****************************************************************
***********SENS. ANALYSIS - INTO POVERTY ONLY***************
*****************************************************************

log using "${filepath}\Into pov.smcl", replace

restore
preserve

*reconfiguring employment variable so that usual files can be used

drop poverty2
rename intopov poverty2
replace incCC = 0 if missing(poverty2)

bootstrap ratio = r(ratio), size(39772) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: PAFanalysesens

log close

*****************************************************************
***********SENS. ANALYSIS - OUT OF POVERTY ONLY***********************
*****************************************************************

log using "${filepath}\Out of pov.smcl", replace

restore
preserve

*reconfiguring employment variable so that usual files can be used

drop poverty2
rename povstatus poverty2
replace incCC = 0 if missing(poverty2)

bootstrap ratio = r(ratio), size(18206) reps(1000) cluster(pidp _mi_m) seed(6837) noisily: PAFanalysesens

log close

*****************************************************************
********************COMPLETE CASE ANALYSIS***********************
*****************************************************************

log using "${filepath}\CCs.smcl", replace

restore
preserve

*dropping those we don't want to include before generating weights

replace incCC = 0 if _mi_m > 0

bootstrap ratio = r(ratio), reps(1000) cluster(pidp) seed(6837) noisily: PAFanalysecc

log close
