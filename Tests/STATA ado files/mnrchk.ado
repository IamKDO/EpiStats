capture program drop mnrchk
program  mnrchk
syntax varlist [ , flag]
local listvar = "`varlist'"

quietly : {

capture drop R_*
foreach var of varlist `listvar'   {
gen R_`var'=0
replace R_`var'=1 if missing(`var') 
}

noisily : display %12s "" _continue
foreach varcol of varlist `listvar' {
      noisily : display "  " %6s "`varcol'" _continue
}
noisily : display " "

foreach varrow of varlist `listvar' {
   noisily : display %12s "NR`varrow'" _continue
   foreach varcol of varlist `listvar'{
   if "`varrow'" == "`varcol'" { 
    noisily : display "   " %5s "-" _continue
	
   } 
   else {
      tab R_`varrow' `varcol', chi2
	if "`flag'" == "flag" {
    noisily : display "   " %5s cond(r(p)<=0.05,"***"," ")   _continue
	}
	else {
    noisily : display "   " %5.3f r(p) _continue
	}
	}
	}
    noisily : display " "
}	

}
	  
end
