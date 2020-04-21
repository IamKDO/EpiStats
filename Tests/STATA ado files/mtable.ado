
capture program drop mtable
program mtable

syntax varlist (min=2 max=2) [if] [in],  by(varname) 

tempvar touse
mark `touse' `if' `in'

tokenize `varlist'

quietly : count if `touse'
if r(N) == 0 {
   exit
}

if "`touse'" != "" {
local touse = " & `touse'"
}

quietly : {
// Get Total  case and controls
count if `1' != . & `2' != . & `by' != .  `touse'
local _Total = r(N)
count if `1' == . | `2' == . | `by' == .  `touse'
local _Missing = r(N)
}

local rowline = "-"
local rowname = ""

local _varsize = 12
local ilength =length("`2'")
if `ilength' > 12 {
local _varsize = `ilength'
}

capture drop _cpc _expc 
capture matrix drop _cpcvalues
preserve
bys `by' : egen _cpc = sum(cond(`1'==0,1,0))
bys `by' : egen _expc = sum(cond(`1'==0,`2',0))
quietly : sum _cpc
local maxcontrol = r(max)

quietly tab `2' _expc if `1' == 1 `touse', matcell(_cpcvalues)
local iRow = r(r)
local iCol = r(c)

if `iRow' < 2 | `iCol' < 2 {
   display in smcl as error "Not able to do a 2 by 2 table"
   exit
}

local stemp : variable label `2' 
if "`stemp'" == "" {
local stemp = "`2'"
}

display in smcl as text "Matched analysis by `by' " "(Maximum controls/cases = " `maxcontrol' ")" 
display in smcl as text "Exposure : `stemp' "
display in smcl as text "Number of obs = " %5.0f `_Total' " , Missing = " %5.0f `_Missing'
display in smcl as text "Distribution of pairs :"

local _num = 0
local _den = 0


forvalues icpc = 1(1)`maxcontrol' { 
	local icpc1 = `icpc' + 1
	// On fait les tables 
	quietly : count if `1' == 1 & _cpc == `icpc' `touse'
	if r(N) == 0 {
	   continue
	}
	display
	display in smcl as text "`icpc' controls per case : " r(N) " pairs"
    // Il faudrait ici remplir une table complete avec zéro !
	// ==> Créer une matrice de zéro et remplacer les valeurs disponibles obtenues avec la table	
    matrix _mvalues = J(2,`icpc1',0)

	quietly tab `2' _expc if `1' == 1 & _cpc == `icpc' `touse', matcell(_cpcvalues) matcol(_colvalues)	
	local iRow = r(r)
	local iCol = r(c)	
	// matrix list _cpcvalues 
	// matrix list _colvalues
	
	if "`iCol'" != "" {
		forvalues imat = 1(1)`iCol' {
		   local icurcol = el(_colvalues,1,`imat') + 1
		   matrix _mvalues[1,`icurcol']=el(_cpcvalues,1,`imat')
		   if `iRow' > 1 {
			matrix _mvalues[2,`icurcol']=el(_cpcvalues,2,`imat')
		   }
		}
	}
	// if `iRow' < 2 | `iCol' < 2 {
	//   display in smcl as text "No data for some cells"
	//   continue
	// }
	
	// préparation du col header
	local _colheader ="Exposed Unexposed"
	local _rowheader ="Exposed Unexposed"
	if `maxcontrol' > 1 {
		local _colheader = ""
		local _rowheader = "+ -"
		// we do the loop reverse as matrice is displayed reverse
		forvalues iloop = `icpc1'(-1)1 {	 
		    // DP = nombre de control exposed
			local DP = `iloop' - 1
			local _s1 = substr("++++++++++++",1,`DP')
			// DP = nombre de controls non exposed = nb control - control exp 
			local DP = `icpc1' - `iloop' 
			local _s2 = substr("------------",1,`DP')
			if "`_s1'" == "" | "`_s2'" == "" { 
			local _colheader = "`_colheader'`_s1'`_s2' "
			}
			else {
			local _colheader ="`_colheader'`_s1'/`_s2' "
			}	
		}
	}
	_table _mvalues ,coltitle("Controls") rowtitle("Cases") colheader("`_colheader'") /// 
		 rowheader("`_rowheader'") reverse 
	forvalues iloop = 1(1)`icpc1' {	 
	// Calcul du numérateur
	    // DP = nombre de controls non exposed = nb control - control exp 
		local DP = `icpc1' - `iloop' 
	    local n1 = el(_mvalues,2,`iloop') * `DP'
		local _num = `_num' + `n1'
	// Calcul du dénominateur 	 
		// DP = nombre de control exposé 
		local DP = `iloop' - 1
	    local n2 = el(_mvalues,1,`iloop') * `DP' 
		local _den = `_den' + `n2'	 
	}
	// that's all
}

display
display  as text %45s "Total discordant pairs for exposed case :" as result %9.0g `_num'
display  as text %45s "Total discordant pairs for unexposed case :" as result %9.0g `_den'

display 

// Rothman method with Mac Nemar Chi2
local chi = (`_num'-`_den')/sqrt(`_num'+`_den')
local chi2 = `chi'^2 
local pchi2 = chi2tail(1,`chi2')
// Mac nemar corrected for continuity
local chicor = (abs(`_num'-`_den')-1)/sqrt(`_num'+`_den')
local chi2cor = `chicor'^2
local pchi2cor = chi2tail(1,`chi2cor')
local or = `_num' / `_den'	
// SD for ln(OR) = sqrt(1/`_num'+1/`_den')  // (Rothman)
local orupb = exp(ln(`or')+1.96*(sqrt(1/`_num'+1/`_den')))   
local orlob = exp(ln(`or')-1.96*(sqrt(1/`_num'+1/`_den')))   


// Exact method (for 2*2 tables only : ci for comparing two proportions) 	
local t=`_num'+`_den'
quietly cii `t' `_num', level(`level')
local orexlb = r(lb)/(1-r(lb))
local orexub = r(ub)/(1-r(ub))
local note "(exact)"
//

di
di  in gr %34s "Matched odds ratio " in ye %9.3g `or' _skip(5) _continue
 
if `maxcontrol' == 1 {	
	di in ye "[ " %10.0g `orlob' " - " %12.0g `orupb' in gr " ]"
	di  in gr %34s "Exact CI " _col(49) in ye /*
	*/ "[ " %10.0g `orexlb' " - " %12.0g `orexub' in gr " ]"
di
di in gr  %34s "McNemar's chi2(1) = " in ye %9.2f `chi2' /*
	*/ _skip(4) in gr "Prob > chi2 =" in ye %7.4f `pchi2'
di in gr  %34s "McNemar's corrected chi2(1) = " in ye %9.2f `chi2cor' /*
	*/ _skip(4) in gr "Prob > chi2 =" in ye %7.4f `pchi2cor'
	
} 
else {
di 
di
clogit `1' `2' , group(`by') or nolog
}
	
// cleanup
restore
capture matrix drop _cpcvalues

end


