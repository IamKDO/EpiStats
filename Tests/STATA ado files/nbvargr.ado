*! ver 1.4    2 apr11
*! ver 1.3.2  22oct08, 1jul04  -- update to version 8
*! ver 1.2    25feb02  -- bug fix
*! ver 1.1    12jun00
program define nbvargr
  version 8.0
  syntax varlist(max=1) [if] [in] [, n(integer 10) *]
  preserve
  tokenize `varlist'
  marksample touse
  quietly keep if `touse'
  local nbvar `1'
  set more off
  
  display 
  display "Obtaining Parameter Estimates"
  display
  quietly summarize `nbvar'
  local mean = r(mean)
  local max = r(max)
  local tot = r(sum)
  local n2 = r(N)
  if `n' > `max' { 
    local n = `max' 
  }
  local n1=`n'+1
  if `n'>`n2' {
    quietly set obs `n1'
  }
  
  quietly nbreg `nbvar' 
  local alpha = e(alpha)
  tempfile nb1 pp1 op1
  
  quietly gen oprob=.
  quietly gen k = _n-1
  forvalues i=0/`n' {
    local j = `i'+1
    quietly count if `nbvar'==`i'
    quietly replace oprob = r(N)/`n2' in `j'
  }
  label variable oprob "observed proportion"
  display as txt "  Observed Proportions"
  list k oprob in 1/`n1', noobs
  
  sort k

  nbprob, m(`mean') a(`alpha') n(`n') saving(`nb1')
  merge k using "`nb1'"
  drop _merge

  pprob, m(`mean') n(`n') saving(`pp1')
  sort k
  merge k using "`pp1'"
  
  local meanstr = string(`mean',"%9.4g")
  local alphstr = string(`alpha', "%9.4g")
  label variable pprob "poisson prob"
  label variable nbprob "neg binom prob"
  quietly keep if k<=`n'
  twoway scatter oprob nbprob pprob k, connect(l l l) msym(i Dh Sh) ///
        ytitle("Proportion") ///
        b1title("mean = `meanstr'; overdispersion = `alphstr' ") `options' 
  
  set more on
end

