capture program drop freq
program freq

syntax varlist [if] [in] [,NOLabel]

if "`nolabel'" != "" {
local nolabel = ", nolabel"
}

quietly tab1 `varlist' `if' `in',missing
local _withmiss = r(N)
tab1 `varlist' `if' `in' `nolabel' 

display
display "Missing values : " `_withmiss' -r(N)

end


