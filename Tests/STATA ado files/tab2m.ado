capture program drop tab2m
program  tab2m
syntax varlist  [if] [in] [ , cs * ]

 tokenize `varlist'
 local first `1'
 macro shift 1
 local rest `*' 
 
 if "`options'" != "" {
 local options = ", "+ "`options'"
 }
 
 
quietly : {
	foreach var of varlist `rest' {
		if "`cs'"=="cs" {
			noisily : tab `var' `first' `if' `in' `options'
		} 
		else {
			noisily : tab `first' `var' `if' `in' `options'
		}	
	}
}
	  
end
