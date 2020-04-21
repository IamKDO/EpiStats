*! version 4  12Feb2010
program esampsi
	version 9, missing

	gettoken a 0 : 0, parse(" ,")
	gettoken b 0 : 0, parse(" ,")
	gettoken c 0 : 0, parse(" ,")
	gettoken d 0 : 0, parse(" ,")
	gettoken e 0 : 0, parse(" ,")


syntax [, cc cs pop pc pcc pcco]

if "`a'"=="" {
    db esampsi
	exit
}

if "`cs'" == "cs" {
local b = `b'/100	
local a = `b' * `a'	
	
quietly sampsi `b' `a' , p(`e') r(`c') alpha(`d')
local tot = r(N_1) + r(N_2)

display as text "Sample size :"
display as text "Number of exposed   : " as result r(N_2)
display as text "Number of unexposed : " as result r(N_1)
display as text "Total number        : " as result `tot'
}

if "`cc'" == "cc" {
local b = `b'/100	
// P2={(OR x P1)/(1+ (OR-1) x P1)}
local a = (`a'*`b')/ ( 1 + (`a'-1) * `b' )	
quietly sampsi `b' `a' , p(`e') r(`c') alpha(`d')
local tot = r(N_1) + r(N_2)

display as text "Sample size :"
display as text "Number of cases    : " as result r(N_1)
display as text "Number of controls : " as result r(N_2)
display as text "Total number       : " as result `tot'
}

if "`pc'" == "pc" {
local b = `b'/100	
local a = `a'/100	
sampsi  `a' `b', n1(`d') n2(`e')

}

if ("`pcc'" == "pcc") | ("`pcco'" == "pcco") {
local a = `a'/100
if ("`pcco'" == "pcco") {
local b = (`a' * `b')/ (1+((`b'-1)*`a'))
}
else {
local b = `b'/100
}
local e = `d' * `e' 	
sampsi  `b' `a' , n1(`d') n2(`e')

}


if "`pop'" == "pop" {
local b = `b'/100	
local a = `a' / 100 	
if `a' > 0.50 {
local b = `a' - `b'
}
else {
local b = `a' + `b'
}
quietly sampsi `a' `b' , onesample p(0.5) alpha(`d')
local n1 = r(N_1)* `c'

local n1 = `n1' * (`e' / (`n1'+`e'-1) )

display as text "Sample size :"
display as text "Total number        : " as result %9.0f `n1'
}

end	

