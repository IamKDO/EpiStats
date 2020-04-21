* decomp - differentiation - trend analysis
program decompdiff

syntax varlist [, diff(string)]

capture drop t_*
graph drop _all

local time = "`2'"
local count = "`1'"

tsset `time'


if "`diff'" == "linear" {
***************************************************
* differentiation trend analysis *
***************************************************

ac `count', lag(104) name(auto)

***************************************************
* Linear trend: first order lag 1 differentiation *
* D.var is the diff between var and var-1
***************************************************
tsline `count' D.`count' , name(line1)
graph export "${PICT}Line plot trend using first order differentiation.png" , replace
ac D.`count' , lag(104) name(auto1)
corrgram D.`count' 

*******************************************************
* Quadratic trend: second order lag 1 differentiation *
* D2.var is order 2 differenciation
*******************************************************
tsline D.`count'  D2.`count' , name(line2)
ac D2.`count', lag(104) name(auto2)
graph export "${PICT}Autocorrelation graph using differentiation.png" , replace
corrgram D2.`count'
}

if "`diff'" == "season" {

***************************************
* First order lag “i” differentiation *
***************************************
gen d1=D.`count'
gen t_d_li=d1-d1[_n-16]-d1[_n-52]
tsline count d1 t_d_li
ac t_d_li, lag(104)
corrgram t_d_li

drop t_*

// gen res2=t_d_li

}


end


