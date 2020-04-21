* isoweek.ado, for Stata version 8.0 or higher
*! version 1  9Feb2009
*
* syntax: isoweek datevar [,generate(varname) format(dmy/mdy)]
* options: generate  - variable to create the handle isoweek; default = _isoweek
*           format   - if datevar is a string, format to be used to transform to date
*
* Adapted by Bernadette Gergonne from 
* a formula developed by Evert van den Heuvel for Excel. 
* http://www.cpearson.com/excel/weeknum.htm
*
* works at least from 1960 to 2010. 
*
* Ado programming by Gilles DESVE
* (Opinions, comments are welcome as well!)
* in case of any questions please contact us: 
* g.desve@epiconcept.fr

program isoweek
version 9

syntax varname [,generate(string) format(string)  fill nb ]

tokenize `varlist'

if "`format'" == ""{
local format = "dmy"
}

local PREFIX = "`generate'"
if "`generate'"=="" {
local PREFIX = "_iso"
}
capture drop `PREFIX'week
capture drop `PREFIX'weeklabel
capture drop `PREFIX'weekyear
capture drop `PREFIX'weeknb
capture drop `PREFIX'wy
capture drop  _thedate

// generate "DDD" a stata-date from a string date ("DATE") of the form (dd/mm/yyyy)
// stata date = the number of days since the 1st January 1960
local datype : type `1'
if substr("`datype'",1,3) == "str" { 
	gen _thedate = date(`1',"`format'") 
}
else {
	gen _thedate = `1'
}
format _thedate %d
if "`fill'" != "" {
    preserve
    contract _thedate, freq(_temp_nb)
	tsset _thedate 
	tsfill
	keep if missing(_temp_nb)
	drop _temp_nb
	save _tempweek ,replace emptyok
	restore 
	append using _tempweek 
	capture erase _tempweek
}
// gen ISO = 1+int((DDD-mdy(1,5,year(DDD+4-(dow(DDD+6)+1)))+(dow(mdy(1,3,
// year(DDD+4-(dow(DDD+6)+1))))+1))/7)


// GENERATION OF THE ISO WEEK !!
gen `PREFIX'week = 1 +int((_thedate- mdy(1,5,year(_thedate+4-(dow(_thedate+6)+1))) /// 
     +(dow(mdy(1,3, year(_thedate+4-(dow(_thedate+6)+1))))+1))/7)

gen `PREFIX'weekyear=year(_thedate)
replace `PREFIX'weekyear=year(_thedate)-1 if month(_thedate)==1 & (`PREFIX'week==52 | `PREFIX'week==53 )
replace `PREFIX'weekyear=year(_thedate)+1 if month(_thedate)==12  & `PREFIX'week==1

// gen `PREFIX'wy = string(`PREFIX'weekyear)+"-"+string(`PREFIX'week,"%02.0f")
gen `PREFIX'day = dow(_thedate)
gen `PREFIX'date = _thedate-cond(dow(_thedate)==0,6 ,dow(_thedate)-1) 
format `PREFIX'date %d

// now we can drop except Monday in case we use fill option
if "`fill'" != "" {
   drop if `1'==. & `PREFIX'day != 1
   // some extra dates will be generated for weeks 52-53    
   bys `PREFIX'date : egen `PREFIX'count = count(`PREFIX'date)
   drop if `1'==. & `PREFIX'count > 1
   drop `PREFIX'count
}
drop `PREFIX'day

gen `PREFIX'weeklab = string(`PREFIX'weekyear) + "w" + string(`PREFIX'week)


if "`nb'"=="nb" {
    sort `PREFIX'date 
	gen `PREFIX'weeknb =.
	replace `PREFIX'weeknb=1 in 1 
	replace  `PREFIX'weeknb = cond(`PREFIX'date == `PREFIX'date[_n-1],`PREFIX'weeknb[_n-1],`PREFIX'weeknb[_n-1]+1) if _n > 1 
}
drop _thedate 
*/

end




