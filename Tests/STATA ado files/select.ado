// select.do 
// simulate EpiInfo select command
// Gilles Desvé, 2007

capture program drop select
program select
syntax [anything(name=select)] 


quietly {
if "`select'" != "" {
if "$preserved" == "" {
   save episelect,replace 
   global preserved = "1"
}
noisily keep if `select' 
} 

if "`select'" == "" {
use episelect, clear
capture macro drop preserved
preserve
}

}

end
 