library(EpiStats)

DF <- read.csv2("bigdata.csv")
#DF <- DF[!is.na(DF$genderbin),]

res <- CC(DF, "sorethroat", "sorethroat")

res
