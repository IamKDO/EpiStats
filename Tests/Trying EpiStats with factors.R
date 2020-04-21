# Trying EpiStats with factors

# Installing the packages I need
library(EpiStats)

# Opening the csv file - Stegen outbreak
DF <- read.csv("stegen1.csv", sep=",", stringsAsFactors=FALSE)

# Creating factors for ill and tira
DF$tira1<-factor(DF$tira, levels = c(0,1), labels = c("no", "yes"))
DF$ill1<-factor(DF$ill,  levels = c(0,1), labels = c("not ill", "ill"), ordered=TRUE)
DF$mousse1<-factor(DF$mousse, levels = c(0,1), labels = c("no", "yes"))

## This works well:
# CC with integers
CC(DF,"ill","tira")
# CC with factors
CC(DF,"ill1","tira1")
# CCTable with integers
CCTable(DF,"ill", exposure=c("tira","mousse","beer"), exact=TRUE)
# CCTable with factors
CCTable(DF,"ill1", exposure=c("tira1","mousse1","beer"), exact=TRUE)

## This does not work with factors
# With integers
CCInter(DF,"ill",exposure="mousse",by="tira")
# With factors
CCInter(DF,"ill1",exposure="mousse1",by="tira1")



# CS with integers
CS(DF,"ill","tira")
# CC with factors
CS(DF,"ill1","tira1")

# CSTable with integers
CSTable(DF,"ill", exposure=c("tira","mousse","beer"), exact=TRUE)
# CCTable with factors
CSTable(DF,"ill1", exposure=c("tira1","mousse1","beer"), exact=TRUE)

# CSInter With integers
CSInter(DF,"ill",exposure="mousse",by="tira")
# CSInter With factors
CSInter(DF,"ill1",exposure="mousse1",by="tira1")

