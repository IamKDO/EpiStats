# Trying EpiStats with factors

# Installing the packages I need
library(EpiStats)

# Opening the csv file - Stegen outbreak
DF <- read.csv("stegen1.csv", sep=",", stringsAsFactors=FALSE)

# Creating factors for ill and tira
DF$tira1<-factor(DF$tira, levels = c(0,1), labels = c("no", "yes"))
DF$ill1<-factor(DF$ill,  levels = c(0,1), labels = c("not ill", "ill"), ordered=TRUE)
DF$mousse1<-factor(DF$mousse, levels = c(0,1), labels = c("no", "yes"))


# CSInter With integers
CSInter(DF,"ill",exposure="mousse",by="tira")
# CSInter With factors
CSInter(DF,"ill1",exposure="mousse1",by="tira1")

CSInter(DF,"ill1",exposure="beer",by="tira1", table = TRUE)
# local _inter = (`_rr10' -1) + (`_rr01' - 1) + 1
# inter = (`_rr11' - 1 ) - (`_rr10' -1) - (`_rr01' - 1)

CSInter(DF,"ill",exposure="mousse",by="tportion", table = TRUE)
