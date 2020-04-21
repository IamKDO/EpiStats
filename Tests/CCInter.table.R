# Trying EpiStats with factors

# Installing the packages I need
library(EpiStats)

# Opening the csv file - Stegen outbreak
DF <- read.csv("stegen1.csv", sep=",", stringsAsFactors=FALSE)

# Creating factors for ill and tira
DF$tira1<-factor(DF$tira, levels = c(0,1), labels = c("no", "yes"))
DF$ill1<-factor(DF$ill,  levels = c(0,1), labels = c("not ill", "ill"), ordered=TRUE)
DF$mousse1<-factor(DF$mousse, levels = c(0,1), labels = c("no", "yes"))
DF$redjelly1<-factor(DF$redjelly, levels = c(0,1), labels = c("no", "yes"))


# CCInter With integers
CCInter(DF,"ill",exposure="redjelly",by="tira")
# CCinter with table option
CCInter(DF,"ill",exposure="redjelly",by="tira", table = TRUE)

# CSInter With factors
CCInter(DF,"ill1",exposure="redjelly1",by="tira1", table = TRUE)

CCInter(DF,"ill1",exposure="beer",by="tira1", table = TRUE)
# local _inter = (`_rr10' -1) + (`_rr01' - 1) + 1
# inter = (`_rr11' - 1 ) - (`_rr10' -1) - (`_rr01' - 1)

CCInter(DF,"ill1",exposure="mousse1",by="tportion", table = TRUE)

