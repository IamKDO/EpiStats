\name{CSInter}
\alias{CSInter}
\alias{csinter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stratified analysis for cohort studies measuring risk}
\description{CSInter is useful to determine the effects of a third variable on the association between an exposure and an outcome.  CSInter produces 2 by 2 tables with stratum specific risk ratios, attributable risk among exposed and population attributable risk. Note that the outcome and exposure variable need to be numeric  and binary and coded as "0" and 1". The third variable needs to be numeric, but may have more categories, such as "0", "1" and "2".
}
\usage{
CSInter(x, cases, exposure, by, table = FALSE, full = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ data.frame }
  \item{cases}{ string: illness binary variable (0 / 1) }
  \item{exposure}{ string: exposure binary variable (0 / 1) }
  \item{by}{ string: stratifying variable (a factor) }
  \item{table}{boolean - TRUE if you need to display interaction table}
  \item{full}{boolean - TRUE if you need to display useful values for formatting}
}
\details{CSInter is useful to determine the effects of a third variable on the association between an exposure and an outcome.  CSInter produces 2 by 2 tables with stratum specific risk ratios, attributable risk among exposed and population attributable risk. Note that the outcome and exposure variable need to be numeric  and binary and coded as "0" and 1". The third variable needs to be numeric, but may have more categories, such as "0", "1" and "2".

CSInter displays a summary with the crude RR, the Mantel Haenszel adjusted RR and the result of a
  "Woolf" test for homogeneity of stratum-specific RR.

The option full = TRUE provides you with useful formatting information, which can be handy if you're using "markdown".
}
\value{list:
  \item{df1}{data.frame - cross-table}
  \item{df2}{data.frame - statistics}
  \item{df1.digits}{integer vector - digit number displayed for kable/xtable}
  \item{df2.digits}{integer vector - digit number displayed for kable/xtable}
}
\references{csinter for Stata by *Gilles Desve*}

\author{jean.pierre.decorps@gmail.com}
\note{
- You can use the lowercase command "csinter" instead of "CSInter"
- The "by" variable (the stratifying variable) can have more than 2 levels
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
CS, CSTable
}

\examples{
library(EpiStats)

data(Tiramisu)
DF <- Tiramisu

# Here you can see the association between wmousse and ill for each stratum of tira:
csinter(DF, "ill", "wmousse", by = "tira")

# By storing the results in the object "res", you can use individual elements
# of the results. For example if you would like to view just the Mantel-Haenszel
# risk ratio for beer adjusted for tportion, you can view it by typing:
res <- CSInter(DF, "ill", "beer", "tportion", full = TRUE)
res$df2$Stats[3]

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~documentation }
