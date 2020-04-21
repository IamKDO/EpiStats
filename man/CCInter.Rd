\name{CCInter}
\alias{CCInter}
\alias{ccinter}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stratified analysis for case control studies}
\description{CCInter is useful to determine the effects of a third variable on the association between an exposure and an outcome.
  CCInter produces 2 by 2 tables with stratum specific odds ratios, attributable risk among
  exposed and population attributable risk.

  Note that the outcome and exposure variable need to be numeric and binary and coded as "0" and 1". The third variable needs to be numeric, but may have more categories, such as "0", "1" and "2".}
\usage{
CCInter(x, cases, exposure, by, table = FALSE, full = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ data.frame }
  \item{cases}{ string: case binary variable (0 / 1) }
  \item{exposure}{ string: exposure binary variable (0 / 1) }
  \item{by}{ string: stratifying variable (a factor) }
  \item{table}{boolean - TRUE if you need to display interaction table}
  \item{full}{boolean - TRUE if you need to display useful values for formatting}
}

\value{list:
  \item{df1}{data.frame - cross-table}
  \item{df2}{data.frame - statistics}
  \item{df1.digits}{integer vector - digit number displayed for kable/xtable}
  \item{df1.align}{character - alignment for kable/xtable}
  \item{df2.digits}{integer vector - digit number displayed for kable/xtable}
  \item{df2.align}{character - alignment for kable/xtable}

}
\details{
CCInter is useful to determine the effects of a third variable on the association between an exposure and an outcome. CCInter produces 2 by 2 tables with stratum specific odds ratios, attributable risk among exposed and population attributable risk.
Note that the outcome and exposure variable need to be numeric and binary and coded as "0" and 1". The third variable needs to be numeric, but may have more categories, such as "0", "1" and "2".
CCInter displays a summary with the crude OR, the Mantel Haenszel adjusted OR and the result of a Woolf test for homogeneity of stratum-specific OR.

The option "full = TRUE" provides you with useful formatting information, which can be handy if you're using "markdown".
}
\references{ccinter for Stata by *Gilles Desve*}

\author{jean.pierre.decorps@gmail.com}

\note{
  - You can use lowercas command "ccinter" instead of "CCInter"
  - The "by" variable (the stratifying variable) can have more than 2 levels
}


\seealso{ CC, CCTable}

\examples{
library(EpiStats)

data(Tiramisu)
DF <- Tiramisu

# Here you can see the association between wmousse and ill for each stratum of tira:
CCInter(DF, "ill", "wmousse", by = "tira")

# By storing the results in the object "res", you can use individual elements of the results.
# For example if you would like to view just the Mantel-Haenszel odds ratio for beer adjusted
# for tportion, you can view it by typing:

res <- CCInter(DF, "ill", "beer", "tportion", full = TRUE)
res$df2$Stats[3]

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~documentation }

