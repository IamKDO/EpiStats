\name{CC}
\alias{CC}
\alias{cc}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{Univariate analysis of case control studies}

\description{
    CC is used with case-control studies to determine the association between an exposure and an outcome.
    Note that all variables need to be numeric and binary and coded as "0" and "1".
    Point estimates and confidence intervals for the odds ratio are calculated, along with
    attributable or prevented fractions for the exposed and total population.

    Additionally you can select if you want to display the Fisher's exact test, by specifying exact = TRUE.

If you specify full = TRUE you can easily access useful statistics from the output tables.}

\usage{
CC(data, cases, exposure, exact = FALSE, full = FALSE, title = "CC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data.frame}
  \item{cases}{character - Case variable}
  \item{exposure}{character - Exposure variable}
  \item{exact}{boolean - TRUE if you would like to display Fisher's exact p-value}
  \item{full}{boolean - TRUE if you need to display useful statistics and values for formatting}
  \item{title}{character - title of tables}
}

\value{list:
  \item{df1}{data.frame -  two by two table}
  \item{df2}{data.frame - statistics}
  \item{df1.align}{character - alignment for kable/xtable}
  \item{df2.align}{character - alignment for kable/xtable}
  \item{df1.digits}{integer vector - digit number displayed for kable/xtable}
  \item{df2.digits}{integer vector - digit number displayed for kable/xtable}
  \item{st}{list - individual statistics}

  The item st returns the odds ratio and its 95 percent confidence intervals, the attributable fraction among the exposed and its 95 percent confidence intervals, the attributable fraction among the population and its 95 percent confidence intervals, the Chi square value, the Chi square p-value and the Fisher's exact test p-value.
}

\references{Stata 13: cc https://www.stata.com/manuals13/stepitab.pdf}
\author{jean.pierre.decorps@gmail.com}

\note{You can use the lowercase command "cc" in place of "CC"

Please note also that when the outcome is frequent the odds ratio will overestimate the risk ratio (if OR>1) or underestimate the risk ratio (OR<1). If the outcome is rare, the risk ratio and the odds ratio are similiar.

In a case control study, the attributable fraction among the exposed and among the population assume that the OR approximates the risk ratio.

Please interpret all measures with caution.
}

\seealso{
CCTable, CCInter, CS, CSTable, CSInter
}

\examples{
library(EpiStats)

# Dataset by Anja Hauri, RKI.
data(Tiramisu)
DF <- Tiramisu

# The CC command looks at the association between the outcome variable "ill"
# and an exposure "mousse"

CC(DF, "ill", "mousse")

# The option exact = TRUE provides Fisher's exact test p-values
CC(DF, "ill", "mousse", exact = TRUE)

# With the option full = TRUE you can easily use individual elements of the results:
result <- CC(DF, "ill", "mousse", full = TRUE)
result$st$odds_ratio$point_estimate

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~documentation }
