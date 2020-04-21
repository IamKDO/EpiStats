\name{CS}
\alias{CS}
\alias{cs}

\title{ Univariate analysis of cohort study measuring risk }
\description{
CS analyses cohort studies with equal follow-up time per subject. The risk (the proportion of individuals who become cases) is calculated overall and among the exposed and unexposed.
Note that all variables need to be numeric and binary and coded as "0" and "1".

Point estimates and confidence intervals for the risk ratio and risk difference are calculated, along with attributable or preventive fractions for the exposed and the total population.

Additionally you can select if you want to display the Fisher's exact test, by specifying exact = TRUE.

If you specify full = TRUE you can easily access useful statistics from the output tables.
}
\usage{
CS(x, cases, exposure, exact = F, full = FALSE, title = "CS")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data.frame}
  \item{cases}{character - Case variable}
  \item{exposure}{character - Exposure variable}
  \item{exact}{boolean - TRUE if you would like to display Fisher's exact p-value}
  \item{full}{boolean - TRUE if you need to display useful statistics and values for formatting}
  \item{title}{character - title of tables}
}

\value{list:
  \item{df1}{data.frame -  two by two table}
  \item{df2}{data.frame - statistics}
  \item{st}{list - individual statistics}
  \item{df1.digits}{integer vector - digit number displayed for kable/xtable}
  \item{df2.digits}{integer vector - digit number displayed for kable/xtable}
  \item{df2.align}{character - alignment for kable/xtable}

The item st returns the risk difference and its 95 percent confidence intervals, the risk ratio and its 95 percent confidence intervals, the attributable fraction among the exposed and its 95 percent confidence intervals, the attributable fraction among the population and its 95 percent confidence intervals, the Chi square value, the Chi square p-value and the Fisher's exact test p-value.
}

\references{Stata 13: cs. https://www.stata.com/manuals13/stepitab.pdf}
\author{jean.pierre.decorps@gmail.com}
\note{You can use the lowercase command "cs" in place of "CS"}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
CSTable, CSInter, CC, CCTable, CCInter
}
\examples{
library(EpiStats)

# Dataset by Anja Hauri, RKI.
# Dataset provided with package.
data(Tiramisu)
DF <- Tiramisu

# The CS command looks at the association between the outcome variable "ill"
# and an exposure "mousse"
CS(DF, "ill", "mousse")

# The option exact = TRUE provides Fisher's exact test p-values
CS(DF, "ill", "mousse", exact = TRUE)

# With the option full = TRUE you can easily use individual elements of the results:
result <- CS(DF, "ill", "mousse", full = TRUE)
result$st$risk_ratio$point_estimate

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~documentation }
