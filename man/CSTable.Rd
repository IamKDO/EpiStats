\name{CSTable}
\alias{CSTable}
\alias{cstable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary table for univariate analysis of cohort studies measuring risk}
\description{
CSTable is used for univariate analysis of cohort studies with several exposures. The results are summarised in one table with one row per exposure making comparisons between exposures easier and providing a useful table for integrating into reports. Note that all variables need to be numeric and binary and coded as "0" and "1".

The results of this function contain: The name of exposure variables, the total number of exposed, the number of exposed cases, the attack rate among the exposed, the total number of unexposed, the number of unexposed cases, the attack rate among the unexposed, risk ratios, 95\% percent confidence intervals, and p-values.

You can optionally choose to display the Fisher's exact p-value instead of the Chi squared p-value, with the option exact = TRUE.

You can specify the sort order, with the option sort="rr" to order by risk ratios. The default sort order is by p-values.

The option full = TRUE provides you with useful formatting information, which can be handy if you're using "markdown".

}
\usage{
CSTable(x, cases, exposure = c(), exact = FALSE, sort = "pvalue", full = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ data.frame}
  \item{cases}{ string - variable containing cases (binary 0 / 1)}
  \item{exposure}{ string vector - names of variables containing exposure (binary 0 / 1)}
  \item{exact}{ boolean - TRUE if you want the Fisher's exact p-value instead of CHI2}
  \item{sort}{character - [pvalue, rr, ar] sort by pvalue (default) or by risk ratio, or by percent of attributable risk}
    \item{full}{boolean - TRUE if you need to display useful values for formatting}
}
\details{
  The results of this function contain: The name of exposure variables, the total number
  of exposed, the number of exposed cases, the attack rate among the exposed, the total number of unexposed, the number of
  unexposed cases, the attack rate among the unexposed, risk ratios, 95%CI intervals, 95% p-values.

  You can optionally choose to display the Fisher's exact p-value instead of the Chi squared p-value, with the option exact = TRUE.

  You can specify the sort order, with the option sort="rr" to order by risk ratios. The default sort order is by p-values.

  The option full = TRUE provides you with useful formatting information, which can be handy if you're using "markdown".
}
\value{list :
  \item{df}{data.frame - results table}
  \item{digits}{integer vector - digit number displayed for kable/xtable}
  \item{align}{character - alignment for kable/xtable}
}
\references{cstable for Stata by *Gilles Desve* and *Peter Makary*}
\author{jean.pierre.decorps@gmail.com}
\note{
  - You can use the lowercase command "cstable" instead of "CSTable"
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
CS, CSInter
}
\examples{
library(EpiStats)

data(Tiramisu)
df <- Tiramisu

# You can see the association between several exposures and being ill.
CSTable(df, "ill", exposure=c("sex", "wmousse", "tira", "beer", "mousse"))

# By storing results in res, you can also use individual elements of the results.
# For example if you would like to view a particular risk ratio,
# you can view it by typing (for example):
res = CSTable(df, "ill", exposure = c("sex", "wmousse", "tira", "beer", "mousse"), exact=TRUE)
res$df$`Risk Ratio`[1]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~documentation }
