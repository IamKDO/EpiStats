\name{CCTable}
\alias{CCTable}
\alias{cctable}

\title{
Summary table for univariate analysis of case control studies
}
\description{
CCTable is used for univariate analysis of case control studies with several exposures. The results are summarised in one table with one row per exposure making comparisons between exposures easier and providing a useful table for integrating into reports. Note that all variables need to be numeric and binary and coded as "0" and "1".

The results of this function contain: The name of exposure variables, the total number of cases, the number of exposed cases, the percentage of exposed among cases, the number of controls, the number of exposed controls, the percentage of exposed among controls, odds ratios, 95\%CI intervals, p-values.

You can optionally choose to display the Fisher's exact p-value instead of the Chi squared p-value, with the option exact = TRUE.

You can specify the sort order, with the option sort = "or" to order by odds ratios. The default sort order is by p-values.

The option full = TRUE provides you with useful formatting information, which can be handy if you're using "markdown".
}
\usage{
CCTable(x, cases, exposure = c(), exact = FALSE, sort = "pvalue", full = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ data.frame}
  \item{cases}{ character - cases binary variable (0 / 1)}
  \item{exposure}{character vector - exposure variables}
  \item{exact}{ boolean - TRUE if you want the Fisher's exact p-value instead of CHI2}
  \item{sort}{ character - [pvalue, or, pe] sort by pvalue (default) or by odds ratio, or by percent
  exposed}
  \item{full}{boolean - TRUE if you need to display useful values for formatting}
}
\details{
The results of this function contain: The name of exposure variables, the total number of cases, the number of exposed cases, the percentage of exposed among cases, the number of controls, the number of exposed controls, the percentage of exposed among controls, odds ratios, 95\%CI intervals, p-values.

You can optionally choose to display the Fisher???s exact p-value instead of the Chi squared p-value, with the option exact = TRUE.

You can specify the sort order, with the option sort=???or??? to order by odds ratios. The default sort order is by p-values.

The option "full = TRUE" provides you with useful formatting information, which can be handy if you're using "markdown".
}
\value{list :
  \item{df}{data.frame - results table}
  \item{digits}{integer vector - digit number displayed for kable/xtable}
  \item{align}{character - alignment for kable/xtable}
}
\note{
  - You can use the lowercase command "cctable" instead of "CCTable"
}

\references{
cctable for Stata by *Gilles Desve* and *Peter Makary*.
}
\author{jean.pierre.decorps@gmail.com}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{CC, CCInter}

\examples{
library(EpiStats)

data(Tiramisu)
df <- Tiramisu

# You can see the association between several exposures and being ill.
cctable(df, "ill", exposure=c("sex", "wmousse", "tira", "beer", "mousse"))

# By storing results in res, you can also use individual elements of the results.
# For example if you would like to view a particular odds ratio,
# you can view it by typing (for example):

res = CCTable(df, "ill", exposure = c("sex", "wmousse", "tira", "beer", "mousse"), exact=TRUE)
res$df$`Odds Ratio`[1]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~documentation }
