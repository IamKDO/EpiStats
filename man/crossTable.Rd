\name{crossTable}
\alias{crossTable}


\title{contingency table of 2 variables}

\description{
   Creates a contingency table of 2 variables. Percentages are optionals by row, column or both.
   It can provides an optional statistic (Fisher or Chisquare).
}


\usage{ crossTable(data, var1, var2, percent="none", statistic="none")}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data.frame}
  \item{var1}{character - first varname - can be unquoted}
  \item{var2}{character - second varname - can be unquoted}
  \item{percent}{character - "none" (default) or ("row", "col", "both") - can be unquoted}
  \item{statistic}{character - "none" (default) or ("fisher", "chi2") - can be unquoted}
}

\value{data.frame - contingency table}

\author{jean.pierre.decorps@gmail.com}


\seealso{
orderFactors, CC, CS
}

\examples{
library(EpiStats)

# Dataset by Anja Hauri, RKI.
data(Tiramisu)
DF <- Tiramisu

# Table with percentagges and statistic on ordered factors
DF \%<>\%
  orderFactors(ill , values = c(1,0), labels = c("YES", "NO")) \%>\%
  orderFactors(sex, values = c("males", "females"), labels = c("Males", "Females"))

crossTable(DF, "ill", "sex", "both", "chi2")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~documentation }
