\name{orderFactors}
\alias{orderFactors}


\title{Generates ordered factors.}

\description{
   Generates ordered factors for a list of columns by name or by index or range.
}


\usage{ orderFactors(data, ..., values, labels=NULL)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data.frame}
  \item{...}{character - first varname - can be unquoted}
  \item{values}{character - second varname - can be unquoted}
  \item{labels}{character - NULL (default) or ("row", "col", "both") - can be unquoted}
}

\value{data.frame - contingency table}

\author{jean.pierre.decorps@gmail.com}


\seealso{ crossTable }

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
