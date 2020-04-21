\name{Tiramisu}
\alias{Tiramisu}
\docType{data}
\title{A foodborne disease outbreak dataset}
\description{
The dataset available with the EpiStats package is from an outbreak investigation carried out in Germany in 1998 by Anja Hauri, Robert Koch Institute.
}
\usage{data(Tiramisu)}
\format{
  A data frame with 291 observations with the following 21 variables.
  \describe{
    \item{\code{ill}}{a numeric vector}
    \item{\code{dateonset}}{a date}
    \item{\code{sex}}{a factor with levels \code{females} \code{males}}
    \item{\code{age}}{a numeric vector}
    \item{\code{tira}}{a numeric vector}
    \item{\code{tportion}}{a numeric vector}
    \item{\code{wmousse}}{a numeric vector}
    \item{\code{dmousse}}{a numeric vector}
    \item{\code{mousse}}{a numeric vector}
    \item{\code{mportion}}{a numeric vector}
    \item{\code{beer}}{a numeric vector}
    \item{\code{uniquekey}}{a numeric vector}
    \item{\code{redjelly}}{a numeric vector}
    \item{\code{fruitsalad}}{a numeric vector}
    \item{\code{tomato}}{a numeric vector}
    \item{\code{mince}}{a numeric vector}
    \item{\code{salmon}}{a numeric vector}
    \item{\code{horseradish}}{a numeric vector}
    \item{\code{chickenwin}}{a numeric vector}
    \item{\code{roastbeef}}{a numeric vector}
    \item{\code{pork}}{a numeric vector}
  }
}
\references{
The dataset available with the EpiStats package is from an outbreak investigation carried out in Germany in 1998 by Anja Hauri, Robert Koch Institute. It is used in case studies by organisations including EPIET, ECDC and EpiConcept. It is provided with this package with Anja's permission.
}
\examples{
data(Tiramisu)
## maybe str(Tiramisu) ; plot(Tiramisu) ...
}
\keyword{datasets}
