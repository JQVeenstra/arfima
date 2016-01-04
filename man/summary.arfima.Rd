\name{summary.arfima}
\alias{summary.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extensive Summary of an Object
}
\description{
Provides a very comprehensive summary of a fitted \code{arfima} object.  Includes correlation and covariance matrices (observed and expected), 
the Fisher Information matrix of those parameters for which it is defined, and more, for each mode.
}
\usage{
\method{summary}{arfima}(object, digits = max(4, getOption("digits") - 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A fitted \code{arfima} object
}
 
  \item{digits}{
The number of digits to print
}
  \item{\dots}{
Optional arguments, currently not used.
}
}

\value{
A list of lists (one for each mode) of all relevant information about the fit that can be passed to \code{print.summary.arfima}.
}

\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}

\author{
Justin Veenstra
}

\seealso{
\code{\link{arfima}}, \code{\link{iARFIMA}}, \code{\link{vcov.arfima}}
}
\examples{
\dontrun{
data(tmpyr)

fit <- arfima(tmpyr, order = c(1, 0, 1))
fit

summary(fit)
}
}
\keyword{ ts }
