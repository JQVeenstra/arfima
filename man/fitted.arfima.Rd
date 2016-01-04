\name{fitted.arfima}
\alias{fitted.arfima}

\title{
Extract Model Fitted Values
}
\description{
Extract fitted values from an \code{arfima} object. 
}
\usage{
\method{fitted}{arfima}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{arfima} object.
}
  \item{\dots}{
Optional parameters.  Currently not used.
}
}

\value{
A list of vectors of fitted values, one for each mode.
}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
Justin Veenstra
}

\seealso{
\code{\link{arfima}}, \code{\link{resid.arfima}}
}
\examples{
\dontrun{
set.seed(8564)
sim <- arfima.sim(1000, model = list(phi = c(0.2, 0.1), dfrac = 0.4, theta = 0.9))
fit <- arfima(sim, order = c(2, 0, 1))

fit

resid <- resid(fit)
par(mfrow = c(1, 3))
fitted <- fitted(fit)
plot(fitted[[1]], resid[[1]])
plot(fitted[[2]], resid[[2]])
plot(fitted[[3]], resid[[3]])
par(mfrow = c(1, 1))
}
}
\keyword{ ts}