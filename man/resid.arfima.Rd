\name{resid.arfima}
\alias{resid.arfima}
\alias{residuals.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the Residuals of a Fitted Object
}
\description{
Extracts the residuals or regression residuals from a fitted \code{arfima} object
}
\usage{
\method{resid}{arfima}(object, reg = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A fitted \code{arfima} object
}
  \item{reg}{
Whether to extract the regression residuals instead.  If \code{TRUE}, throws an error if no regression was done. 
}

  \item{\dots}{
Optional parameters.  Currently not used. 
}
}

\value{
A list of vectors of residuals, one for each mode.
}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
Justin Veenstra
}

\seealso{
\code{\link{arfima}}, \code{\link{fitted.arfima}}
}
\examples{
\dontrun{
set.seed(8564)
sim <- arfima.sim(1000, model = list(phi = c(0.2, 0.1), dfrac = 0.4, theta = 0.9))
fit <- arfima(sim, order = c(2, 0, 1))

fit

resid <- resid(fit)
par(mfrow = c(1, 3))
plot(resid[[1]])
plot(resid[[2]])
plot(resid[[3]])
fitted <- fitted(fit)
plot(fitted[[1]], resid[[1]])
plot(fitted[[2]], resid[[2]])
plot(fitted[[3]], resid[[3]])
par(mfrow = c(1, 1))
}
}

\keyword{ ts}