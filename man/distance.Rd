\name{distance}
\alias{distance}

\title{
The distance between modes of an \code{arfima} fit.
}
\description{
The distance between modes of an \code{arfima} fit.
}
\usage{
distance(ans, p = 2, digits = 4)
}

\arguments{ 
  \item{ans}{
An object of class "arfima".
}
  \item{p}{
	The p in the p-norm to be used.
}
  \item{digits}{
The number of digits to print.
}
}

\value{
A list of two data frames: one with distances in operator space, the second with distances in the transformed (PACF) space.
}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
Justin Veensta
}



\examples{
\dontrun{
set.seed(8564)
sim <- arfima.sim(1000, model = list(phi = c(0.2, 0.1), dfrac = 0.4, theta = 0.9))
fit <- arfima(sim, order = c(2, 0, 1))

fit

distance(fit)
}
}

\keyword{ ts}
