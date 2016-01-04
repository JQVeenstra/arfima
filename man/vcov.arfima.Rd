\name{vcov.arfima}
\alias{vcov.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts the Variance-Covariance Matrix
}
\description{
Extracts the variance-covariance matrices (one or two for each mode) from a fitted \code{arfima} object.
}
\usage{
\method{vcov}{arfima}(object, type = c("b", "o", "e"), cor = FALSE, 
digits = max(4, getOption("digits") - 3), tapprox = FALSE, summ = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A fitted \code{arfima} object
}
  \item{type}{
Which type of covariance matrix to return:  "o" is the observed matrix (from solving the Hessian), "e" is the expected matrix (from solving the information matrix), 
and "b" is both.
}
  \item{cor}{
Whether or not the correlation matrix should be returned instead.
}
  \item{digits}{
The number of digits to print. 
}
  \item{tapprox}{
Whether or not to use an approximation to find the expected matrix.  Highly recommended to be FALSE, as it takes much longer, and is an approximation.
}
  \item{summ}{
Whether the call is from the \code{summary.arfima} function.  Should not be used except internally.
}
  \item{\dots}{
Optional arguments, currently not used. 
}
}

\value{
A list of lists (one for each mode) with components \code{observed} and/or \code{expected}.
}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
Justin Veenstra
}

\seealso{
\code{\link{summary.arfima}}, \code{\link{arfima}}
}
\examples{
\dontrun{
set.seed(1234)
sim <- arfima.sim(1000, model = list(dfrac = 0.4, phi = .8, theta = -0.5))
fit1 <- arfima(sim, order = c(1, 0, 1)) 
fit2 <- arfima(sim, order = c(1, 0, 1), lmodel = "g")
fit3 <- arfima(sim, order = c(1, 0, 1), lmodel = "h")
fit1
fit2
fit3
vcov(fit1)
vcov(fit2)
vcov(fit2)
}
}

\keyword{ ts}
