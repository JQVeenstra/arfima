\name{predict.arfima}
\alias{predict.arfima}

\title{
Predicts from a fitted object.
}
\description{
Performs prediction of a fitted \code{arfima} object. Includes prediction for each mode, bootstrap predictions and intervals, and exact and limiting prediction error standard deviations.
}
\usage{
\method{predict}{arfima}(object, n.ahead = 1, 
newxreg = NULL, predint = 0.95, bootpred = TRUE, B = if (bootpred) 1000 else 0,
trex = FALSE, seed = NA, setmuhat0 = FALSE, cpus = 1, trend = NULL, ...)
}

\arguments{
  \item{object}{
A fitted \code{arfima} object
}
  \item{n.ahead}{
The number of steps ahead to predict
}
  \item{newxreg}{
If a regression fit, the new regressors
}
  \item{predint}{
For bootstrap prediction intervals, the percentiles to use
}
  \item{bootpred}{
Whether or not to generate bootstrap predictions and intervals.  Uses a bootstrap of the residuals as innovations to \code{\link{arfima.sim}}.
}
  \item{B}{
The number of bootstrap replicates to perform 
}
  \item{trex}{
Whether or not to use the exact formula for prediction error variances if there is integer differencing. \code{FALSE} is the default, and is much faster, with very little difference
}
  \item{seed}{
The seed or seeds (if multiple, must have length equal to \code{B}) for the bootstrap replicates
}
  \item{setmuhat0}{
Experimental. Sets muhat equal to zero
}
  \item{cpus}{
The number of CPUs to use for prediction. Currently not implemented
}
  \item{trend}{
An optional vector the length of \code{n.ahead} or longer to add to the predictions 
}
  \item{\dots}{
Optional arguments. Currently not used
}
}

\value{
A list of lists, one for each mode with relavent details about the prediction
}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
Justin Veenstra
}

\seealso{
\code{\link{arfima}}, \code{\link{plot.predarfima}}, \code{\link{print.predarfima}}
}
\examples{
\dontrun{
set.seed(82365)
sim <- arfima.sim(1000, model = list(dfrac = 0.4, theta=0.9, dint = 1))
fit <- arfima(sim, order = c(0, 1, 1))
fit
pred <- predict(fit, n.ahead = 5)
pred
plot(pred)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
