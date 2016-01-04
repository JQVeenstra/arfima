\name{iARFIMA}
\alias{iARFIMA}

\title{
	The Fisher information matrix of an ARFIMA process
}
\description{
	Computes the approximate or (almost) exact Fisher information matrix of an ARFIMA process
}
\usage{
iARFIMA(phi = numeric(0), theta = numeric(0),
phiseas = numeric(0), thetaseas = numeric(0), 
period = 0, dfrac = TRUE, dfs = FALSE, exact = TRUE)
}

\arguments{
    \item{phi}{
	The autoregressive parameters in vector form. 
}
  \item{theta}{
	The moving average parameters in vector form.  See Details for differences from \code{\link{arima}}.
}

  \item{phiseas}{
	The seasonal autoregressive parameters in vector form. 
}

  \item{thetaseas}{
	The seasonal moving average parameters in vector form.  See Details for differences from \code{\link{arima}}.
}

  \item{period}{
	The periodicity of the seasonal components.  Must be >= 2.
}

  \item{dfrac}{
	TRUE if we include the fractional d parameter, FALSE otherwise
}
  \item{dfs}{
	TRUE if we include the seasonal fractional d parameter, FALSE otherwise
}
  \item{exact}{
	If FALSE, calculate the approximate information matrix via psi-weights.  Otherwise the (almost) exact information matrix will be calculated.  See "Details".
}


}
\details{
	The matrices are calculated as outlined in Veenstra and McLeod (2012), which draws on many references.  The psi-weights approximation has a fixed maximum lag for the
	weights as 2048 (to be changed to be adaptable.)  The fractional difference(s) by AR/MA components have a fixed maximum lag of 256, also to be changed.  Thus
	the exact matrix has some approximation to it.  Also note that the approximate method takes much longer than the "exact" one.

	The moving average parameters are in the Box-Jenkins convention: they are the negative of the parameters given by \code{\link{arima}}. 
}
\value{
	The information matrix of the model.
}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
	Justin Veenstra
}

\seealso{
\code{\link{IdentInvertQ}}
}
\examples{
\dontrun{
tick <- proc.time()
exactI <- iARFIMA(phi = c(.4, -.2), theta = c(.7), phiseas = c(.8, -.4),
	d = TRUE, dfs = TRUE, period = 12)
proc.time() - tick
tick <- proc.time()
approxI <- iARFIMA(phi = c(.4, -.2), theta = c(.7), phiseas = c(.8, -.4), 
	d = TRUE, dfs = TRUE, period = 12, exact = FALSE)
proc.time() - tick
exactI
max(abs(exactI - approxI))
}
}

\keyword{ ts }