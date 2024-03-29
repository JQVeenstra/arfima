% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lARFIMA.R
\name{lARFIMA}
\alias{lARFIMA}
\title{Exact log-likelihood of a long memory model}
\usage{
lARFIMA(
  z,
  phi = numeric(0),
  theta = numeric(0),
  dfrac = numeric(0),
  phiseas = numeric(0),
  thetaseas = numeric(0),
  dfs = numeric(0),
  H = numeric(0),
  Hs = numeric(0),
  alpha = numeric(0),
  alphas = numeric(0),
  period = 0,
  useC = 3
)
}
\arguments{
\item{z}{A vector or (univariate) time series object, assumed to be (weakly)
stationary.}

\item{phi}{The autoregressive parameters in vector form.}

\item{theta}{The moving average parameters in vector form.  See Details for
differences from \code{\link{arima}}.}

\item{dfrac}{The fractional differencing parameter.}

\item{phiseas}{The seasonal autoregressive parameters in vector form.}

\item{thetaseas}{The seasonal moving average parameters in vector form.  See
Details for differences from \code{\link{arima}}.}

\item{dfs}{The seasonal fractional differencing parameter.}

\item{H}{The Hurst parameter for fractional Gaussian noise (FGN).  Should
not be mixed with \code{dfrac} or \code{alpha}: see "Details".}

\item{Hs}{The Hurst parameter for seasonal fractional Gaussian noise (FGN).
Should not be mixed with \code{dfs} or \code{alphas}: see "Details".}

\item{alpha}{The decay parameter for power-law autocovariance (PLA) noise.
Should not be mixed with \code{dfrac} or \code{H}: see "Details".}

\item{alphas}{The decay parameter for seasonal power-law autocovariance
(PLA) noise.  Should not be mixed with \code{dfs} or \code{Hs}: see
"Details".}

\item{period}{The periodicity of the seasonal components.  Must be >= 2.}

\item{useC}{How much interfaced C code to use: an integer between 0 and 3.
The value 3 is strongly recommended. See "Details".}
}
\value{
The exact log-likelihood of the model given with respect to z, up to
an additive constant.
}
\description{
Computes the exact log-likelihood of a long memory model with respect to a
given time series.
}
\details{
The log-likelihood is computed for the given series z and the parameters.
If two or more of \code{dfrac}, \code{H} or \code{alpha} are present and/or
two or more of \code{dfs}, \code{Hs} or \code{alphas} are present, an error
will be thrown, as otherwise there is redundancy in the model.  Note that
non-seasonal and seasonal components can be of different types: for example,
there can be seasonal FGN with FDWN at the non-seasonal level.

The moving average parameters are in the Box-Jenkins convention: they are
the negative of the parameters given by \code{\link{arima}}.

For the useC parameter, a "0" means no C is used; a "1" means C is only used
to compute the log-likelihood, but not the theoretical autocovariance
function (tacvf); a "2" means that C is used to compute the tacvf and not
the log-likelihood; and a "3" means C is used to compute everything.

Note that the time series is assumed to be stationary: this function does
not do any differencing.
}
\examples{

set.seed(3452)
sim <- arfima.sim(1000, model = list(phi = c(0.3, -0.1)))
lARFIMA(sim, phi = c(0.3, -0.1))

}
\references{
Box, G. E. P., Jenkins, G. M., and Reinsel, G. C. (2008) Time
Series Analysis: Forecasting and Control.  4th Edition. John Wiley and Sons,
Inc., New Jersey.

Veenstra, J.Q. Persistence and Antipersistence:  Theory and
Software (PhD Thesis)
}
\seealso{
\code{\link{arfima}}

\code{\link{lARFIMAwTF}}

\code{\link{tacvfARFIMA}}
}
\author{
Justin Veenstra
}
\keyword{ts}
