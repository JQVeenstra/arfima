\name{Boot.arfima}
\alias{Boot.arfima}

\title{
A parametric bootstrap based on an arfima fit.
}
\description{
Simulates a parametric bootstrap from a fitted model using either standard normal innovations or residuals from the fit.
}
\usage{
\method{Boot}{arfima}(obj, R = 1, pred = FALSE, seed = NA, ...)
}

\arguments{
  \item{obj}{
An object of class "arfima".
}
  \item{R}{
The number of bootstrap replicates to be generated (for each mode).
  }
  \item{pred}{
If TRUE, use a bootstrap sample of the residuals instead of standard normal innovations. This is for predictions.
}
  \item{seed}{
The seed to use when computing the bootstrap replicates.
}
  \item{\dots}{
Optional arguments. Not currently used.
}
}

\value{
A simulated time series with the same length as the original fitted time series is produced (for each mode) when R=1.
When R>1, a matrix with R columns is produced with each column a separate bootstrap realization (for each mode).  This is usually used by predict.

}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
Justin Veenstra and A. I. McLeod
}



\seealso{
\code{\link{Boot}}
}

\keyword{ ts }
