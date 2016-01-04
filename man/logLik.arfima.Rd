\name{logLik.arfima}
\alias{logLik.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Log-Likelihood Values
}
\description{
Extracts log-likelihood values from a \code{arfima} fit.
}
\usage{
\method{logLik}{arfima}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A fitted \code{arfima} object
}
  \item{\dots}{
Optional arguments not currently used. 
}
}
\details{
Uses the function \code{DLLoglikelihood} from the package \code{\link{ltsa}}.  The log-likelihoods returned are exact up to an additive constant.  
}
\value{
A vector of log-likelihoods, one for each mode, is returned, along with the degrees of freedom. 
}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
Justin Veenstra
}
\seealso{
\code{\link{AIC.arfima}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }