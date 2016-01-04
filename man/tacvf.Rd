\name{tacvf}
\alias{tacvf}

\title{
Extracts the tacvfs of a fitted object
}

\description{
Extracts the theoretical autocovariance functions (tacvfs) from a fitted \code{arfima} or one of its modes (an \code{ARFIMA}) object.
}

\usage{
tacvf(obj, xmaxlag = 0, forPred = FALSE, n.ahead = 0, ...)
}

\arguments{
\item{obj}{
An object of class "arfima" or "ARFIMA".  The latter class is a mode of the former.
}
\item{xmaxlag}{
The number of extra points to be added on to the end.  That is, if the original series has length 300, and xmaxlag = 5, the tacvfs will go from lag 0 to lag 304.  
}  

\item{forPred}{
Should only be \code{TRUE} from a call to \code{predict.arfima}.
}  
\item{n.ahead}{
Only used internally.
}  
  \item{\dots}{
Optional arguments, currently not used. 
}
}

\value{
A list of tacvfs, one for each mode, the length of the time series.
}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
Justin Veenstra
}
 \seealso{
\code{\link{plot.tacvf}}, \code{\link{print.tacvf}}, \code{\link{tacfplot}}, \code{\link{arfima}}
 }

\keyword{ ts }
