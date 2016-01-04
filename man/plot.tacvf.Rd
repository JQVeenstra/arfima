\name{plot.tacvf}
\alias{plot.tacvf}

\title{
Plots the output from a call to \code{tacvf}
}
\description{
Plots the theoretical autocovariance functions of the modes for a fitted \code{arfima} object
}
\usage{
\method{plot}{tacvf}(x, type = "o", pch = 20, 
xlab = NULL, ylab = NULL, main = NULL, xlim = NULL, 
ylim = NULL, tacf = FALSE, maxlag = NULL, lag0 = !tacf, ...)
}

\arguments{
  \item{x}{
A \code{tacvf} object from a call to said function
}
  \item{type}{
See \code{\link{plot}}. The default is recommended for short \code{maxlag}
}
  \item{pch}{
See \code{\link{plot}}
}
  \item{xlab}{
See \code{\link{plot}}
}
  \item{ylab}{
See \code{\link{plot}}
}
  \item{main}{
See \code{\link{plot}}
}
  \item{xlim}{
See \code{\link{plot}}
}
  \item{ylim}{
See \code{\link{plot}}
}
  \item{tacf}{
If \code{TRUE}, plots the theoretical autocorellations instead 
}
  \item{maxlag}{
The maximum lag for the plot
}
  \item{lag0}{
Whether or not to plot lag 0 of the tacvfs/tacfs.  Default \code{!tacf}. Used by \code{\link{tacfplot}}. 
}
  \item{\dots}{
Currently not used
} 
}
\details{
Only plots up to nine tacvfs. It is highly recommended that the \code{arfima} object be weeded before calling \code{tacvf}
}
\value{
None. There is a plot as output.
}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
Justin Veenstra
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tacvf}}
}
\examples{
set.seed(1234)
sim <- arfima.sim(1000, model = list(theta = 0.99, dfrac = 0.49))
fit <- arfima(sim, order = c(0, 0, 1))
plot(tacvf(fit))
plot(tacvf(fit), tacf = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }