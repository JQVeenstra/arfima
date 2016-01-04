\name{print.summary.arfima}
\alias{print.summary.arfima}

\title{
Prints the output of a call to \code{summary} on an \code{arfima} object
}
\description{
Prints the output of a call to \code{summary} on an \code{arfima} object
}

\usage{
\method{print}{summary.arfima}(x, digits = max(6, getOption("digits") - 3),
signif.stars = getOption("show.signif.stars"), ...)
}

\arguments{
  \item{x}{
A \code{summary.arfima} object
}
  \item{digits}{
The number of digits to print
}
  \item{signif.stars}{
Whether to print stars on significant output
}
  \item{\dots}{
Currently not used
}
}

\value{
Returns the object \code{x} invisibly
}
\references{
Veenstra, J. and McLeod, A. I. (Working Paper).
The arfima R package: Exact Methods for Hyperbolic Decay Time Series

}
\author{
Justin Veenstra
}

\seealso{
\code{\link{arfima}}, \code{\link{print.arfima}}, \code{\link{summary.arfima}}, \code{\link{print}}
}
\examples{
\dontrun{
set.seed(54678)
sim <- arfima.sim(1000, model = list(phi = 0.9, H = 0.3))
fit <- arfima(sim, order = c(1, 0, 0), lmodel = "g")
summary(fit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts}
