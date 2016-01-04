\name{bestModes}
\alias{bestModes}

\title{
Finds the best modes of an \code{arfima} fit.
}

\description{
Finds the best modes of an \code{arfima} fit with respect to log-likelihood.
}

\usage{
bestModes(object, bestn)
}

\arguments{
  \item{object}{
An object of class "arfima".
}
  \item{bestn}{
The top number of modes to keep with respect to the log-likelihood.
}
}
\details{
This is the easiest way to remove modes with lower log-likelihoods.
}
\value{
The \code{bestn} "best" modes.
}

\author{
Justin Veenstra
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{arfima}}
}
\examples{
\dontrun{
set.seed(8765)
sim <- arfima.sim(1000, model = list(phi = 0.4, theta = 0.9, dfrac = 0.4))
fit <- arfima(sim, order = c(1, 0, 1))
fit
fit <- bestModes(fit, 2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }

