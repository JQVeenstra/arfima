\name{removeMode}
\alias{removeMode}

\title{
Removes a mode from an \code{arfima} fit.
}

\description{
Removes a mode from an \code{arfima} fit.
}

\usage{
removeMode(object, num)
}

\arguments{
  \item{object}{
An object of class "arfima".
}
  \item{num}{
The number of the mode as in the printed value of the object.
}
}
\details{
This function is useful if one suspects a mode is spurious and does not want to call the weed function.
}
\value{
The original object with the mode removed.
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
fit <- removeMode(fit, 3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }

