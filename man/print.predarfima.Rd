\name{print.predarfima}
\alias{print.predarfima}

\title{
Prints predictions and prediction intervals
}
\description{
Prints the output of \code{predict} on an \code{arfima} object
} 
\usage{
\method{print}{predarfima}(x, digits = 
max(6, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{
An object of class "predarfima"
}
  \item{digits}{
The number of digits to print
}
  \item{\dots}{
Currently not used
}
}
\details{
Prints all the relavent output of the prediction function of the \code{arfima} package
}
\value{
\code{x} is returned invisibly
}

\author{
Justin Veenstra
}

\seealso{
\code{\link{arfima}}, \code{\link{predict.arfima}}, \code{\link{predict}}, \code{\link{plot.predarfima}}
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
