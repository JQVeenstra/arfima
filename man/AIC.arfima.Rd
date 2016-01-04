\name{AIC.arfima}
\alias{AIC.arfima}
\alias{BIC}
\alias{BIC.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Information criteria for \code{arfima} objects
}
\description{
Computes information criteria for \code{arfima} objects.  See \code{\link{AIC}} for more details. 
}
\usage{
\method{AIC}{arfima}(object, ..., k = 2)
\method{BIC}{arfima}(object, ...)
}

\arguments{
  \item{object}{
An object of class "arfima".  Note these functions can only be called on one object at a time because of possible multimodality. 
}
  \item{\dots}{
Other models fit to data for which to extract the AIC/BIC.  Not recommended, as an \code{arfima} object can be multimodal. 
}

  \item{k}{
The penalty term to be used.  See \code{\link{AIC}}.
}
}

\value{
The information criteria for each mode in a vector.
}

\author{
Justin Veenstra
}


\examples{
\dontrun{
set.seed(34577)
sim <- arfima.sim(500, model = list(theta = 0.9, phi = 0.5, dfrac = 0.4))
fit1 <- arfima(sim, order = c(1, 0, 1), cpus = 2)
fit2 <- arfima(sim, order = c(1, 0, 1), cpus = 2, lmodel = "g")
fit3 <- arfima(sim, order = c(1, 0, 1), cpus = 2, lmodel = "h")

AIC(fit1)
AIC(fit2)
AIC(fit3)
}
}
\keyword{ ts}
