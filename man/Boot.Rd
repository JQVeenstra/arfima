\name{Boot}
\alias{Boot}
\title{Generic Bootstrap Function}
\description{
Generic function to bootstrap a fitted model.
}
\usage{
Boot(obj, R=1, ...)
}
\arguments{
  \item{obj}{fitted object }
  \item{R}{number of bootstrap replicates }
  \item{...}{optional arguments }
}
\details{
At present, the only function implemented is
\code{\link{Boot.arfima}}.
}
\value{
Parametric bootstrap simulation
}
\author{ A.I. McLeod and Justin Veenstra}


\seealso{ 
\code{\link{Boot.arfima}} 
}

\keyword{ ts}
