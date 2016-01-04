\name{ARToPacf}
\alias{ARToPacf}

\title{
Converts AR/MA coefficients from operator space to the PACF space
}
\description{
Converts AR/MA coefficients from operator space to the PACF box-space;  usually for internal use
}
\usage{
ARToPacf(phi)
}

\arguments{
  \item{phi}{
The AR/MA coefficients in operator space
}
}

\value{
The AR/MA coefficients in the PACF space
}
\references{
Barndorff-Nielsen O. E., Schou G. (1973). "On the parametrization of autoregressive models
by partial autocorrelations." Journal of Multivariate Analysis, 3, 408-419

McLeod A. I., Zhang Y (2006).  "Partial autocorrelation parameterization for subset autore-
gression." Journal of Time Series Analysis, 27(4), 599-612
}

\author{
A. I. McLeod
}

\keyword{ ts }
