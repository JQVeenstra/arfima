#' Converts AR/MA coefficients from operator space to the PACF space
#' 
#' Converts AR/MA coefficients from operator space to the PACF box-space;
#' usually for internal use
#'
#'
#' @param phi The AR/MA coefficients in operator space
#' @return The AR/MA coefficients in the PACF space
#' @author A. I. McLeod
#' @references Barndorff-Nielsen O. E., Schou G. (1973). "On the
#' parametrization of autoregressive models by partial autocorrelations."
#' Journal of Multivariate Analysis, 3, 408-419
#'
#' McLeod A. I., Zhang Y (2006).  "Partial autocorrelation parameterization for
#' subset autore- gression." Journal of Time Series Analysis, 27(4), 599-612
#' @keywords ts
#' @export ARToPacf
"ARToPacf" <- function(phi) {
    phik <- phi
    L <- length(phi)
    if (L == 0)
        return(numeric(0))
    pi <- numeric(L)
    for (k in 1:L) {
        LL <- L + 1 - k
        pi[L + 1 - k] <- a <- phik[LL]
        phikp1 <- phik[-LL]
        if (abs(a) == 1)
            break
        phik <- (phikp1 + a * rev(phikp1))/(1 - a^2)
    }
    pi
}



#' Converts AR/MA coefficients from the PACF space to operator space
#'
#' Converts AR/MA coefficients from PACF box-space to operator space; usually
#' for internal use
#'
#'
#' @param pi The AR/MA coefficients in PACF box-space
#' @return The AR/MA coefficients in operator space.
#' @author A. I. McLeod
#' @references Barndorff-Nielsen O. E., Schou G. (1973). "On the
#' parametrization of autoregressive models by partial autocorrelations."
#' Journal of Multivariate Analysis, 3, 408-419
#'
#' McLeod A. I. , Zhang Y (2006).  "Partial autocorrelation parameterization
#' for subset autore- gression." Journal of Time Series Analysis, 27(4),
#' 599-612
#' @keywords ts
#' @export PacfToAR
"PacfToAR" <- function(pi) {
    L <- length(pi)
    if (L == 0)
        return(numeric(0))
    if (L == 1)
        return(pi)
    phik <- pi[1]
    for (k in 2:L) {
        phikm1 <- phik
        phik <- c(phikm1 - pi[k] * rev(phikm1), pi[k])
    }
    phik
}
