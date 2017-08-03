#' Exact log-likelihood of a long memory model
#' 
#' Computes the exact log-likelihood of a long memory model with respect to a
#' given time series.
#' 
#' The log-likelihood is computed for the given series z and the parameters.
#' If two or more of \code{dfrac}, \code{H} or \code{alpha} are present and/or
#' two or more of \code{dfs}, \code{Hs} or \code{alphas} are present, an error
#' will be thrown, as otherwise there is redundancy in the model.  Note that
#' non-seasonal and seasonal components can be of different types: for example,
#' there can be seasonal FGN with FDWN at the non-seasonal level.
#' 
#' The moving average parameters are in the Box-Jenkins convention: they are
#' the negative of the parameters given by \code{\link{arima}}.
#' 
#' For the useC parameter, a "0" means no C is used; a "1" means C is only used
#' to compute the log-likelihood, but not the theoretical autocovariance
#' function (tacvf); a "2" means that C is used to compute the tacvf and not
#' the log-likelihood; and a "3" means C is used to compute everything.
#' 
#' Note that the time series is assumed to be stationary: this function does
#' not do any differencing.
#' 
#' @param z A vector or (univariate) time series object, assumed to be (weakly)
#' stationary.
#' @param phi The autoregressive parameters in vector form.
#' @param theta The moving average parameters in vector form.  See Details for
#' differences from \code{\link{arima}}.
#' @param dfrac The fractional differencing parameter.
#' @param phiseas The seasonal autoregressive parameters in vector form.
#' @param thetaseas The seasonal moving average parameters in vector form.  See
#' Details for differences from \code{\link{arima}}.
#' @param dfs The seasonal fractional differencing parameter.
#' @param H The Hurst parameter for fractional Gaussian noise (FGN).  Should
#' not be mixed with \code{dfrac} or \code{alpha}: see "Details".
#' @param Hs The Hurst parameter for seasonal fractional Gaussian noise (FGN).
#' Should not be mixed with \code{dfs} or \code{alphas}: see "Details".
#' @param alpha The decay parameter for power-law autocovariance (PLA) noise.
#' Should not be mixed with \code{dfrac} or \code{H}: see "Details".
#' @param alphas The decay parameter for seasonal power-law autocovariance
#' (PLA) noise.  Should not be mixed with \code{dfs} or \code{Hs}: see
#' "Details".
#' @param period The periodicity of the seasonal components.  Must be >= 2.
#' @param useC How much interfaced C code to use: an integer between 0 and 3.
#' The value 3 is strongly recommended. See "Details".
#' @return The exact log-likelihood of the model given with respect to z, up to
#' an additive constant.
#' @author Justin Veenstra
#' @seealso \code{\link{arfima}}
#' 
#' \code{\link{lARFIMAwTF}}
#' 
#' \code{\link{tacvfARFIMA}}
#' @references Box, G. E. P., Jenkins, G. M., and Reinsel, G. C. (2008) Time
#' Series Analysis: Forecasting and Control.  4th Edition. John Wiley and Sons,
#' Inc., New Jersey.
#' 
#' Veenstra, J.Q. Persistence and Antipersistence:  Theory and
#' Software (PhD Thesis)
#' @keywords ts
#' @examples
#' 
#' set.seed(3452)
#' sim <- arfima.sim(1000, model = list(phi = c(0.3, -0.1)))
#' lARFIMA(sim, phi = c(0.3, -0.1))
#' 
#' @export lARFIMA
"lARFIMA" <- function(z, phi = numeric(0), theta = numeric(0), dfrac = numeric(0), phiseas = numeric(0), 
    thetaseas = numeric(0), dfs = numeric(0), H = numeric(0), Hs = numeric(0), alpha = numeric(0), 
    alphas = numeric(0), period = 0, useC = 3) {
    
    if (is.null(phi) || any(is.na(phi))) 
        phi <- numeric(0)
    if (is.null(theta) || any(is.na(theta))) 
        theta <- numeric(0)
    if (is.null(phiseas) || any(is.na(phiseas))) 
        phiseas <- numeric(0)
    if (is.null(thetaseas) || any(is.na(thetaseas))) 
        thetaseas <- numeric(0)
    
    if (round(useC) != useC) 
        stop("non-integer useC")
    n <- length(z)
    if (useC == 0) 
        useC <- useCt <- F else if (useC == 1) {
        useC <- T
        useCt <- F
    } else if (useC == 2) {
        useC <- F
        useCt <- T
    } else if (useC == 3) 
        useC <- useCt <- T else stop("invalid useC")
    r <- as.double(tacvfARFIMA(phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, 
        thetaseas = thetaseas, dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, 
        period = period, maxlag = (n - 1), useCt = useCt))
    
    if (useC) {
        logl <- tryCatch(DLLoglikelihood(r, z), error = function(err) NULL)
        if (is.null(logl)) 
            logl <- -1e+10
    } else {
        error <- numeric(n)
        sigmasq <- numeric(n)
        error[1] <- z[1]
        sigmasq[1] <- r[1]
        phee <- as.vector(r[2]/r[1])
        error[2] <- z[2] - phee * z[1]
        sigmasqkm1 <- r[1] * (1 - phee^2)
        sigmasq[2] <- sigmasqkm1
        g <- r[1] * sigmasqkm1
        
        for (k in 2:(n - 1)) {
            phikk <- (r[k + 1] - phee %*% rev(r[2:k]))/sigmasqkm1
            sigmasqk <- sigmasqkm1 * (1 - phikk^2)
            phinew <- phee - as.vector(phikk) * rev(phee)
            phee <- c(phinew, phikk)
            sigmasqkm1 <- sigmasqk
            g <- g * sigmasqk
            error[k + 1] <- z[k + 1] - crossprod(phee, rev(z[1:k]))
            sigmasq[k + 1] <- sigmasqk
        }
        S <- sum((error * error)/sigmasq)
        logl <- -n/2 * log(S/n) - log(g)/2
    }
    
    return(logl)
} 
