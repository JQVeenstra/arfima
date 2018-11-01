#' Exact log-likelihood of a long memory model with a transfer function model
#' and series included
#'  
#' Computes the exact log-likelihood of a long memory model with respect to a
#' given time series as well as a transfer fucntion model and series. This
#' function is not meant to be used directly.
#' 
#' Once again, this function should not be used externally.
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
#' @param xr The regressors in vector form
#' @param r The order of the delta(s)
#' @param s The order of the omegas(s)
#' @param b The backshifting to be done
#' @param delta Transfer function parameters as in Box, Jenkins, and Reinsel.
#' Corresponds to the "autoregressive" part of the dynamic regression.
#' @param omega Transfer function parameters as in Box, Jenkins, and Reinsel.
#' Corresponds to the "moving average" part of the dynamic regression: note
#' that omega_0 is not restricted to 1.  See "Details" for issues.
#' @param meanval If the mean is to be estimated dynamically, the mean.
#' @return A log-likelihood value
#' @author Justin Veenstra
#' @references Veenstra, J.Q. Persistence and Antipersistence:  Theory and
#' Software (PhD Thesis)
#' @keywords ts
#' @export lARFIMAwTF
"lARFIMAwTF" <- function(z, phi = numeric(0), theta = numeric(0), dfrac = numeric(0), phiseas = numeric(0), 
    thetaseas = numeric(0), dfs = numeric(0), H = numeric(0), Hs = numeric(0), alpha = numeric(0), 
    alphas = numeric(0), xr = numeric(0), r = numeric(0), s = numeric(0), b = numeric(0), 
    delta = numeric(0), omega = numeric(0), period = 0, useC = 3, meanval = 0) {
    
    
    if (is.null(phi) || any(is.na(phi)) || all(phi == 0)) 
        phi <- numeric(0)
    if (is.null(theta) || any(is.na(theta)) || all(theta == 0)) 
        theta <- numeric(0)
    if (is.null(phiseas) || any(is.na(phiseas)) || all(phiseas == 0)) 
        phiseas <- numeric(0)
    if (is.null(thetaseas) || any(is.na(thetaseas)) || all(thetaseas == 0)) 
        thetaseas <- numeric(0)
    if (is.null(delta) || any(is.na(delta))) 
        delta <- numeric(0)
    
    n <- length(z)
    
    
    if (length(delta) != sum(r)) 
        stop("delta wrong length")
    if (length(omega) != sum(s)) 
        stop("omega wrong length")
    if (length(r) != length(s) || length(r) != length(b)) 
        stop("model specifications incorrect.")
    
    z <- funcTF(y = z, x = xr, delta = delta, omega = omega, b = b, rx = r, sx = s, nx = n, 
        meanval = meanval)$y
    
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
    r <- tacvfARFIMA(phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
        dfs = dfs, H = H, alpha = alpha, alphas = alphas, Hs = Hs, period = period, maxlag = (n - 
            1), useCt = useCt)
    if (is.null(r)) 
        return(-1e+10)
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

funcTF <- function(y, x, delta, omega, b, rx, sx, nx, meanval) {
    
    .C("tfcalc", y = as.numeric(y), n = as.integer(nx), x = as.numeric(x), delta = as.numeric(delta), 
        r = as.integer(rx), omega = as.numeric(omega), s = as.integer(sx), b = as.integer(b), 
        num = as.integer(length(sx)), meanval = as.double(meanval))
} 
