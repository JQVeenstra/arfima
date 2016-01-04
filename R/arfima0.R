#' Exact MLE for ARFIMA
#' 
#' The time series is corrected for the sample mean and then exact MLE is used
#' for the other parameters. This is a simplified version of the arfima()
#' function that may be useful in simulations and bootstrapping.
#' 
#' The sample mean is asymptotically efficient.
#' 
#' @param z time series
#' @param order (p,d,q) where p=order AR, d=regular difference, q=order MA
#' @param lmodel type of long-memory component: FD, FGN, PLA or NONE
#' @return list with components: \item{bHat}{transformed optimal parameters}
#' \item{alphaHat}{estimate of alpha} \item{HHat}{estimate of H}
#' \item{dHat}{estimate of d} \item{phiHat}{estimate of phi}
#' \item{thetaHat}{estimate of theta} \item{wLL}{optimized value of Whittle
#' approximate log-likelihood} \item{LL}{corresponding exact log-likelihood}
#' \item{convergence}{convergence indicator}
#' @author JQ (Justin) Veenstra and A. I. McLeod
#' @keywords ts
#' @examples
#' 
#' z <- rnorm(100)
#' arfima0(z, lmodel="FGN")
#' 
#' @export arfima0
arfima0 <- function(z, order = c(0, 0, 0), lmodel = c("FD", "FGN", "PLA", "NONE")) {
    lmodel <- match.arg(lmodel)
    p <- order[1]
    d <- order[2]
    q <- order[3]
    stopifnot(p >= 0, q >= 0, d >= 0)
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
    stopifnot(is.wholenumber(p), is.wholenumber(d), is.wholenumber(q))
    w <- if (d > 0) 
        diff(z, differences = d) else z
    w <- w - mean(w)
    n <- length(w)
    alg <- 1
    binit <- numeric(p + q + 1)
    binit[1] <- 1
    penaltyLoglikelihood <- (-n/2 * log(sum(w^2)/n)) * 0.01
    Entropy <- function(beta, p, q) {
        alpha <- beta[1]
        phi <- theta <- numeric(0)
        if (alpha <= 0 || alpha >= 2 || (p > 0 || q > 0) && abs(beta[-1]) >= 1) 
            LL <- -penaltyLoglikelihood else {
            if (p > 0) 
                phi <- PacfToAR(beta[2:(p + 1)])
            if (q > 0) 
                theta <- PacfToAR(beta[(p + 2):(p + q + 1)])
            r <- switch(lmodel, FD = tacvfARFIMA(phi = phi, theta = theta, dfrac = 0.5 - 
                alpha/2, maxlag = n - 1), FGN = tacvfARFIMA(phi = phi, theta = theta, H = 1 - 
                alpha/2, maxlag = n - 1), PLA = tacvfARFIMA(phi = phi, theta = theta, alpha = alpha, 
                maxlag = n - 1), NONE = tacvfARFIMA(phi = phi, theta = theta, maxlag = n - 
                1))
            LL <- -DLLoglikelihood(r, w)
        }
        LL
    }
    if (p + q > 0 || lmodel != "NONE") {
        ans <- optim(par = binit, fn = Entropy, p = p, q = q, method = "L-BFGS-B", lower = c(0.01, 
            rep(-0.99, p + q)), upper = c(1.99, rep(0.99, p + q)), control = list(trace = 0))
        if (ans$convergence != 0) {
            # convergence problem. Use Nelder-Mead with penalty function
            alg <- 2
            ans <- optim(par = binit, fn = Entropy, method = "Nelder-Mead")
        }
        bHat <- ans$par
        LL <- -ans$value
        convergence <- ans$convergence
    } else {
        bHat <- numeric(0)
        LL <- -Entropy(1, 0, 0)
        convergence <- 0
    }
    alphaHat <- bHat[1]
    HHat <- 1 - alphaHat/2
    dHat <- HHat - 0.5
    phiHat <- thetaHat <- numeric(0)
    if (p > 0) 
        phiHat <- PacfToAR(bHat[2:(p + 1)])
    if (q > 0) 
        thetaHat <- PacfToAR(bHat[(p + 2):(p + q + 1)])
    ans <- list(bHat = bHat, alphaHat = alphaHat, HHat = HHat, dHat = dHat, phiHat = phiHat, 
        thetaHat = thetaHat, LL = LL, convergence = convergence)
    unlist(ans)
} 
