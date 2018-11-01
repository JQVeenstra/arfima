InvertibleQ <- function(phi) {
    if (length(phi) == 0) 
        return(T)
    return(all(abs(ARToPacf(phi)) < 1))
}

InvertibleD <- function(d) {
    if (length(d) == 0) 
        return(T)
    return(d > -1 && d < 0.5)
}

InvertibleH <- function(H) {
    if (length(H) == 0) 
        return(T)
    return(H > 0 && H < 1)
}

InvertibleAlpha <- function(alpha) {
    if (length(alpha) == 0) 
        return(T)
    return(alpha > 0 && alpha < 3)
}

IdentifiableQ <- function(phi = numeric(0), theta = numeric(0)) {
    p <- length(phi)
    q <- length(theta)
    k <- p + q
    A <- matrix(0, nrow = k, ncol = k)
    if (p > 0) 
        for (j in 1:p) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i, j] <- -1
                if (i - j <= q && i - j > 0) 
                  A[i, j] <- theta[i - j]
            }
        }
    if (q > 0) 
        for (j in 1:q) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i, j + p] <- 1
                
                if (i - j <= p && i - j > 0) 
                  A[i, j + p] <- -phi[i - j]
            }
        }
    return(det(A) > 0)
}

IdentifiableQQ <- function(phi = numeric(0), theta = numeric(0), phiseas = numeric(0), thetaseas = numeric(0), 
    negphi = T) {
    p <- length(phi)
    q <- length(theta)
    ps <- length(phiseas)
    qs <- length(thetaseas)
    
    k <- p + ps + q + qs
    A <- matrix(0, nrow = k, ncol = k)
    if (p > 0) 
        for (j in 1:p) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i, j] <- -1
                if (i - j <= q && i - j > 0) 
                  A[i, j] <- theta[i - j]
            }
        }
    if (q > 0) 
        for (j in 1:q) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i, j + p] <- if (negphi) 
                    1 else -1
                if (i - j <= p && i - j > 0) 
                  A[i, j + p] <- if (negphi) 
                    -phi[i - j] else phi[i - j]
            }
        }
    if (ps > 0) 
        for (j in 1:ps) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i + p + q, j + p + q] <- -1
                if (i - j <= qs && i - j > 0) 
                  A[i + p + q, j + p + q] <- thetaseas[i - j]
            }
        }
    if (qs > 0) 
        for (j in 1:qs) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i + p + q, j + p + q + ps] <- if (negphi) 
                    1 else -1
                if (i - j <= ps && i - j > 0) 
                  A[i + p + q, j + p + q + ps] <- if (negphi) 
                    -phiseas[i - j] else phiseas[i - j]
            }
        }
    return(det(A) > 0)
}



#' Checks invertibility, stationarity, and identifiability of a given set of
#' parameters
#' 
#' Computes whether a given long memory model is invertible, stationary, and
#' identifiable.
#' 
#' This function tests for identifiability via the information matrix of the
#' ARFIMA process.  Whether the process is stationary or invertible amounts to
#' checking whether all the variables fall in correct ranges.
#' 
#' The moving average parameters are in the Box-Jenkins convention: they are
#' the negative of the parameters given by \code{\link{arima}}.
#' 
#' If \code{dfrac}/\code{H}/\code{alpha} are mixed and/or
#' \code{dfs}/\code{Hs}/\code{alphas} are mixed, an error will not be thrown,
#' even though only one of these can drive the process at either level. Note
#' also that the FGN or PLA have no impact on the identifiability of the model,
#' as information matrices containing these parameters currently do not have
#' known closed form.  These two parameters must be within their correct ranges
#' (0<H<1 for FGN and 0 < alpha < 3 for PLA.)
#' 
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
#' @param delta The delta parameters for transfer functions.
#' @param period The periodicity of the seasonal components.  Must be >= 2.
#' @param debug When TRUE and model is not stationary/invertible or
#' identifiable, prints some helpful output.
#' @param ident Whether to test for identifiability.
#' @return TRUE if the model is stationary, invertible and identifiable.  FALSE
#' otherwise.
#' @author Justin Veenstra
#' @seealso \code{\link{iARFIMA}}
#' @references McLeod, A.I. (1999) Necessary and sufficient condition for
#' nonsingular Fisher information matrix in ARMA and fractional ARMA models The
#' American Statistician 53, 71-72.
#' 
#' Veenstra, J. and McLeod, A. I. (2012, Submitted) Improved Algorithms for
#' Fitting Long Memory Models: With R Package
#' @keywords ts
#' @examples
#' 
#' IdentInvertQ(phi = 0.3, theta = 0.3)
#' IdentInvertQ(phi = 1.2)
#' 
#' @export IdentInvertQ
IdentInvertQ <- function(phi = numeric(0), theta = numeric(0), phiseas = numeric(0), thetaseas = numeric(0), 
    dfrac = numeric(0), dfs = numeric(0), H = numeric(0), Hs = numeric(0), alpha = numeric(0), 
    alphas = numeric(0), delta = numeric(0), period = 0, debug = FALSE, ident = TRUE) {
    if (!(InvertibleQ(phi) && InvertibleQ(theta) && InvertibleQ(phiseas) && InvertibleQ(thetaseas) && 
        InvertibleD(dfrac) && InvertibleD(dfs) && InvertibleH(H) && InvertibleH(Hs) && InvertibleQ(delta) && 
        InvertibleAlpha(alpha) && InvertibleAlpha(alphas))) {
        if (debug) 
            warning("Model is non-stationary or non-invertible")
        return(FALSE)
    }
    if ((length(H) + length(dfrac) + length(alpha) > 1) || (length(Hs) + length(dfs) + length(alphas) > 
        1)) {
        if (debug) 
            warning("Model is contains fractional d or H in either seasonal or non-seasonal components")
        return(FALSE)
    }
    if (ident && (length(phi) > 0 || length(theta) > 0 || length(phiseas) > 0 || length(thetaseas) > 
        0 || length(dfrac) > 0 || length(dfs) > 0)) {
        if (length(dfrac) > 0) 
            d <- T else d <- F
        if (length(dfs) > 0) 
            dd <- T else dd <- F
        I <- iARFIMA(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            dfrac = d, dfs = dd, period = period)
        if (!recdet(I)) {
            if (debug) 
                warning("Information matrix of SARMA or ARFIMA process is not positive definite")
            return(FALSE)
        }
    }
    
    TRUE
}


recdet = function(a) {
    p <- sqrt(length(a))
    if (p != round(p)) 
        stop("not a square matrix")
    a <- as.matrix(a, nrow = p)
    for (i in 1:p) {
        if (det(as.matrix(a[1:i, 1:i], nrow = i)) <= 0) 
            return(F)
    }
    return(T)
} 
