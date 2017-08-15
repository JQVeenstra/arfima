arfima.sim.work <- function(n, phi = numeric(0), theta = numeric(0), dfrac = numeric(0),
    H = numeric(0), alpha = numeric(0), phiseas = numeric(0), thetaseas = numeric(0), dfs = numeric(0),
    Hs = numeric(0), alphas = numeric(0), period = 0, useC = T, useCt = T, sigma2 = 1, rand.gen = rnorm,
    innov = NULL, ...) {

    r <- tacvfARFIMA(phi = phi, theta = theta, dfrac = dfrac, H = H, alpha = alpha, phiseas = phiseas,
        thetaseas = thetaseas, dfs = dfs, Hs = Hs, alphas = alphas, period = period, maxlag = n -
            1, useCt = useCt, sigma2 = sigma2)

    if (length(innov) == 0)
        z <- DLSimulate(n, r, rand.gen = rand.gen, ...) else z <- DLSimForNonPar(n, innov, r, dint = 0, dseas = 0, period = 0, muHat = 0, zinit = NULL)  #set to 0 and null because arfima.sim already adds mean and integrates.

    return(z)
}

DLSimForNonPar <- function(n, a, r, dint = 0, dseas = 0, period = 0, muHat = 0, zinit = NULL) {


    nr <- length(r)
    if (min(n, nr) <= 0)
        stop("input error")
    if (nr < n)
        r <- c(r, rep(0, n - nr))
    r <- r[1:n]
    EPS <- .Machine$double.eps  # 1+EPS==1, machine epsilon
    z <- numeric(n)
    out <- .C("durlevsim", z = as.double(z), as.double(a), as.integer(n), as.double(r),
        as.double(EPS), fault = as.integer(1))
    fault <- out$fault
    if (fault == 1)
        stop("error: sequence not p.d.")
    z <- out$z + muHat

    z
}




#' Simulate an ARFIMA time series.
#'
#' This function simulates an long memory ARIMA time series, with one of
#' fractionally differenced white noise (FDWN), fractional Gaussian noise
#' (FGN), power-law autocovariance (PLA) noise, or short memory noise and
#' possibly seasonal effects.
#'
#' A suitably defined stationary series is generated, and if either of the
#' dints (non-seasonal or seasonal) are greater than zero, the series is
#' integrated (inverse-differenced) with zinit equalling a suitable amount of
#' 0s if not supplied.  Then a suitable amount of points are taken out of the
#' beginning of the series (i.e. dint + period * seasonal dint = the length of
#' zinit) to obtain a series of length n.  The stationary series is generated
#' by calculating the theoretical autovariance function and using it, along
#' with the innovations to generate a series as in McLeod et. al. (2007).
#'
#' @param n The number of points to be generated.
#' @param model The model to be simulated from.  The phi and theta arguments
#' should be vectors with the values of the AR and MA parameters. Note that
#' Box-Jenkins notation is used for the MA parameters: see the "Details"
#' section of \code{\link{arfima}}.  The dint argument indicates how much
#' differencing should be required to make the process stationary.  The dfrac,
#' H, and alpha arguments are FDWN, FGN and PLA values respectively; note that
#' only one (or none) of these can have a value, or an error is returned. The
#' seasonal argument is a list, with the same parameters, and a period, as the
#' model argument. Note that with a seasonal model, we can have mixing of
#' FDWN/FGN/HD noise: one in the non-seasonal part, and the other in the
#' seasonal part.
#' @param useC How much interfaced C code to use: an integer between 0 and 3.
#' The value 3 is strongly recommended. See the "Details" section of
#' \code{\link{arfima}}.
#' @param sigma2 The desired variance for the innovations of the series.
#' @param rand.gen The distribution of the innovations.  Any distribution
#' recognized by \code{R} is possible
#' @param muHat The theoretical mean of the series before integration (if
#' integer integration is done)
#' @param zinit Used for prediction; not meant to be used directly.  This
#' allows a start of a time series to be specified before inverse differencing
#' (integration) is applied.
#' @param innov Used for prediction; not meant to be used directly.  This
#' allows for the use of given innovations instead of ones provided by
#' \code{rand.gen}.
#' @param \dots Other parameters passed to the random variate generator;
#' currently not used.
#' @return A sample from a multivariate normal distribution that has a
#' covariance structure defined by the autocovariances generated for given
#' parameters.  The sample acts like a time series with the given parameters.
#' @author JQ (Justin) Veenstra
#' @seealso \code{\link{arfima}}
#' @references McLeod, A. I., Yu, H. and Krougly, Z. L. (2007) Algorithms for
#' Linear Time Series Analysis: With R Package Journal of Statistical Software,
#' Vol. 23, Issue 5
#'
#' Veenstra, J.Q. Persistence and Antipersistence:  Theory and
#' Software (PhD Thesis)
#'
#' P. Borwein (1995) An efficient algorithm for Riemann Zeta function Canadian
#' Math. Soc. Conf. Proc., 27, pp. 29-34.
#' @keywords fit ts
#' @examples
#'
#' set.seed(6533)
#' sim <- arfima.sim(1000, model = list(phi = .2, dfrac = .3, dint = 2))
#'
#' fit <- arfima(sim, order = c(1, 2, 0))
#' fit
#'
#' @export arfima.sim
arfima.sim <- function(n, model = list(phi = numeric(0), theta = numeric(0), dint = 0, dfrac = numeric(0),
    H = numeric(0), alpha = numeric(0), seasonal = list(phi = numeric(0), theta = numeric(0),
        dint = 0, period = numeric(0), dfrac = numeric(0), H = numeric(0), alpha = numeric(0))),
    useC = 3, sigma2 = 1, rand.gen = rnorm, muHat = 0, zinit = NULL, innov = NULL, ...) {

    H <- model$H
    dfrac <- model$dfrac
    dint <- model$dint
    alpha <- model$alpha

    if (length(H) == 0)
        H <- numeric(0)
    if (length(dfrac) == 0)
        dfrac <- numeric(0)
    if (length(alpha) == 0)
        alpha <- numeric(0)
    if (length(H) + length(dfrac) + length(alpha) > 1) {
        warning("two or more of dfrac, H and alpha have been specified: using dfrac if it exists, otherwise H")
      if(length(dfrac)==0)  {
        alpha <- numeric(0)
      }
      else {
        H <- alpha <- numeric(0)
      }

    }
    if (useC == 0)
        useC <- useCt <- F else if (useC == 1) {
        useC <- T
        useCt <- F
    } else if (useC == 2) {
        useC <- F
        useCt <- T
    } else useC <- useCt <- T
    if (length(dint) == 0)
        dint <- 0
    if (!is.null(dint) && (round(dint) != dint || dint < 0))
        stop("dint must be an integer >= 0")

    phi <- model$phi
    theta <- model$theta
    seasonal <- model$seasonal
    if (length(phi) == 0)
        phi <- numeric(0)
    if (length(theta) == 0)
        theta <- numeric(0)
    if (!is.null(seasonal)) {
        if ((length(seasonal$period) == 0) || (seasonal$period == 0)) {
            dseas <- 0
            period <- 0
            dfs <- numeric(0)
            phiseas <- numeric(0)
            thetaseas <- numeric(0)
            Hs <- numeric(0)
            alphas <- numeric(0)
        } else {
            dseas <- seasonal$dint
            period <- seasonal$period
            dfs <- seasonal$dfrac
            Hs <- seasonal$H
            alphas <- seasonal$alpha

            if (length(period) == 0 || period < 2 || period != round(period))
                stop("if using a periodic model, must have an integer period >=2")
            if (length(dseas) == 0)
                dseas <- 0
            if (round(dseas) != dseas || dseas < 0)
                stop("seasonal d must be an integer >= 0")
            if (length(dfs) == 0)
                dfs <- numeric(0)
            if (length(Hs) == 0)
                Hs <- numeric(0)
            if (length(alphas) == 0)
                alphas <- numeric(0)
            if (length(Hs) + length(dfs) + length(alphas) > 1) {
                warning("two or more of seasonal dfrac, H and alpha have been specified: using only seasonal dfrac")
                Hs <- numeric(0)
            }
            phiseas <- seasonal$phi
            thetaseas <- seasonal$theta
            if (length(phiseas) == 0)
                phiseas <- numeric(0)
            if (length(thetaseas) == 0)
                thetaseas <- numeric(0)
        }

    } else {
        seasonal <- NULL
        dseas <- 0
        period <- 0
        dfs <- numeric(0)
        phiseas <- numeric(0)
        thetaseas <- numeric(0)
        Hs <- numeric(0)
        alphas <- numeric(0)
    }

    if (!IdentInvertQ(phi = phi, theta = theta, dfrac = dfrac, H = H, alpha = alpha, phiseas = phiseas,
        thetaseas = thetaseas, dfs = dfs, Hs = Hs, alphas = alphas, period = period, ident = TRUE)) {
        stop("Model is non-causal, non-invertible or non-identifiable")
    }


    z <- arfima.sim.work(n, phi = phi, theta = theta, dfrac = dfrac, H = H, alpha = alpha,
        phiseas = phiseas, thetaseas = thetaseas, dfs = dfs, Hs = Hs, alphas = alphas, period = period,
        useC = useC, useCt = useCt, sigma2 = sigma2, rand.gen = rand.gen, innov = innov,
        ...)
    z <- z - mean(z) + muHat
    if (dint + dseas > 0) {
        icap <- dint + dseas * period
        if (is.null(zinit))
            zinit <- rep(0, icap)
        z <- integ(z, zinit, dint, dseas, period)[-c(1:icap)]
    }
    as.ts(z)
}
