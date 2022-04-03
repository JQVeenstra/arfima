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
#' \emph{Note:} if you would like to fit a function from a fitted arfima model,
#' the function \code{sim_from_fitted} can be used.
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
#' @seealso \code{\link{arfima}}, \code{\link{sim_from_fitted}}
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

    if(class(model)=='arfima') {
      warning('Model was of class arfima, using the function sim_from_fitted')
      return(sim_from_fitted(n, model))
    }
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

#' Simulate an ARFIMA time series from a fitted arfima object.
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
#' \emph{Note:} if you would like to fit from parameters, use the funtion,
#'  \code{arfima.sim}.
#'
#' @param n The number of points to be generated.
#' @param model The model to be simulated from.  The phi and theta arguments
#' should be vectors with the values of the AR and MA parameters. Note that
#' Box-Jenkins notation is used for the MA parameters: see the "Details"
#' section of \code{\link{arfima}}.
#' @param X The xreg matrix to add to the series, \emph{required} if there is an xreg
#' argument in \code{model}.  An error will be thrown if there is a mismatch
#' between this argument and whether \code{model} was called with a external
#' regressor
#' @param seed An optional seed that will be set before the simulation.  If 
#' \code{model} is multimodal, a seed will be chosen randomly if not provided, 
#' and all modes will simulate a time series with said seed set.
#' @return A sample (or list of samples) from a multivariate normal distribution that has a
#' covariance structure defined by the autocovariances generated for given
#' parameters.  The sample acts like a time series with the given parameters.
#' The returned value will be a list if the fit is multimodal.
#' @author JQ (Justin) Veenstra
#' @seealso \code{\link{arfima}}, \code{\link{arfima.sim}}
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
#' sim2 <- sim_from_fitted(100, fit)
#' 
#' fit2 <- arfima(sim2, order = c(1, 2, 0))
#' fit2
#' 
#' set.seed(2266)
#' #Fairly pathological series to fit for this package
#' series = arfima.sim(500, model=list(phi = 0.98, dfrac = 0.46))
#' 
#' X = matrix(rnorm(1000), ncol = 2)
#' colnames(X) <- c('c1', 'c2')
#' series_added <- series + X%*%c(2, 5)
#' 
#' fit <- arfima(series, order = c(1, 0, 0), numeach = c(2, 2))
#' fit_X <- arfima(series_added, order=c(1, 0, 0), xreg=X, numeach = c(2, 2))
#' 
#' from_series <- sim_from_fitted(1000, fit)
#'  
#' fit1a <- arfima(from_series[[1]], order = c(1, 0, 0), numeach = c(2, 2))
#' fit1a
#' fit1 <- arfima(from_series[[1]], order = c(1, 0, 0))
#' fit1
#' fit2 <- arfima(from_series[[1]], order = c(1, 0, 0))
#' fit2
#' fit3 <- arfima(from_series[[1]], order = c(1, 0, 0))
#' fit3
#' fit4 <- arfima(from_series[[1]], order = c(1, 0, 0))
#' fit4
#' 
#' Xnew = matrix(rnorm(2000), ncol = 2)
#' from_series_X <- sim_from_fitted(1000, fit_X, X=Xnew)
#' 
#' fit_X1a <- arfima(from_series_X[[1]], order=c(1, 0, 0), xreg=Xnew, numeach = c(2, 2))
#' fit_X1a
#' fit_X1 <- arfima(from_series_X[[1]], order=c(1, 0, 0), xreg=Xnew)
#' fit_X1
#' fit_X2 <- arfima(from_series_X[[2]], order=c(1, 0, 0), xreg=Xnew)
#' fit_X2
#' fit_X3 <- arfima(from_series_X[[3]], order=c(1, 0, 0), xreg=Xnew)
#' fit_X3
#' fit_X4 <- arfima(from_series_X[[4]], order=c(1, 0, 0), xreg=Xnew)
#' fit_X4
#' 
#' @export sim_from_fitted
sim_from_fitted <- function(n, model, X = NULL, seed = NULL) {
  
  if(class(model)!='arfima') {
    stop("Cannot use this function with any other input than a fitted arfima function")
  }
  
  if(any(model$r > 0) || any(model$s>1)) {
    stop('Cannot currently simulate from a transfer function model fit.')
  }
  
  if(!is.null(model$xreg) && is.null(X))
    stop('Cannot simulate from a model with xreg and no X')
  
  if(is.null(model$xreg) && !is.null(X))
    stop('Cannot simulate from a model with X and no xreg')
  
  if(model$regOnly) 
    warning('This is a pure regression model.')
  
  
  
  if(!is.null(X)) {
    d <- dim(X)
    if (is.null(d)) {
      X <- as.matrix(X)
      d <- dim(X)
    }
    if(d[1]!=n)
      stop('X does not have n rows')
    if(d[2]!=dim(model$xreg)[2])
      stop('X does not match xreg size')
    
  }
  mmodal<- FALSE
  if(length(model$modes)>1) {
    warning('Model has more than one mode:  what is returned will be a list containing all simulations.')
    mmodal <- TRUE
    if(is.null(seed)) {
      seed <- sample(999999, 1)
      p <- paste0('Since no seed is set, and since multiple simluations are being returned, the seed ', seed)
      p <- paste0(p, '\n has been set for you.  If you prefer different behaviour, let me know on GitHub.')
      warning(p)
    }
  }
  
  
  
  
  
  if(mmodal)
    ret <- list()
  period <- model$period
  dint <- model$dint
  dseas <- model$dseas
  for (i in 1:length(model$modes)) {
    phi <- model$modes[[i]]$phi
    theta <- model$modes[[i]]$theta
    phiseas <- model$modes[[i]]$phiseas
    thetaseas <- model$modes[[i]]$thetaseas
    
    dfrac <- model$modes[[i]]$dfrac
    dfs <- model$modes[[i]]$dfs
    H <- model$modes[[i]]$H
    Hs <- model$modes[[i]]$Hs
    alpha <- model$modes[[i]]$alpha
    alphas <- model$modes[[i]]$alphas
    
    omega <- model$modes[[i]]$omega
    muHat <- model$modes[[i]]$muHat
    sigma2 <- model$modes[[i]]$sigma2
    if(length(muHat)==0) muHat <- 0
    
    if(model$regOnly) {
      series <- rnorm(n)
      series <- series + X %*% omega + muHat
    }
    else {
      series <- arfima.sim(n, model = list(phi=phi, theta=theta, dint = dint, dfrac = dfrac, H = H, alpha = alpha, 
                                           seasonal = list(period=period, phi=phiseas, theta=thetaseas,
                                                           dint = dseas, dfrac = dfs, H = Hs, alpha = alphas)),
                           muHat = muHat, sigma2 = sigma2)
      if(!is.null(model$xreg)) {
        series <- series + X%*%omega
      }
    }
    if(!is.null(seed)) set.seed(seed)
    if(mmodal)
      ret[[i]] <- as.ts(series)
    else
      ret <- as.ts(series)
  }
    
  ret
}
