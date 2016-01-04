#' A parametric bootstrap based on an arfima fit.
#'
#' Simulates a parametric bootstrap from a fitted model using either standard
#' normal innovations or residuals from the fit.
#'
#'
#' @param obj An object of class "arfima".
#' @param R The number of bootstrap replicates to be generated (for each mode).
#' @param pred If TRUE, use a bootstrap sample of the residuals instead of
#' standard normal innovations. This is for predictions.
#' @param seed The seed to use when computing the bootstrap replicates.
#' @param \dots Optional arguments. Not currently used.
#' @return A simulated time series with the same length as the original fitted
#' time series is produced (for each mode) when R=1. When R>1, a matrix with R
#' columns is produced with each column a separate bootstrap realization (for
#' each mode).  This is usually used by predict.
#' @author JQ (Justin) Veenstra and A. I. McLeod
#' @seealso \code{\link{Boot}}
#' @references Veenstra, J.Q. Persistence and Antipersistence:  Theory and
#' Software (PhD Thesis)
#' @keywords ts
Boot.arfima <- function(obj, R = 1, pred = FALSE, seed = NA, ...) {

    if (length(seed) > 1 && length(seed) != R) {
        warning("set seed length is not equal to R: using only first seed")
        seed <- seed[1]
    }

    n <- length(obj$z)

    m <- length(obj$modes)
    zinit <- NULL
    dint <- obj$dint
    dseas <- obj$dseas
    res <- vector("list", m)

    for (i in 1:m) {

        cmode <- obj$modes[[i]]
        z <- Boot.ARFIMA(cmode)
        res[[i]] <- z
    }
    names(res) <- 1:m
    res
}



Boot.ARFIMA <- function(obj, dint, dseas, period, R = 1, n.ahead = NULL, n = 1, zinit = NULL,
    lastpoint = 0, pred = FALSE, seed = NA, ...) {
    if (length(seed) > 1 && length(seed) != R) {
        warning("set seed length is not equal to R: using only first seed")
        seed <- seed[1]
    }
    if (is.null(n.ahead))
        n.ahead <- n

    phi <- obj$phi
    theta <- obj$theta
    phiseas <- obj$phiseas
    thetaseas <- obj$thetaseas
    dfrac <- obj$dfrac
    dfs <- obj$dfs
    logl <- obj$logl
    sigma2 <- obj$sigma2
    muHat <- obj$muHat
    H <- obj$H
    Hs <- obj$Hs
    H <- obj$H
    Hs <- obj$Hs

    if (R == 1) {
        if (length(seed) == 1 && !is.na(seed))
            set.seed(seed)
        if (pred) {
            indexes <- sample(1:n, replace = TRUE)
            innov <- obj$residuals[indexes]
        } else innov <- NULL
        if (!is.null(innov))
            sigma2 <- 1
        z <- arfima.sim(n.ahead, model = list(phi = phi, theta = theta, dint = dint, dfrac = dfrac,
            H = H, seasonal = list(phi = phiseas, theta = thetaseas, dint = dseas, period = period,
                dfrac = dfs, H = Hs)), sigma2 = sigma2, useC = 3, zinit = zinit, muHat = muHat,
            innov = innov)
        if(dint>0) z <- z - muHat

    } else {
        z <- matrix(0, nrow = n.ahead, ncol = R)
        if (length(seed) == 1 && !is.na(seed))
            set.seed(seed)
        for (i in 1:R) {
            if (length(seed) == R)
                set.seed(seed[i])
            if (pred) {
                indexes <- sample(1:n, replace = TRUE)
                innov <- obj$residuals[indexes]
            } else innov <- NULL
            if (!is.null(innov))
                sigma2 <- 1
            zz <- arfima.sim(n.ahead, model = list(phi = phi, theta = theta, dint = dint,
                dfrac = dfrac, H = H, seasonal = list(phi = phiseas, theta = thetaseas,
                  dint = dint, period = period, dfrac = dfs, H = Hs)), sigma2 = sigma2,
                useC = 3, zinit = zinit, muHat = muHat, innov = innov)

            z[, i] <- zz
        }
    }

    z
}
