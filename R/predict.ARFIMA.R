#' Predicts from a fitted object.
#'
#' Performs prediction of a fitted \code{arfima} object. Includes prediction
#' for each mode and exact and limiting prediction error standard deviations.
#' \strong{NOTE:  the standard errors in beta are currently not taken into
#' account in the prediction intervals shown.  This will be updated as soon
#' as possible.}
#'
#'
#' @param object A fitted \code{arfima} object
#' @param n.ahead The number of steps ahead to predict
#' @param prop.use The proportion (between 0 and 1) or percentage (between
#' >1 and 100) of data points to use for prediction.  Defaults to the string
#' "default", which sets the number of data points \code{n.use} to the minimum
#' of the series length and 1000.  Overriden by \code{n.use}.
#' @param n.use Directly set the number mentioned in \code{prop.use}.
#' @param newxreg If a regression fit, the new regressors
#' @param predint The percentile to use for prediction intervals assuming normal
#' deviations.
#' @param exact Controls whether exact (based on the theoretical autocovariance
#' matrix) prediction variances are calculated (which is recommended), as well
#' as whether the exact prediction formula is used when the process is
#' differenced (which can take a fair amount of time if the length of the series
#' used to predict is large).  Defaults to the string "default", which is
#' \code{TRUE} for the first and \code{FALSE} for the second.  A Boolean value
#' (\code{TRUE} or \code{FALSE}) will set both to this value.
#' @param setmuhat0 Experimental. Sets muhat equal to zero
#' @param cpus The number of CPUs to use for prediction. Currently not
#' implemented
#' @param trend An optional vector the length of \code{n.ahead} or longer to
#' add to the predictions
#' @param xreg Alias for newxreg
#' @param \dots Optional arguments. Currently not used
#' @return A list of lists, ceiling(prop.use * n)one for each mode with relavent details about the
#' prediction
#' @author JQ (Justin) Veenstra
#' @seealso \code{\link{arfima}}, \code{\link{plot.predarfima}},
#' \code{\link{print.predarfima}}
#' @references Veenstra, J.Q. Persistence and Antipersistence:  Theory and
#' Software (PhD Thesis)
#' @keywords ts
#' @examples
#'
#' \donttest{
#' set.seed(82365)
#' sim <- arfima.sim(1000, model = list(dfrac = 0.4, theta=0.9, dint = 1))
#' fit <- arfima(sim, order = c(0, 1, 1), back=T)
#' fit
#' pred <- predict(fit, n.ahead = 5)
#' pred
#' plot(pred, numback=50)
#' #Predictions aren't really different due to the
#' #series.  Let's see what happens when we regress!
#'
#' set.seed(23524)
#' #Forecast 5 ahead as before
#' #Note that we need to integrate the regressors, since time series regression
#' #usually assumes that regressors are of the same order as the series.
#' n.fore <- 5
#' X <- matrix(rnorm(3000+3*n.fore), ncol = 3)
#' X <- apply(X, 2, cumsum)
#' Xnew <- X[1001:1005,]
#' X <- X[1:1000,]
#' beta <- matrix(c(2, -.4, 6), ncol = 1)
#' simX <- sim + as.vector(X%*%beta)
#' fitX <- arfima(simX, order = c(0, 1, 1), xreg = X, back=T)
#' fitX
#' #Let's compare predictions.
#' predX <- predict(fitX, n.ahead = n.fore, xreg = Xnew)
#' predX
#' plot(predX, numback = 50)
#' #With the mode we know is really there, it looks better.
#' fitX <- removeMode(fitX, 2)
#' predXnew <- predict(fitX, n.ahead = n.fore, xreg = Xnew)
#' predXnew
#' plot(predXnew, numback=50)
#'
#' #
#' }
#'
predict.arfima <- function(
  object, n.ahead = 1, prop.use = "default", newxreg = NULL, predint = 0.95,
  exact = c("default", T, F), setmuhat0 = FALSE, cpus = 1,
  trend = NULL, n.use = NULL, xreg = NULL,...)
  {
    if(!is.null(xreg)) {
      if(!is.null(newxreg)) stop('Please only use xreg or newxreg!')
      newxreg <- xreg
    }
    n <- object$n
    if(is.null(n.use)) {
      if(prop.use == "default") n.use <- 1000
      else if(is.numeric(prop.use)) {
        if(prop.use > 0 && prop.use <= 1)
          n.use <- ceiling(prop.use * n)
        else if(prop.use > 1 &&  prop.use <= 100)
          n.use <- ceiling(prop.use * n/100)
        else stop('Improper prop.use:  must be between 0 and 1 or 1 and 100')
      }
      else {
        if(is.numeric(n.use)) {
          n.use <- ceiling(n.use)
        }
        else {
          stop(paste('Invalid n.use value of ', n.use))
        }
      }
    }
    stopifnot(n.use>0)
    if(n>=1000 && n.use < 1000)
      warning(paste('For longer time series it is best to use at least\n
                    1000 points in prediction:', n.use, 'being used here.'))

    if(exact[1] == "default") {
      trex <- FALSE
      exact <- !object$differencing
    }
    else if(is.logical(exact) && exact) {
      trex <- exact <- TRUE
    }
    else if(is.logical(exact) && !exact) {
      trex <- exact <- FALSE
    }
    else stop(paste('Got', exact, 'in exact; expecting T, F, or "default"'))

    if (!is.null(object$s) && any(object$s > 1))
        stop("predict only takes static regression parameters, not full transfer functions/dynamic regressions")
    if (!is.null(object$r) && any(object$r > 1))
        stop("predict only takes static regression parameters, not full transfer functions/dynamic regressions")
    myNROW <- function(x) if (is.null(x))
        0 else nrow(x)

    myNCOL <- function(x) if (is.null(x))
        0 else ncol(x)


    nrxreg <- myNROW(newxreg)
    if (nrxreg > 0)
      if(is.null(object$xreg))
        stop("Can't predict with (new)xreg when no xreg in fit")
      else if(!object$strReg)
        stop("Only predict with regular regression at this time (no transfer functions)")

    if (nrxreg > 0 && is.data.frame(newxreg)) {
      namexreg <- setdiff(object$namexreg, "Intercept")
      if(length(union(colnames(newxreg), namexreg)) != length(colnames(newxreg)))
        stop("Named arguments in xreg and (new)xreg are different")
      if(length(union(colnames(newxreg), namexreg)) != length(namexreg))
        stop("Named arguments in xreg and (new)xreg are different")
      newxreg <- newxreg[namexreg]
    }


    if (nrxreg > 0 && is.null(object$xreg))
        stop("no xreg in input to arfima, but (new)xreg input to predict")

    if (nrxreg==0 && !is.null(object$xreg))
        stop("xreg in input to arfima, but no (new)xreg input to predict")
    if (myNCOL(object$xreg) != myNCOL(newxreg))
        stop("unconformable newxreg and xreg")

    if (myNROW(newxreg) && (myNROW(newxreg) != n.ahead))
        stop("unconformable 'newxreg' and 'n.ahead'")


    m <- length(object$modes)
    tacvfs <- tacvf(obj = object, xmaxlag = n.ahead, forPred = TRUE, n.ahead = n.ahead +
        1)
    if (nrxreg > 0) {
        res <- resid(object, reg = TRUE)
        for (i in 1:m) res[[i]] <- as.vector(res[[i]])
        y <- NULL
    } else y <- as.vector(object$z)
    dint <- object$dint
    dseas <- object$dseas
    period <- object$period

    ## multicore!!  call function inside function!
    preds <- vector("list", m)
    limiting <- FALSE
    for (i in 1:m) {
        zy <- if (!is.null(y))
            y else res[[i]]
        nz <- length(zy)
        zz <- c(zy, rep(0, n.ahead))
        muHat <- tacvfs[[i + 1]]$muHat

        if (object$differencing) {
            if (setmuhat0)
                muHat <- 0
            icap <- dint + dseas * period
            zinit <- zy[(nz - icap + 1):nz]
            if (dint > 0)
                zy <- diff(zy, differences = dint)
            if (dseas > 0)
                zy <- diff(zy, differences = dseas, lag = period)
            if (nrxreg) {
              newxreg <- as.matrix(newxreg)
              xreginit <- object$xreg[(nz - icap + 1):nz,]
              newxreg <- rbind(xreginit, newxreg)
              if (dint > 0)
                newxreg <- diff(newxreg, differences = dint)
              if (dseas > 0)
                newxreg <- diff(newxreg, differences = dseas, lag = period)
            }
        } else {
            icap <- 0
            zinit <- NULL
        }

        nz <- length(zy)
        startzy <- if(nz-n.use<=1) 1 else nz-n.use
        znew <- predictWork(y = zy[startzy:nz], r = tacvfs[[i + 1]]$tacvf, dint = dint, dseas = dseas,
            period = period, muHat = muHat, n.ahead = n.ahead, trex = trex, exact=exact)

        znew$Forecast <- znew$Forecast


        if (!is.null(trend)) {
            if (length(trend) < n.ahead)
                stop("length of trend less than n.ahead")
            if (length(trend) != n.ahead)
                warning("trend is not the same length as n.ahead")
            meanwx <- trend[1:n.ahead]
        } else meanwx <- rep(0, n.ahead)


        if (nrxreg) {
            meanwx <- meanwx + as.numeric(newxreg %*% object$modes[[i]]$omega)
        }

        znew$Forecast <- znew$Forecast + meanwx
        if (icap > 0) {
            znew$Forecast <- integ(z = znew$Forecast, zinit = zinit, dint = dint, dseas = dseas,
                period = period)
            znew$Forecast <- znew$Forecast[-c(1:icap)]
        }



        znew$lowernp <- znew$uppernp <- znew$medvalnp <- znew$lowerp <- znew$upperp <- znew$medvalp <- NULL

        if (length(tacvfs[[i + 1]]$psis) > 0) {
            limiting <- TRUE
            sigmas <- cumsum(tacvfs[[i + 1]]$psis^2)[1:n.ahead]
            znew$limitVar <- sigmas * tacvfs[[i + 1]]$sigma2
            znew$limitSD <- sqrt(znew$limitVar)
        } else {
            znew$limitVar <- NULL
            znew$limitSD <- NULL
        }
        znew$sigma2 <- tacvfs[[i + 1]]$sigma2
        preds[[i]] <- znew
    }
    preds$z <- object$z
    preds$limiting <- limiting
    preds$predint <- predint
    preds$name <- deparse(substitute(object))
    preds$m <- m
    class(preds) <- "predarfima"
    preds
}

predictWork <- function(y, r, dint, dseas, period, muHat, n.ahead, exact=TRUE, trex = FALSE) {


    ny <- length(y)
    znew <- TrFore(y, r, muHat, maxLead = n.ahead, getV = exact)
    coeffs <- NULL
    if (dint + dseas > 0) {
        if (trex)
            coeffs <- rev(Z(l = n.ahead, d = dint, ds = dseas, s = period)) else coeffs <- wtsforexact(dint = dint, dseas = dseas, period = period, len = n.ahead)
    }

    znew$Forecast <- znew$Forecasts[1, ]
    znew$Forecasts <- NULL

    if (length(coeffs) > 0) {
        exactsig <- rep(0, n.ahead)
        P <- exactVals(r = r, n = ny, maxLead = n.ahead)
        for (i in 1:n.ahead) {
            exactsig[i] <- sum(sapply(1:i, function(u) sapply(1:i, function(v) coeffs[i -
                u + 1] * coeffs[i - v + 1] * (r[abs(u - v) + 1] - P[u, v]))))
        }
    } else exactsig <- as.numeric(znew$SDForecasts[1, ]^2)

    znew$exactVar <- exactsig
    znew$exactSD <- sqrt(exactsig)

    znew
}


Z <- function(l, d, ds, s) {

    if ((d == 0) && (ds == 0))
        return(numeric(0))
    if (ds > 0 && s == 0)
        stop("No period supplied")
    worker1 <- function(m, value, val) {
        if (m > 0) {
            val[m] <- val[m] + value
            if (d > 0) {
                for (i in 1:d) {
                  val <- worker1(m - i, value * choose(d, i) * (-1)^(i + 1), val)
                }
            }
            if (ds > 0) {
                for (j in 1:ds) {
                  val <- worker1(m - j * s, value * choose(ds, j) * (-1)^(j + 1), val)
                }
            }
            if (d > 0 && ds > 0) {
                for (i in 1:d) {
                  for (j in 1:ds) {
                    val <- worker1(m - i - j * s, value * choose(d, i) * choose(ds, j) *
                      (-1)^(i + j + 1), val)
                  }
                }
            }
        }
        val
    }

    val <- numeric(l)
    val <- worker1(l, 1, val)

    val
}

print.predARFIMA <- function(x, digits = max(6, getOption("digits") - 3), ...) {
    n <- length(x$z)
    forecasts <- x$Forecast
    exactSD <- x$exactSD
    approxSD <- x$approxSD
    n.ahead <- length(forecasts)
    ans <- rbind(forecasts, exactSD, approxSD)
    namer <- c("Forecasts", "Exact SD", if (!is.null(approxSD)) "Approximate SD" else NULL)
    rownames(ans) <- namer
    colnames(ans) <- 1:n.ahead
    ret <- list(`Forecasts and SDs` = ans)
    print(ret, digits = digits, ...)
}


TrFore <- function(z, r, zm, maxLead, getV = TRUE) {
    n <- length(z)

    if (length(r) < (n + maxLead - 1))
        stop("error: length(r) must be >= ", n + maxLead - 1)
    if (!(maxLead > 0))
        stop("maxLead must be > 0")

    zk <- z - zm
    zf <- vf <- matrix(numeric(maxLead), ncol = maxLead)
    gk <- t(matrix(r[n + 1 + outer(1:maxLead, 1:n, "-")], ncol = maxLead, byrow = TRUE))

    GI <- TrenchInverse(toeplitz(r[1:n]))
    gkGI <- crossprod(t(gk), GI)
    zf[1, ] <- zm + gkGI %*% zk
    if (getV) {
        for (j in 1:maxLead) vf[1, j] <- r[1] - sum(gkGI[j, ] * gk[j, ])
    }
    dimnames(zf) <- dimnames(vf) <- list(n, 1:maxLead)
    if (getV)
        ans <- list(Forecasts = zf, SDForecasts = sqrt(vf))
    else
      ans <- list(Forecasts = zf)
    ans
}

exactVals <- function(r, n, maxLead, way = 5) {
    if (length(r) < n + maxLead)
        stop("r is too short")
    GI <- TrenchInverse(toeplitz(r[1:n]))
    vfs <- matrix(numeric(maxLead^2), ncol = maxLead)
    if (way == 1) {
        for (j in 1:(maxLead)) {
            gkj <- matrix(r[(n + j):(j + 1)], nrow = 1)
            for (h in 1:(maxLead)) {
                gkh <- matrix(r[(n + h):(h + 1)], ncol = 1)
                vfs[j, h] <- gkj %*% GI %*% gkh
            }
        }
    } else if (way == 2) {
        for (j in 1:(maxLead)) {
            gkj <- matrix(r[(n + j):(j + 1)], nrow = 1)
            for (h in 1:(maxLead)) {
                gkh <- matrix(r[(n + h):(h + 1)], ncol = 1)
                vfs[j, h] <- vfs[h, j] <- gkj %*% GI %*% gkh
            }
        }
    } else if (way == 3) {
        mat <- mat1 <- NULL
        for (j in 1:(maxLead)) {
            gkj <- matrix(r[(n + j):(j + 1)], nrow = 1)
            mat <- rbind(mat, gkj)
        }
        vfs <- mat %*% GI %*% t(mat)
    } else if (way == 4) {
        mat <- t(n + 1 + outer(1:maxLead, 1:n, "-"))
        gk <- matrix(r[mat], ncol = maxLead)
        vfs <- crossprod(t(crossprod(gk, GI)), gk)
    } else {
        ## seems to be fastest
        mat <- t(n + 1 + outer(1:maxLead, 1:n, "-"))
        gk <- t(matrix(r[mat], ncol = maxLead))
        vfs <- gk %*% GI %*% t(gk)
    }

    vfs
}
