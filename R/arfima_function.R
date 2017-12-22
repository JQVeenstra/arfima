
#' Fit ARFIMA, ARIMA-FGN, and ARIMA-PLA (multi-start) models
#'
#' Fits ARFIMA/ARIMA-FGN/ARIMA-PLA multi-start models to times series data.
#' Options include fixing parameters, whether or not to fit fractional noise,
#' what type of fractional noise (fractional Gaussian noise (FGN), fractionally
#' differenced white noise (FDWN), or the newly introduced power-law
#' autocovariance noise (PLA)), etc.  This function can fit regressions with
#' ARFIMA/ARIMA-FGN/ARIMA-PLA errors via the xreg argument, including dynamic
#' regression (transfer functions).
#'
#' Fits by direct optimization using optim.  The optimizer choices are: 0 -
#' BFGS; 1 - Nealder-Mead; 2 - SANN; otherwise CG.
#'
#' A word of warning: it is generally better to use the default, and only use
#' Nelder-Mead to check for spurious modes.  SANN takes a long time (and may
#' only find one mode), and CG may not be stable.
#'
#' If using Nelder-Mead, it must be stressed that Nelder-Mead can take out
#' non-spurious modes or add spurious modes: we have checked visually where we
#' could.  Therefore it is wise to use BFGS as the default and if there are
#' modes close to the boundaries, check using Nelder-Mead.
#'
#' The moving average parameters are in the Box-Jenkins convention: they are
#' the negative of the parameters given by \code{\link{arima}}.  That is, the
#' model to be fit is, in the case of a non-seasonal ARIMA model, phi(B)
#' (1-B)^d z[t] = theta(B) a[t], where phi(B) = 1 - phi(1) B - ... - phi(p) B^p
#' and theta(B) = 1 - theta(1) B - ... - theta(q) B^q.
#'
#' For the useC parameter, a "0" means no C is used; a "1" means C is only used
#' to compute the log-likelihood, but not the theoretical autocovariance
#' function (tacvf); a "2" means that C is used to compute the tacvf and not
#' the log-likelihood; and a "3" means C is used to compute everything.
#'
#' @param z The data set (time series)
#' @param order The order of the ARIMA model to be fit: c(p, d, q).  We have
#' that p is the number of AR parameters (phi), d is the amount of integer
#' differencing, and q is the number of MA parameters (theta).  Note we use the
#' Box-Jenkins convention for the MA parameters, in that they are the negative
#' of \code{\link{arima}}: see "Details".
#' @param numeach The number of starts to fit for each parameter.  The first
#' argument in the vector is the number of starts for each AR/MA parameter,
#' while the second is the number of starts for the fractional parameter.  When
#' this is set to 0, no fractional noise is fit.  Note that the number of
#' starts in total is multiplicative: if we are fitting an ARFIMA(2, d, 2), and
#' use the older number of starts (c(2, 2)), we will have 2^2 * 2 * 2^2 = 32
#' starting values for the fits.  \strong{Note that the default has changed
#' from c(2, 2) to c(1, 1) since package version 1.4-0}
#' @param back Setting this to true will restore the defaults in numeach.
#' @param dmean Whether the mean should be fit dynamically with the optimizer.
#' Note that the likelihood surface will change if this is TRUE, but this
#' is usually not worrisome.
#' See the referenced thesis for details.
#' @param whichopt Which optimizer to use in the optimization: see "Details".
#' @param itmean This option is under investigation, and will be set to FALSE
#' automatically until it has been decided what to do.
#'
#' Whether the mean should be fit iteratively using the function
#' \code{\link{TrenchMean}}.  Currently itmean, if set to TRUE, has higher
#' priority that dmean: if both are TRUE, dmean will be set to FALSE, with a
#' warning.
#' @param fixed A list of parameters to be fixed.  If we are to fix certain
#' elements of the AR process, for example, fixed$phi must have length equal to
#' p.  Any numeric value will fix the parameter at that value; for example, if
#' we are modelling an AR(2) process, and we wish to fix only the first
#' autoregressive parameter to 0, we would have fixed = list(phi = c(0, NA)).
#' NA corresponds to that parameter being allowed to change in the optimization
#' process.  We can fix the fractional parameters, and unlike
#' \code{\link{arima}}, can fix the seasonal parameters as well. Currently,
#' fixing regression/transfer function parameters is disabled.
#' @param lmodel The long memory model (noise type) to be used: "d" for FDWN,
#' "g" for FGN, "h" for PLA, and "n" for none (i.e. ARMA short memory models).
#' Default is "d".
#' @param seasonal The seasonal components of the model we wish to fit, with
#' the same components as above.  The period must be supplied.
#' @param useC How much interfaced C code to use: an integer between 0 and 3.
#' The value 3 is strongly recommended. See "Details".
#' @param cpus The number of CPUs used to perform the multi-start fits.  A
#' small number of fits and a high number of cpus (say both equal 4) with n not
#' large can actually be slower than when cpus = 1.  The number of CPUs should
#' not exceed the number of threads available to R.
#' @param rand Whether random starts are used in the multistart method.
#' Defaults to FALSE.
#' @param numrand The number of random starts to use.
#' @param seed The seed for the random starts.
#' @param eps3 How far to start from the boundaries when using a grid for the
#' multi-starts (i.e. when rand is FALSE.)
#' @param xreg A matrix, data frame, or vector of regressors for regression or
#' transfer functions.
#' @param reglist A list with the following elements:
#' \itemize{
#' \item regpar -
#' either NA or a list, matrix, data frame, or vector with 3 columns.  If
#' regpar is a vector, the matrix xreg must have one row or column only.  In
#' order, the elements of regpar are: r, s, and b.  The values of r are the the
#' orders of the delta parameters as in Box, Jenkins and Reinsel, the values of
#' s are the orders of omega parameters, and the values of b are the
#' backshifting to be done.
#'
#' \item minn - the minimum value for the starting value of the search, if
#' reglist$numeach > 1.
#' \item maxx - the maximum value for the starting value
#' of the search, if reglist$numeach > 1.
#' \item numeach - the number of starts
#' to try for each regression parameter.
#'
#' }
#' @param check If TRUE, checks at each optim iteration whether the model is
#' identifiable.  This makes the optimization much slower.
#' @param autoweed Whether to automatically (before the fit is returned) weed
#' out modes found that are found that are close together (usually the same
#' point.)
#' @param weedeps The maximum distance between modes that are close together
#' for the mode with the lower log-likelihood to be weeded out.  If adapt is
#' TRUE (default) this value changes.
#' @param adapt If TRUE, if dim is the dimensionality of the search, weedeps is
#' changed to \eqn{(1 + weedeps)^{dim} - 1}.
#' @param weedtype The type of weeding to be done.  See \code{\link{weed}}.
#' @param weedp The p in the p-norm to be used in the weeding.  p = 2 (default)
#' is Euclidean distance.
#' @param quiet If TRUE, no auxiliary output is generated. The default (FALSE)
#' has information of fits being proformed.
#' @param startfit Meant primarily for debugging (for now), allows starting places
#' for the fitting process.  Overrides \code{numeach}.
#' @return An object of class "arfima".  In it, full information on the fit is
#' given, though not printed under the print.arfima method.  The phis are the
#' AR parameters, and the thetas are the MA parameters.  Residuals, regression
#' residuals, etc., are all available, along with the parameter values and
#' standard errors.  Note that the muHat returned in the arfima object is
#' of the \strong{differenced} series, if differencing is applied.
#'
#' Note that if multiple modes are found, they are listed in order of
#' log-likelihood value.
#' @author JQ (Justin) Veenstra
#' @seealso \code{\link{arfima.sim}}, \code{\link{SeriesJ}},
#' \code{\link{arfima-package}}
#' @references McLeod, A. I., Yu, H. and Krougly, Z. L. (2007) Algorithms for
#' Linear Time Series Analysis: With R Package Journal of Statistical Software,
#' Vol. 23, Issue 5
#'
#' Veenstra, J.Q. Persistence and Antipersistence:  Theory and
#' Software (PhD Thesis)
#'
#' P. Borwein (1995) An efficient algorithm for Riemann Zeta function Canadian
#' Math. Soc. Conf. Proc., 27, pp. 29-34.
#' @keywords ts
#' @examples
#'
#' \donttest{
#' set.seed(8564)
#' sim <- arfima.sim(1000, model = list(phi = c(0.2, 0.1),
#' dfrac = 0.4, theta = 0.9))
#' fit <- arfima(sim, order = c(2, 0, 1), back=T)
#'
#' fit
#'
#' data(tmpyr)
#'
#' fit <- arfima(tmpyr, order = c(1, 0, 1), numeach = c(3, 3))
#' fit
#'
#' plot(tacvf(fit), maxlag = 30, tacf = TRUE)
#'
#' data(SeriesJ)
#' attach(SeriesJ)
#'
#' fitTF <- arfima(YJ, order= c(2, 0, 0), xreg = XJ, reglist =
#' list(regpar = c(2, 2, 3)), lmodel = "n")
#' fitTF
#'
#' detach(SeriesJ)
#'
#' }
#'
#' @export arfima
"arfima" <- function(z, order = c(0, 0, 0), numeach = c(1, 1), dmean = TRUE, whichopt = 0,
    itmean = FALSE, fixed = list(phi = NA, theta = NA, frac = NA, seasonal = list(phi = NA,
        theta = NA, frac = NA), reg = NA), lmodel = c("d", "g", "h", "n"), seasonal = list(order = c(0,
        0, 0), period = NA, lmodel = c("d", "g", "h", "n"), numeach = c(1, 1)), useC = 3,
    cpus = 1, rand = FALSE, numrand = NULL, seed = NA, eps3 = 0.01, xreg = NULL, reglist = list(regpar = NA,
        minn = -10, maxx = 10, numeach = 1), check = F, autoweed = TRUE, weedeps = 0.01,
    adapt = TRUE, weedtype = c("A", "P", "B"), weedp = 2, quiet = FALSE, startfit = NULL, back = FALSE) {

    if (itmean) {
        warning("itmean is currently under examination, and is not useable at the moment.\nSetting itmean to FALSE")
        itmean <- FALSE
    }
    if(back) {
      numeach <- c(2, 2)
      seasonal$numeach <- c(2, 2)
      dmean <- TRUE
    }
    if(!is.null(startfit)) {
      stopifnot(is.matrix(startfit))
      starts <- startfit
      num <- nrow(starts)
      startflag <- TRUE
    }
    else startflag <- FALSE

    if (any(is.na(z))) {
      stop("Missing data not supported at the moment.")
    }

    cpus <- abs(round(cpus))
    if (cpus < 1)
        cpus <- 1

    getHess <- TRUE  ## to save on options but not change code

    inclbdries <- FALSE  ## to save on options but not change code


    tse <- FALSE ## to save on options but not change code

    regpar <- reglist$regpar

    straightRegress <- FALSE
    namexreg <- NULL
    intindex <- 0
    intname <- NULL

    if (!is.null(regpar) && all(is.na(regpar)))
        regpar <- NULL
    if (!is.null(xreg)) {

        if (!(is.data.frame(xreg) || is.matrix(xreg) || is.vector(xreg) || is.ts(xreg)))
            stop("Invalid xreg type")


        if (any(is.na(xreg)))
            stop("no missing values allowed in xreg: if dynamic regression needs padding, please pad with zeroes")

        constantInds <- numeric(0)
        if (is.null(regpar)) {
            if (!is.data.frame(xreg) && !is.matrix(xreg)) {
              xreg <- data.frame(xreg)
              colnames(xreg) <- as.character(as.list((match.call()))$xreg)
            }

            namexreg <- colnames(xreg)

            if (is.null(namexreg))
              namexreg <- paste(rep("X", ncol(xreg)), 1:ncol(xreg), sep = "")
            regpar <- matrix(0, nrow = 3, ncol = dim(xreg)[2])
            regpar <- as.data.frame(t(regpar))
            names(regpar) <- c("r", "s", "b")
            constantInds <- c(1:dim(xreg)[2])[apply(xreg, 2, var) == 0]

            straightRegress <- ncol(xreg)
            if(length(constantInds)) {
              for(i in constantInds) {
                if(all(xreg[, i] == 1)) {
                  constantInds <- setdiff(constantInds, i)

                  if(i > 1) {
                    xreg <- xreg[, c(i, setdiff(1:ncol(xreg), i))]
                    namexreg <- namexreg[c(i, setdiff(1:ncol(xreg), i))]
                    intindex <- 1i
                    if(length(constantInds)>0)
                      constantInds <- c(1:dim(xreg)[2])[apply(xreg, 2, var) == 0][-1]
                  }
                  intindex <- 1
                  break
                }
              }
              straightRegress <- straightRegress - length(constantInds) - 1
            }
            straightRegress <- (straightRegress>0)

            if(!straightRegress) {
              warning("All regression terms constant.  Setting xreg to NULL")
              xreg <- NULL
            }
            else if(straightRegress && dmean != TRUE) {
              warning("Regression with time series errors requires dmean=T, setting")
              dmean <- TRUE
            }
            else if(straightRegress && intindex && is.numeric(dmean)) {
              warning("with an intercept term, dmean must be TRUE.  Setting it so.")
              dmean <- TRUE
            }
          }

         else {
            if (is.matrix(regpar) || is.data.frame(regpar)) {
                if (dim(regpar)[2] != dim(xreg)[2]) {
                  stop("xreg and regpar sizes do not match")
                }
                namexreg <- colnames(xreg)
                if (is.null(namexreg))
                  namexreg <- paste(rep("X", ncol(xreg)), 1:ncol(xreg), sep = "")
                if (any(is.na(regpar)))
                  stop("no missing values allowed in regpar")
                if (dim(regpar)[1] != 3)
                  stop("regpar does not have 3 rows")
                if (!is.null(names(regpar)) && names(regpar) != c("r", "s", "b"))
                  stop("names(regpar), if not NULL, must be r, s, and b")
                regpar <- as.data.frame(t(regpar))
                names(regpar) <- c("r", "s", "b")
            } else if (is.vector(xreg) || is.ts(xreg)) {
                pars <- data.frame(r = regpar[1], s = regpar[2], b = regpar[3])
                regpar <- pars
                xreg <- data.frame(xreg)
                colnames(xreg) <- as.character(as.list((match.call()))$xreg)
                namexreg <- colnames(xreg)
                if(is.null(namexreg))
                  namexreg <- "X1"
            } else stop("invalid regpar type")

        }
        xreg <- as.matrix(xreg)

        if (length(constantInds) > 0) {

          warning(paste("columns", paste(constantInds, collapse = ", "), "are constant and will be removed"))

          xreg <- xreg[, -constantInds]

          namexreg <- namexreg[-constantInds]

        }


        if(intindex && is.logical(dmean) && dmean) {
          regpar <- regpar[-1, ]
          xreg <- xreg[, -1]
          intname <- namexreg[1]
          namexreg <- namexreg[-1]
        }
        else if(intindex) {
          warning(paste(paste("there is a constant column named", namexreg[1]),
                  "take to be the intercept, but dmean is not TRUE.  This column",
                  "will be removed and ignored.", sep = "\n"))
          regpar <- regpar[-1, ]
          xreg <- xreg[, -1]
          intindex <- 0
        }
        else if(is.logical(dmean) && dmean) {
          intname <- "Intercept"
        }

        s <- regpar$s <- regpar$s + 1
        r <- regpar$r
        b <- regpar$b

        if (any(r > 0) || any(s>1)) {
            #cat("note: transfer functions do not work with dynamic mean:  setting dmean to FALSE\n")
            dmean <- FALSE
        }
        numvarreg <- sum(s) + sum(r)
        regeach <- if(straightRegress) 1 else reglist$numeach
        if (length(regeach) == 0)
            regeach <- 1
        if (length(regeach) != 1 || regeach != round(regeach) || regeach <= 0) {
            warning("invalid reglist$numeach: resetting to 1")
            regeach <- 1
        }

        regmin <- reglist$minn
        regmax <- reglist$maxx
        if (length(regmin) == 0 || is.na(regmin))
            regmin <- if(regeach>1) -10 else 0
        if (length(regmax) == 0 || is.na(regmax))
            regmax <- if(regeach>1) 10 else 0
        if (itmean)
            stop("Please note that for regression problems, the mean is calculated separately\n via optimization or on the regression parameters or by the mean of the series;\n therefore itmean is not valid")
      }
    else {
      r <- s <- b <- 0
      numvarreg <- regeach <- 0
      namexreg <- reglist <- regpar <- NULL
    }

    if (is.numeric(dmean) && !is.double(dmean))
        stop("dmean must be logical or a double")

    if (itmean && ((is.logical(dmean) && dmean) || is.numeric(dmean))) {
        warning("itmean trumps dmean: setting dmean to FALSE")
        dmean <- FALSE
    }

    info <- NULL
    dint <- order[2]
    if (length(dint) == 0)
        dint <- 0
    lmodel <- lmodel[1]
    if (is.null(numeach) || (length(numeach) != 2 && (lmodel != "n" && length(numeach) !=
        1)) || !is.vector(numeach))
        stop("numeach must be a vector of length 1 or 2:\n \n\t\tposition 1 is the values for ARMA\n parameters, position 2 is the value for the fractional parameter\n\n\t\tif a vector of length 1, we must have lmodel = \"n\"")
    if (numeach[1] != round(numeach[1]) || numeach[1] <= 0) {
        warning("incorrect numeach ARMA: resetting to 1")
        numeach[1] <- 1
    }
    if (numeach[2] != round(numeach[2]) || numeach[2] < 0) {
        warning("incorrect numeach fractional: resetting to 1")
        numeach[2] <- 1
    }

    if (numeach[2] == 0 && lmodel != "n") {
        cat("lmodel = \"", lmodel, "\" and numeach[2] = 0; setting lmodel = \"n\"\n", sep = "")
        lmodel <- "n"
    }

    if (lmodel == "n") {
        if (length(numeach) == 1)
            numeach <- c(numeach, 0) else numeach[2] <- 0
    }

    if (!is.null(seasonal) && !is.na(seasonal$period) && !is.null(seasonal$period) && (seasonal$period !=
        0) && ((!is.null(seasonal$order) || !all(seasonal$order == 0)) || (!is.null(seasonal$lmodel) &&
        seasonal$lmodel != "n"))) {
        dseas <- seasonal$order[2]
        if (length(dseas) == 0)
            dseas <- 0
        pseas <- seasonal$order[1]
        qseas <- seasonal$order[3]
        if (length(pseas) == 0)
            pseas <- 0
        if (length(qseas) == 0)
            qseas <- 0

        if (!is.null(seasonal$lmodel))
            slmodel <- seasonal$lmodel[1] else slmodel <- "d"
        if (!is.null(seasonal$numeach))
            snumeach <- seasonal$numeach else {
            warning("no seasonal$numeach:  setting to numeach")
            snumeach <- numeach
        }
        if (is.null(snumeach) || (length(snumeach) != 2 && (slmodel != "n" && length(snumeach) !=
            1)) || !is.vector(snumeach))
            stop("seasonal$numeach ust be a vector of length 1 or 2:\n \n\t\tposition 1 is the values for ARMA\n parameters, position 2 is the value for the fractional parameter\n\n\t\tif a vector of length 1, we must have seasonal$lmodel = \"n\"")
        period <- seasonal$period

        if (snumeach[1] != round(snumeach[1]) || snumeach[1] <= 0) {
            warning("incorrect seasonal$numeach ARMA: resetting to 1")
            snumeach[1] <- 1
        }

        if (snumeach[2] != round(snumeach[2]) || snumeach[2] < 0) {
            warning("incorrect seasonal$numeach fractional: resetting to 1")
            snumeach[2] <- 1
        }

        if (snumeach[2] == 0 && slmodel != "n") {
            cat("seasonal lmodel =", slmodel, "and seasonal numeach[2] = 0; setting seasonal lmodel = \"n\"\n")
            slmodel <- "n"
        }

        if (slmodel == "n") {
            if (length(snumeach) == 1)
                snumeach <- c(snumeach, 0) else snumeach[2] <- 0
        }
        slmr <- if (slmodel != "n")
            1 else 0
    } else {
        dseas <- pseas <- qseas <- period <- slmr <- 0
        snumeach <- c(0, 0)
        slmodel <- "n"
    }

    p <- order[1]
    q <- order[3]
    lmr <- if (lmodel != "n")
        1 else 0


    numvars <- p + q + pseas + qseas + lmr + slmr + numvarreg +
      if (dmean == TRUE)
        1 else 0

    if (length(fixed) > 0) {
        if (!(is.null(fixed$reg)) && !(all(is.na(fixed$reg))) && length(fixed$reg) > 0) {
            warning("no fixed regression parameters for now.")
            fixed$reg <- NA
        }
        arfixed <- fixed$phi
        mafixed <- fixed$theta
        fracfix <- fixed$frac
        if((length(arfixed)==1 && is.na(arfixed))||length(arfixed)==0)
          arfixed <- rep(NA, p)
        if((length(mafixed)==1 && is.na(mafixed))||length(mafixed)==0)
          mafixed <- rep(NA, q)
        if(length(fracfix)==0 && lmodel !="n") {
          fracfix <- NA
        }
        warnstring <- "In the fixed argument:"
        if (length(arfixed) != p) {
          warnstring <- paste0(warnstring, "\n length phi is ", length(arfixed),
                               " but p = ", p, " unfixing phi")
          arfixed <- if (p > 0)
            rep(NA, p) else numeric(0)

        }
        else if(p>0 && length(arfixed)==p&&any(!is.na(arfixed))) {
          if(sum(is.na(arfixed))==0) {
            if(!InvertibleQ(arfixed)) {
              stop("phi parameters fixed at nonstationary parameters")
            }
          }
        }
        if (length(mafixed) != q) {
          warnstring <- paste0(warnstring, "\n length theta is ", length(mafixed),
                               " but q = ", q, " unfixing theta")
          mafixed <- if (q > 0)
            rep(NA, q) else numeric(0)

        }
        else if(any(!is.na(mafixed))) {
          if(sum(is.na(mafixed))==0) {
            if(!InvertibleQ(mafixed)) {
              stop("theta parameters fixed at noninvertible parameters")
            }
          }
        }

        if (!is.null(fracfix) && length(fracfix) != 1 && !is.na(fracfix)) {
          warnstring <- paste0(warnstring, "\n fractional parameter is fixed greater",
                               "\n than length 1, resetting fractional prameter ")
          if(lmodel=="n") {
            warnstring<- paste0(warnstring,"\n lmodel == 'n' and fractional parameter fixed, resetting")
            fracfix <- numeric(0)
          }
          else fracfix <- NA
        }
        if (lmodel == "n")
            fracfix <- numeric(0)
        if(lmodel!="n"&&(!is.null(fracfix)&&!is.na(fracfix))) {
          fix_stop <- F
          if(lmodel=="d") {
            if(!InvertibleD(fracfix)) fix_stop <- T
          }
          else if(lmodel=="g") {
            if(!InvertibleH(fracfix)) fix_stop <- T
          }
          else {
            if(!InvertibleAlpha(fracfix)) fix_stop <- T
          }
          if(fix_stop){
            stop("fractional parameter fixed at nonstationary/noninvertible value")
          }

        }

        fixx <- c(arfixed, mafixed)
        if (period > 0) {
          sarfixed <- fixed$seasonal$phi
          smafixed <- fixed$seasonal$theta
          sfracfix <- fixed$seasonal$frac

          if((length(sarfixed)==1 && is.na(sarfixed))||length(sarfixed)==0)
            sarfixed <- rep(NA, pseas)
          if((length(smafixed)==1 && is.na(smafixed))||length(smafixed)==0)
            smafixed <- rep(NA, qseas)
          if(length(sfracfix)==0 && slmodel !="n") {
            sfracfix <- NA
          }
          if (length(sarfixed) != pseas) {
            warnstring <- paste0(warnstring, "\n length seasonal phi is ", length(sarfixed),
                                 " but pseas = ", pseas, " unfixing s_phi")
            sarfixed <- if (pseas > 0)
              rep(NA, pseas) else numeric(0)

          }
          else if(pseas>0 && length(sarfixed)==pseas&&any(!is.na(sarfixed))) {
            if(sum(is.na(sarfixed))==0) {
              if(!InvertibleQ(sarfixed)) {
                stop("seasonal phi parameters fixed at nonstationary parameters")
              }
            }
          }

          if (length(smafixed) != qseas) {
            warnstring <- paste0(warnstring, "\n length seasonal theta is ", length(smafixed),
                                 " but qseas = ", qseas, " unfixing s_theta")
            smafixed <- if (qseas > 0)
              rep(NA, qseas) else numeric(0)

          }
          else if(qseas>0 && length(smafixed)==qseas&&any(!is.na(smafixed))) {
            if(sum(is.na(smafixed))==0) {
              if(!InvertibleQ(smafixed)) {
                stop("seasonal theta parameters fixed at noninvertible parameters")
              }
            }
          }

          if (!is.null(sfracfix) && length(sfracfix) != 1 && !is.na(sfracfix)) {
            warnstring <- paste0(warnstring, "\n seasonal fractional parameter is",
                                 "\n fixed greater than length 1, resetting seasonal",
                                 "\n fractional prameter ")
            if(slmodel=="n") {
              warnstring("\n slmodel == 'n' and seasonal fractional parameter fixed",
                        "\n resetting")
              sfracfix <- numeric(0)
            }
            else sfracfix <- NA
          }
          if(slmodel!="n"&&(!is.null(sfracfix)&&!is.na(sfracfix))) {
            fix_stop <- F
            if(slmodel=="d") {
              if(!InvertibleD(sfracfix)) fix_stop <- T
            }
            else if(slmodel=="g") {
              if(!InvertibleH(sfracfix)) fix_stop <- T
            }
            else {
              if(!InvertibleAlpha(sfracfix)) fix_stop <- T
            }
            if(fix_stop) {
              stop("seasonal fractional parameter fixed at nonstationary/noninvertible value")
            }
          }
          sfixx <- c(sarfixed, smafixed)
        } else {
            sfixx <- NULL
            sfracfix <- NULL
        }

        if(nchar(warnstring)>31) warning(warnstring)

        fixx <- c(fixx, sfixx, fracfix, sfracfix, rep(NA, sum(r)+sum(s)))
        # if((sum(r)+sum(s)==0)&&dmean==T) {
        #   fixx <- c(fixx, NA)
        # }
        indfixx <- !is.na(fixx)

        if (is.logical(dmean) && dmean) {
            fixx <- c(fixx, NA)
            indfixx <- c(indfixx, FALSE)
        } else if (is.logical(itmean) && itmean) {
            fixx <- c(fixx, NA)
            indfixx <- c(indfixx, FALSE)
        }

    } else {
        fixx <- rep(NA, numvars)
        indfixx <- rep(FALSE, numvars)

        if (is.logical(dmean) && dmean) {
            fixx <- c(fixx, NA)
            indfixx <- c(indfixx, FALSE)
        } else if (is.logical(itmean) && itmean) {
            fixx <- c(fixx, NA)
            indfixx <- c(indfixx, FALSE)
        }
    }

    numtries <- 1

    if (p + q > 0)
        numtries <- numtries * numeach[1]^(p + q)
    if (numeach[2] > 1)
        numtries <- numtries * numeach[2]
    if (pseas + qseas > 0)
        numtries <- numtries * snumeach[1]^(pseas + qseas)
    if (snumeach[2] > 1)
        numtries <- numtries * snumeach[2]
    if (regeach > 1)
        numtries <- numtries * regeach

    if (sum(as.integer(indfixx)) > 0) {

        if (p > 0 && any(!is.na(arfixed)))
            numtries <- numtries/(numeach[1]^(sum(!is.na(arfixed))))
        if (q > 0 && any(!is.na(mafixed)))
            numtries <- numtries/(numeach[1]^(sum(!is.na(mafixed))))
        if (lmodel != "n" && (is.null(fracfix)||!is.na(fracfix)))
            numtries <- numtries/numeach[2]
        if (pseas > 0 && any(!is.na(sarfixed)))
            numtries <- numtries/(snumeach[1]^(sum(!is.na(sarfixed))))
        if (qseas > 0 && any(!is.na(smafixed)))
            numtries <- numtries/(snumeach[1]^(sum(!is.na(smafixed))))
        if (slmodel != "n" && (is.null(sfracfix)||!is.na(sfracfix)))
            numtries <- numtries/snumeach[2]
    }

    # if(numvars == 0 && numtries > 1) stop('invalid specification of variables')
    if (round(dint) != dint || dint < 0)
        stop("d must be an integer >= 0")
    if (round(dseas) != dseas || dseas < 0)
        stop("dseas must be an integer >= 0")

    if (!quiet) {
      if(!startflag) {
        if (autoweed && numtries > 1 && !startflag)
          cat("Note: autoweed is ON.  It is possible, but not likely,\n that unique modes may be lost.\n") else if (!autoweed && numtries > 1)
            cat("Note: autoweed is OFF.  Running weed after the\n fit is highly recommended.\n") else {
              cat("Note: only one starting point.  Only one mode can be found -- this is now the default behavior.\n")
            }
      }
      else {
        if(autoweed)
          cat("Note: autoweed is ON.  It is possible, but not likely,\n that unique modes may be lost.\n")
        else
          cat("Note: autoweed is OFF.  Running weed after the\n fit may be required.\n")
        cat("Starting fit with given start matrix.\n")
      }
    }
    y <- z

    xr <- xreg

    if (dint > 0) {
        y <- diff(y, differences = dint)
        if (!is.null(xreg))
            xreg <- diff(xreg, differences = dint)
    }
    if (dseas > 0) {
        y <- diff(y, differences = dseas, lag = period)
        if (!is.null(xreg))
            xreg <- diff(xreg, differences = dseas, lag = period)
    }
    D <- dint + dseas * period


    ## Should this go before differencing?
    if (any(r > 0) && any(abs(colMeans(xreg)) > 1e-14)) {
        cat("Please note for transfer functions the means of each X variable must be 0: subtracting mean from each X\n")
        for (i in 1:ncol(xreg)) {
            xreg[, i] <- xreg[, i] - colMeans(xreg)[i]
        }
    }


    flag0 <- regOnly <- FALSE
    if (p + q + numeach[2] + pseas + qseas + snumeach[2] == 0) {
      if(is.null(xreg))
        flag0 <- TRUE
      else
        regOnly <- TRUE
    }

    differencing <- if (D > 0) TRUE else FALSE



    if (is.logical(dmean) && !dmean && !straightRegress)
      meanval <- mean(y)
    else if (is.numeric(dmean))
      meanval <- dmean
    else
      meanval <- 0


    sumfixed <- sum(indfixx)

    numvars1 <- numvars - sumfixed
    if (numvars1 == 0)
        getHess <- FALSE


    if (flag0) {
        if(sum(numeach)>2||sum(seasonal$numeach)>2||
           (!is.null(numrand)) && numrand>0)
          warning('With only the mean to be fit, since the problem is convex,
                  it has a unique minimum and thus only one fit will be performed.')
        if ((is.logical(dmean) && dmean) || itmean) {
            fun <- function(m) {
                lARFIMA(y - m)
            }
            fit <- optim(mean(y), fn = fun, method = "Brent", lower = mean(y) - 10, upper = mean(y) +
                10, control = list(trace = 0, fnscale = -1), hessian = getHess)
            hess <- fit$hessian[1, 1]
            se <- tryCatch(sqrt(diag(solve(-hess))), error = function(err) rep(NA, dim(hess)[1]))
            res <- list(maximum = fit[[1]][1], objective = fit[[2]][1], hess = hess, se = se)
        } else if (is.logical(dmean)) {
            res <- list(maximum = mean(y), objective = lARFIMA(y - mean(y)), hess = NA,
                se = -Inf)
        } else {
            res <- list(maximum = dmean, objective = lARFIMA(y - dmean), hess = NA, se = -Inf)
        }
        allpars <- vector("list", 1)
        allpars[[1]]$muHat <- res$maximum
        allpars[[1]]$loglik <- res$objective
        allpars[[1]]$hessian <- res$hess
        allpars[[1]]$se <- res$se
        allpars[[1]]$phi <- allpars[[1]]$theta <- allpars[[1]]$phiseas <-
          allpars[[1]]$thetaseas <- allpars[[1]]$dfrac <- allpars[[1]]$dfs <-
          allpars[[1]]$phip <- allpars[[1]]$thetap <- allpars[[1]]$H <-
          allpars[[1]]$alpha <- allpars[[1]]$alphas <- allpars[[1]]$delta <-
          allpars[[1]]$omega <- allpars[[1]]$Hs <- numeric(0)
        allpars[[1]]$pars <-
          c(allpars[[1]]$phi, allpars[[1]]$theta, allpars[[1]]$phiseas,
            allpars[[1]]$thetaseas, allpars[[1]]$dfrac, allpars[[1]]$H,
            allpars[[1]]$alpha, allpars[[1]]$dfs, allpars[[1]]$Hs,
            allpars[[1]]$alphas, allpars[[1]]$delta, allpars[[1]]$omega)
        sigma2muse <- getsigma2muhat(y - allpars[[1]]$muHat)
        allpars[[1]]$sigma2 <- sigma2muse
        rr <- tacvfARFIMA(maxlag = length(y) - 1)
        res <- DLResiduals(rr, y - allpars[[1]]$muHat)
        allpars[[1]]$residuals <- res
        allpars[[1]]$fitted <- y - allpars[[1]]$muHat - res
        class(allpars[[1]]) <- "ARFIMA"
        ans <- list(z = z, differencing = differencing, dint = dint, dseas = dseas, period = period,
                    modes = allpars, dmean = dmean, itmean = itmean, p = p, q = q, pseas = pseas,
                    qseas = qseas, lmodel = lmodel, slmodel = slmodel, weeded = TRUE, getHess = getHess,
                    numvars = numvars, numcut = 0, n = length(z), xreg = xr, r = r, s = s, b = b,
                    call = match.call(), flag0 = flag0, numeach = numeach, strReg = F,
                    regOnly = regOnly, intindex = intindex, intname = intname, namexreg=namexreg,
                    fixx = fixx, indfixx = indfixx)
        # ans <- list(z = z, differencing = differencing, dint = dint, dseas = dseas, period = period,
        #     modes = allpars, dmean = dmean, itmean = itmean, p = p, q = q, pseas = pseas,
        #     qseas = qseas, lmodel = lmodel, slmodel = slmodel, weeded = TRUE, getHess = getHess,
        #     numvars = numvars, numcut = 0, n = length(z), xreg = xr, r = r, s = s, b = b,
        #     call = match.call(), flag0 = flag0)
    } else {
      num <- if (!rand)
        numtries else numrand
      if(startflag) num <- nrow(starts)

      if (rand && !startflag) {
        if (length(num) == 0) {
          warning("random starts selected and no number of starts selected: setting to 8")
          num <- 8
        }
        starts <- matrix(0, nrow = num, ncol = numvars)
        pl1 <- 0
        if (!is.na(seed) && !is.null(seed))
          set.seed(seed) else warning("Random start seed not selected:  results may not be reproducible")
        if (p + q + pseas + qseas > 0) {
          starts[, 1:(p + q + pseas + qseas)] <- matrix(runif(num * (p + q + pseas +
                                                                       qseas), min = -1, max = 1), ncol = (p + q + pseas + qseas))
        }

        if (lmodel != "n") {
          pl1 <- 1
          if (lmodel == "d")
            starts[, (1 + p + q + pseas + qseas)] <- runif(1, min = -1, max = 0.5)
          else if (lmodel == "g")
            starts[, (1 + p + q + pseas + qseas)] <- runif(1, min = 0, max = 1)
          else
            starts[, (1 + p + q + pseas + qseas)] <- runif(1, min = 0, max = 3)
        }
        if (slmodel != "n") {
          if (slmodel == "d")
            starts[, (1 + p + q + pl1 + pseas + qseas)] <- runif(1, min = -1, max = 0.5)
          else if (slmodel == "g")
            starts[, (1 + p + q + pl1 + pseas + qseas)] <- runif(1, min = 0, max = 1)
          else
            starts[, (1 + p + q + pl1 + pseas + qseas)] <- runif(1, min = 0, max = 3)
        }
        if (!is.null(xreg)) {
          if(!straightRegress) {
            starts[, (1 + p + q + pl1 + pseas + qseas):(1 + p + q + pl1 + pseas + qseas +
                                                          numvarreg)] <- matrix(runif(num * (numvarreg), min = regmin,
                                                                                      max = regmax), nrow = num)
          }
          else if(is.logical(dmean)&&dmean) {
            coeffs <- lsfit(xreg, y)$coef
            for(ii in 1:numvarreg) starts[, (1 + p + q + pl1 + pseas + qseas+ii)] <- coeffs
          }
          else {
            coeffs <- lsfit(xreg, y, intercept = FALSE)$coef
            for(ii in 1:numvarreg) starts[, (1 + p + q + pl1 + pseas + qseas+ii)] <- coeffs
          }
        }


      } else if (num > 1 && !startflag) {
        pl1 <- 0

        li <- vector("list", numvars)
        if (p + q > 0) {
          if (inclbdries)
            seqq <- seq(-1 + eps3, 1 - eps3, length.out = numeach[1]) else {
              seqq <- seq(-1 + eps3, 1 - eps3, length.out = numeach[1] + 2)
              seqq <- seqq[-1]
              seqq <- seqq[-length(seqq)]
            }

        }
        if (pseas + qseas > 0) {
          if (inclbdries)
            seqq <- seq(-1 + eps3, 1 - eps3, length.out = snumeach[1]) else {
              seqq <- seq(-1 + eps3, 1 - eps3, length.out = snumeach[1] + 2)
              seqq <- seqq[-1]
              seqq <- seqq[-length(seqq)]
            }
        }
        for (i in 1:(p + q + pseas + qseas))
          li[[i]] <- seqq
        if (lmodel != "n") {
          pl1 <- 1
          if (lmodel == "d") {
            minner <- -1
            maxxer <- 0.5
          } else if (lmodel == "g") {
            minner <- 0
            maxxer <- 1
          } else {
            minner <- 0
            maxxer <- 3
          }
          if (inclbdries)
            seqq <- seq(minner + eps3, maxxer - eps3, length.out = numeach[2]) else {
              seqq <- seq(minner + eps3, maxxer - eps3, length.out = numeach[2] +
                            2)
              seqq <- seqq[-1]
              seqq <- seqq[-length(seqq)]
            }
          li[[(p + q + pl1 + pseas + qseas)]] <- seqq
        }
        pl1_1 <- 0
        if (slmodel != "n") {
          pl1_1 <- 1
          if (slmodel == "d") {
            minner <- -1
            maxxer <- 0.5
          } else if (slmodel == "g") {
            minner <- 0
            maxxer <- 1
          } else {
            minner <- 0
            maxxer <- 3
          }
          if (inclbdries)
            seqq <- seq(minner + eps3, maxxer - eps3, length.out = snumeach[2]) else {
              seqq <- seq(minner + eps3, maxxer - eps3, length.out = snumeach[2] +
                            2)
              seqq <- seqq[-1]
              seqq <- seqq[-length(seqq)]
            }
          li[[(p + q + 1 + pseas + qseas + pl1)]] <- seqq
        }
        if (!is.null(xreg)) {
          if(!straightRegress) {
            if (inclbdries)
              seqq <- seq(regmin + eps3, regmax - eps3, length.out = regeach)
            else {
              seqq <- seq(regmin + eps3, regmax - eps3, length.out = regeach + 2)
              seqq <- seqq[-1]
              seqq <- seqq[-length(seqq)]
            }

            for (i in 1:(numvarreg)) li[[(p + q + pseas + qseas + pl1 + pl1_1) + i]] <- seqq
          }
          else if(is.logical(dmean)&&dmean){
            coeffs <- lsfit(xreg, y)$coef[c(2:(ncol(xreg)+1), 1)]
            for (i in 1:(numvarreg+1)) li[[(p + q + pseas + qseas + pl1 + pl1_1) + i]] <- coeffs[i]
          }
          else {
            coeffs <- lsfit(xreg, y, intercept = FALSE)$coef
            for (i in 1:(numvarreg)) li[[(p + q + pseas + qseas + pl1 + pl1_1) + i]] <- coeffs[i]
          }
        }
        else if(is.logical(dmean) && dmean)
          li[[length(li)]] <- mean(y)

        starts <- expand.grid(li, KEEP.OUT.ATTRS = FALSE)

        if(sumfixed>0) {
          for(i in 1:numvars) {
            if(indfixx[i]) starts[, i] <- rep(fixx[i], nrow(starts))
          }
          starts <- unique(starts)
        }
      } else if(!startflag) {

        starts <- if (numvars > 0)
          matrix(0, nrow=1,ncol = numvars) else matrix(0)
        if(straightRegress && dmean == TRUE) {
          coeffs <- lsfit(xreg, y)$coef[c(2:(ncol(xreg)+1), 1)]
          starts[, (numvars-length(coeffs)+1):numvars] <- coeffs
        }
        else if(straightRegress) {
          coeffs <- lsfit(xreg, y, intercept = FALSE)$coef
          starts[, (numvars-length(coeffs)+1):numvars] <- coeffs
        }
        else if(dmean == TRUE)
          starts[, numvars] <- mean(y)
        if(sumfixed>0) {
          for(i in 1:numvars) {
            if(indfixx[i]) starts[, i] <- fixx[i]
          }
        }
      }

      if(dmean==F||is.numeric(dmean)) {
        indfixxM <- c(indfixx, TRUE)
        numvarsM <- numvars + 1
      }
      else {
        indfixxM <- indfixx
        numvarsM <- numvars
      }

      if (nrow(starts) != num)
        stop("error in starts")
      if (rand)
        strg <- " the random start " else strg <- " the "
      if (!quiet) {
        cat("Beginning", strg, "fits with ", num, " starting values.\n", sep = "")
        cat("\n")
      }
      fit <- vector("list", num)

      if(!straightRegress)
        xreg <- as.vector(xreg)


        if (cpus == 1) {
            for (i in 1:num) {

                fit[[i]] <- arfimaFit(y = y, p = p, q = q, pseas = pseas, qseas = qseas,
                  period = period, lmodel = lmodel, slmodel = slmodel, whichopt = whichopt,
                  dmean = dmean, parstart = as.double(starts[i, ]), getHess = getHess, itmean = itmean,
                  indfixx = indfixx, fixx = fixx, xreg = xreg, r = r, s = s, b = b, straightRegress = straightRegress)
            }
        } else {
            cl <- makePSOCKcluster(cpus, methods = FALSE)
            numms <- 1:num
            clusterExport(cl, c("lARFIMA", "PacfToAR", "DLLoglikelihood", "lARFIMAwTF"),
                envir = environment(arfima))
            fit <- clusterApplyLB(cl, numms, arfimaFitpar, parstart = starts, y = y, p = p,
                q = q, pseas = pseas, qseas = qseas, period = period, lmodel = lmodel, slmodel = slmodel,
                whichopt = whichopt, useC = useC, dmean = dmean, getHess = getHess, itmean = itmean,
                indfixx = indfixx, fixx = fixx, xreg = xreg, r = r, s = s, b = b, straightRegress = straightRegress)
            stopCluster(cl)
            setDefaultCluster(NULL)
        }

        allpars <- vector("list", num)
        logl <- rep(0, num)

        for (i in 1:num) {

            pars <- fit[[i]][[1]]

            pars[indfixx] <- fixx[indfixx]
            allpars[[i]]$phi <- if (p > 0)
                pars[1:p] else numeric(0)
            allpars[[i]]$theta <- if (q > 0)
                pars[p + (1:q)] else numeric(0)
            allpars[[i]]$phiseas <- if (pseas > 0)
                pars[p + q + (1:pseas)] else numeric(0)
            allpars[[i]]$thetaseas <- if (qseas > 0)
                pars[p + q + pseas + (1:qseas)] else numeric(0)
            allpars[[i]]$phip <- if (p > 0)
                ARToPacf(pars[1:p]) else numeric(0)
            allpars[[i]]$thetap <- if (q > 0)
                ARToPacf(pars[p + (1:q)]) else numeric(0)
            allpars[[i]]$phiseasp <- if (pseas > 0)
                ARToPacf(pars[p + q + (1:pseas)]) else numeric(0)
            allpars[[i]]$thetaseasp <- if (qseas > 0)
                ARToPacf(pars[p + q + pseas + (1:qseas)]) else numeric(0)

            allpars[[i]]$dfrac <- if (lmodel == "d")
                pars[p + q + pseas + qseas + 1] else numeric(0)
            allpars[[i]]$H <- if (lmodel == "g")
                pars[p + q + pseas + qseas + 1] else numeric(0)
            allpars[[i]]$alpha <- if (lmodel == "h")
                pars[p + q + pseas + qseas + 1] else numeric(0)

            allpars[[i]]$dfs <- if (slmodel == "d")
                pars[p + q + pseas + qseas + length(allpars[[i]]$dfrac) + length(allpars[[i]]$H) +
                  length(allpars[[i]]$alpha) + 1] else numeric(0)
            allpars[[i]]$Hs <- if (slmodel == "g")
                pars[p + q + pseas + qseas + length(allpars[[i]]$dfrac) + length(allpars[[i]]$H) +
                  length(allpars[[i]]$alpha) + 1] else numeric(0)
            allpars[[i]]$alphas <- if (slmodel == "h")
                pars[p + q + pseas + qseas + length(allpars[[i]]$dfrac) + length(allpars[[i]]$H) +
                  length(allpars[[i]]$alpha) + 1] else numeric(0)

            allpars[[i]]$delta <- if (!is.null(xreg) && sum(r) > 0)
                pars[p + q + pseas + qseas + length(allpars[[i]]$dfrac) + length(allpars[[i]]$H) +
                  length(allpars[[i]]$dfs) + length(allpars[[i]]$Hs) + length(allpars[[i]]$alpha) +
                  length(allpars[[i]]$alphas) + 1:(sum(r))] else numeric(0)
            allpars[[i]]$omega <- if (!is.null(xreg) && sum(s) > 0)
                pars[p + q + pseas + qseas + length(allpars[[i]]$dfrac) + length(allpars[[i]]$H) +
                  length(allpars[[i]]$dfs) + length(allpars[[i]]$Hs) + length(allpars[[i]]$alpha) +
                  length(allpars[[i]]$alphas) + sum(r) + 1:sum(s)] else numeric(0)
            allpars[[i]]$muHat <- if (is.logical(dmean) && dmean)
                pars[p + q + pseas + qseas + length(allpars[[i]]$dfrac) + length(allpars[[i]]$H) +
                  length(allpars[[i]]$dfs) + length(allpars[[i]]$Hs) + length(allpars[[i]]$alpha) +
                  length(allpars[[i]]$alphas) + sum(r) + sum(s) + 1] else meanval

            rr <- tacvfARFIMA(phi = allpars[[i]]$phi, theta = allpars[[i]]$theta, phiseas = allpars[[i]]$phiseas,
                thetaseas = allpars[[i]]$thetaseas, dfrac = allpars[[i]]$dfrac, dfs = allpars[[i]]$dfs,
                H = allpars[[i]]$H, Hs = allpars[[i]]$Hs, alpha = allpars[[i]]$alpha, alphas = allpars[[i]]$alphas,
                period = period, maxlag = length(y) - 1)


            if (itmean) {
                allpars[[i]]$muHat <- TrenchMean(rr, y)
            }
            if (is.null(fit[[i]][2])) {
                cat("mode ", i, " has no loglikelihood value", sep = "")
                allpars[[i]]$loglik <- logl[i] <- -3e+10
            } else {
                allpars[[i]]$loglik <- logl[i] <- fit[[i]][[2]]
            }

            allpars[[i]]$pars <- c(allpars[[i]]$phi, allpars[[i]]$theta, allpars[[i]]$phiseas,
                allpars[[i]]$thetaseas, allpars[[i]]$dfrac, allpars[[i]]$H, allpars[[i]]$alpha,
                allpars[[i]]$dfs, allpars[[i]]$Hs, allpars[[i]]$alphas, allpars[[i]]$delta,
                allpars[[i]]$omega)

            allpars[[i]]$parsp <- c(allpars[[i]]$phip, allpars[[i]]$thetap, allpars[[i]]$phiseasp,
                allpars[[i]]$thetaseasp, allpars[[i]]$dfrac, allpars[[i]]$H, allpars[[i]]$alpha,
                allpars[[i]]$dfs, allpars[[i]]$Hs, allpars[[i]]$alphas, allpars[[i]]$delta,
                allpars[[i]]$omega)

            if (getHess) {
                hess <- fit[[i]]$hessian[!indfixx, !indfixx]

                se <- rep(-Inf, numvarsM)
                se[!indfixxM] <- tryCatch(sqrt(diag(solve(-hess))),
                              error = function(err) rep(NA, max(dim(hess)[1], 1)))

                allpars[[i]]$se <- se
                allpars[[i]]$hess <- hess

            } else hess <- numeric(0)

            if (!is.null(xreg)) {
                omega <- allpars[[i]]$omega
                delta <- allpars[[i]]$delta

                xreglist <- list(s = s, r = r, b = b, omega = vector("list", length(r)),
                  delta = vector("list", length(r)))
                if (!straightRegress) {
                  for (j in 1:length(r)) {
                    if (r[j] == 0) {
                      xreglist$delta[[j]] <- NA
                      if (s[j] == 1) {
                        if (j == 1)
                          xreglist$omega[[j]] <- omega[1:s[1]] else {
                          ss <- cumsum(s)[j - 1]
                          xreglist$omega[[j]] <- omega[(ss + 1):(ss + s[j])]
                          }
                        names(xreglist$omega[[j]]) <- namexreg[j]
                      } else {
                        if (j == 1)
                          xreglist$omega[[j]] <- omega[1:s[1]] else {
                          ss <- cumsum(s)[j - 1]
                          xreglist$omega[[j]] <- omega[(ss + 1):(ss + s[j])]
                        }
                        names(xreglist$omega[[j]]) <- paste("omega(", 0:(s[j] - 1), ").",
                          namexreg[j], sep = "")
                      }
                    } else {
                      if (j == 1)
                        xreglist$delta[[j]] <- delta[1:r[1]] else {
                        rrr <- cumsum(r)[j - 1]
                        xreglist$delta[[j]] <- delta[(rrr + 1):(rrr + r[j])]
                      }
                      names(xreglist$delta[[j]]) <- paste("delta(", 1:r[j], ").", namexreg[j],
                        sep = "")
                      if (s[j] == 1) {
                        if (j == 1)
                          xreglist$omega[[j]] <- omega[1:s[1]] else {
                          ss <- cumsum(s)[j - 1]
                          xreglist$omega[[j]] <- omega[(ss + 1):(ss + s[j])]
                        }
                        names(xreglist$omega[[j]]) <- paste("omega(0).", namexreg[j], sep = "")
                      } else {
                        if (j == 1)
                          xreglist$omega[[j]] <- omega[1:s[1]] else {
                          ss <- cumsum(s)[j - 1]
                          xreglist$omega[[j]] <- omega[(ss + 1):(ss + s[j])]
                        }
                        names(xreglist$omega[[j]]) <- paste("omega(", 0:(s[j] - 1), ").",
                          namexreg[j], sep = "")
                      }
                    }
                  }
                }
                else if (straightRegress) {
                  xreglist$omega <- as.list(omega)
                  names(xreglist$omega) <- namexreg
                }

                allpars[[i]]$xreglist <- xreglist

                mmu <- allpars[[i]]$muHat

                if(straightRegress) {
                  resR <- z - xr %*% omega
                  resRdiff <- y - xreg %*% omega - mmu
                }
                else {
                  resR <- funcTF(y = z, x = as.vector(xr), delta = delta, omega = omega, b = b,
                                 rx = r, sx = s, nx = length(z), meanval = 0)$y  ##regression residuals

                  resRdiff <- funcTF(y = y, x = as.vector(xreg), delta = delta, omega = omega,
                                     b = b, rx = r, sx = s, nx = length(y), meanval = 0)$y

                }

                res <- DLResiduals(rr, resRdiff)

                sigma2muse <- getsigma2muhat(resRdiff, phi = allpars[[i]]$phi, theta = allpars[[i]]$theta,
                                             thetaseas = allpars[[i]]$thetaseas, phiseas = allpars[[i]]$phiseas, period = period,
                                             dfrac = allpars[[i]]$dfrac, dfs = allpars[[i]]$dfs, H = allpars[[i]]$H,
                                             Hs = allpars[[i]]$Hs, alpha = allpars[[i]]$alpha, alphas = allpars[[i]]$alphas)
                allpars[[i]]$fitted <- resRdiff - res

            } else {

                yy <- y - allpars[[i]]$muHat
                resR <- NULL
                sigma2muse <- getsigma2muhat(yy, phi = allpars[[i]]$phi, theta = allpars[[i]]$theta,
                  thetaseas = allpars[[i]]$thetaseas, phiseas = allpars[[i]]$phiseas, period = period,
                  dfrac = allpars[[i]]$dfrac, dfs = allpars[[i]]$dfs, H = allpars[[i]]$H,
                  Hs = allpars[[i]]$Hs, alpha = allpars[[i]]$alpha, alphas = allpars[[i]]$alphas)
                res <- DLResiduals(rr, yy)
                allpars[[i]]$fitted <- yy - res
            }
            allpars[[i]]$regResiduals <- resR
            allpars[[i]]$residuals <- res
            allpars[[i]]$sigma2 <- sigma2muse

            class(allpars[[i]]) <- "ARFIMA"
        }

        ord <- order(logl, decreasing = TRUE)

        allpars <- allpars[ord]
        index <- 1:num

        weeded <- (num == 1)

        ans <- list(z = z, differencing = differencing, dint = dint, dseas = dseas, period = period,
            modes = allpars, dmean = dmean, itmean = itmean, p = p, q = q, pseas = pseas,
            qseas = qseas, lmodel = lmodel, slmodel = slmodel, weeded = weeded, getHess = getHess,
            numvars = numvars, numcut = 0, n = length(z), xreg = xr, r = r, s = s, b = b,
            call = match.call(), flag0 = flag0, numeach = numeach, strReg = straightRegress,
            regOnly = regOnly, intindex = intindex, intname = intname, namexreg=namexreg,
            fixx = fixx, indfixx = indfixx)
    }

    class(ans) <- "arfima"
    if (autoweed && !ans$weeded) {
        ans <- weed(ans, type = weedtype[1], eps2 = weedeps, adapt = adapt, pn = weedp)
    }
    ans
}

# Tells whether points share a circle in the p-norm
pcircle <- function(x1, x2, rad, pn = 2) {
    zz <- x1 - x2
    z <- (sum(abs(zz)^pn))^(1/pn)
    return(z > rad)
}


# Weeds out similar (or the same) modes


#' Weeds out fits from a call to arfima that are too close to each other.
#'
#' Weeds out fits from a call to arfima that are too close to each other.
#'
#'
#' @param ans The result from a call to arfima.
#' @param type The space to perform the weeding in. "A" is for operating
#' parameters. "P" is in the PACF space. "B" performs weeding in both. "N"
#' performs no weeding and is only used internally.
#' @param walls If more than one mode is on a wall in the PACF space, all modes
#' but the one with the highest log-likelihood on the same wall are deleted.
#' @param eps2 The maximum distance between modes that are close together for
#' the mode with the lower log-likelihood to be weeded out.  If adapt is TRUE
#' (default) this value changes.
#' @param eps3 The minimum distance from a wall for a secondary mode to be
#' weeded out, if walls are TRUE.
#' @param adapt If TRUE, if dim is the dimensionality of the search, eps2 is
#' changed to \eqn{(1 + eps2)^{dim} - 1}.
#' @param pn The p in the p-norm to be used in the weeding.  p = 2 (default) is
#' Euclidean distance.
#' @return An object of class "arfima" with modes possibly weeded out.
#' @author JQ (Justin) Veenstra
#' @seealso \code{\link{arfima}}, \code{\link{distance}}
#' @keywords ts
#' @examples
#'
#' \donttest{
#' set.seed(1234)
#' sim <- arfima.sim(1000, model = list(theta = 0.9, dfrac = 0.4))
#' fit <- arfima(sim, order = c(0, 0, 1), autoweed = FALSE, back=T)
#' fit
#' distance(fit)
#' fit1 <- weed(fit)
#' fit1
#' distance(fit1)
#' }
#'
#' @export weed
weed <- function(ans, type = c("A", "P", "B", "N"), walls = FALSE, eps2 = 0.025, eps3 = 0.01,
    adapt = TRUE, pn = 2) {
    if (class(ans) != "arfima")
        stop("weed only defined for arfima objects")

    allpars <- ans$modes

    if(length(allpars)==1) {
      ans$weeded <- TRUE
      return(ans)
    }

    if (ans$flag0){
      stopifnot(length(allpars[[1]]$pars)==0)
      allpars <- allpars[[1]]
      ans$modes <- allpars
      ans$weeded <- TRUE
      return(ans)
    }

    if (length(pn) == 0 || pn <= 0 || !is.finite(pn)) {
        warning("invalid pn in weed: setting to 2")
        pn <- 2
    }
    if (length(eps2) == 0 || eps2 < 0) {
        warning("invalid epsilon in weed: setting to 0.025")
        eps2 <- 0.025
    }

    type <- type[1]

    if (!(type %in% c("A", "P", "B", "N"))) {
        warning("invalid type in weed: setting to \"A\"")
        type <- "A"
    }
    if (length(eps3) == 0 || eps3 <= 0) {
        warning("invalid wall epsilon in weed: setting to 0.01")
        eps3 <- 0.01
    }
    if (!is.logical(walls)) {
        warning("invalid walls selection in weed: setting to FALSE")
        walls <- FALSE
    }

    modes <- vector("list", 1)
    p <- ans$p
    q <- ans$q
    pseas <- ans$pseas
    qseas <- ans$qseas
    lmodel <- ans$lmodel
    slmodel <- ans$slmodel

    numvars <- p + q + qseas + pseas + if (lmodel != "n")
        1 else 0
    numvars <- numvars + if (slmodel != "n")
        1 else 0
    numvars <- numvars + sum(ans$r) + sum(ans$s)

    if(numvars==0) {
      stop("Should not get here.  Please contact the package maintainer.")
    }

    if (adapt)
        eps2 <- (1 + eps2)^numvars - 1

    num <- length(allpars)
    logl <- rep(0, num)

    for (i in 1:num) logl[i] <- allpars[[i]]$loglik

    ord <- order(logl, decreasing = TRUE)

    allpars <- allpars[ord]


    if (type %in% c("A", "B")) {
        modes <- vector("list", 1)
        modes[[1]] <- allpars[[1]]
        flag <- NULL
        if (num > 1)
            for (i in 2:num) {
                flag <- NULL
                for (j in 1:length(modes)) {
                  flag <- c(flag, pcircle(x1 = allpars[[i]]$pars, x2 = modes[[j]]$pars,
                    rad = eps2, pn = pn))
                }
                if (all(flag)) {
                  modes[[length(modes) + 1]] <- allpars[[i]]
                }
            }

    } else modes <- allpars

    num <- length(modes)

    logl <- rep(0, num)

    for (i in 1:num) logl[i] <- modes[[i]]$loglik

    ord <- order(logl, decreasing = TRUE)

    modes <- modes[ord]

    if (type %in% c("P", "B")) {
        modes1 <- modes
        modes <- vector("list", 1)
        modes[[1]] <- modes1[[1]]

        flag <- NULL
        if (!ans$regOnly && num > 1) {
            for (i in 2:num) {
                flag <- NULL
                for (j in 1:length(modes)) {
                  flag <- c(flag, pcircle(x1 = modes1[[i]]$parsp, x2 = modes[[j]]$parsp,
                    rad = eps2, pn = pn))
                }
                if (all(flag)) {
                  modes[[length(modes) + 1]] <- modes1[[i]]
                }
            }
        }
        else if(num>1) {
          for (i in 2:num) {
            flag <- NULL
            for (j in 1:length(modes)) {
              flag <- c(flag, pcircle(x1 = modes1[[i]]$omega, x2 = modes[[j]]$omega,
                                      rad = eps2, pn = pn))
            }
            if (all(flag)) {
              modes[[length(modes) + 1]] <- modes1[[i]]
            }
          }
        }
    }

    wall <- inds <- vector("list", numvars * 2)
    pl1 <- p2  <- 0
    arma_varnum <- p + q + pseas + qseas
    for (i in 1:numvars) {
      if(ans$regOnly) break
      if (i <= arma_varnum) {
        minner = -1
        maxxer = 1
      } else if (lmodel != "n" && !pl1) {
        pl1 <- 1
        if (lmodel == "d") {
          minner = -1
          maxxer = 0.5
        } else if (lmodel == "g") {
          minner = 0
          maxxer = 1
        } else if (lmodel == "h") {
          minner = 0
          maxxer = 3
        }
      } else {
        if (slmodel == "d") {
          minner = -1
          maxxer = 0.5
        } else if (slmodel == "g") {
          minner = 0
          maxxer = 1
        } else if (slmodel == "h") {
          minner = 0
          maxxer = 3
        }
      }
      for (j in 1:length(modes)) {
        wall[[2 * i - 1]] <- c(wall[[2 * i - 1]], abs(maxxer - modes[[j]]$parsp[i]) <
                                 eps3)
        wall[[2 * i]] <- c(wall[[2 * i]], abs(minner - modes[[j]]$parsp[i]) < eps3)
      }
      inds[[2 * i - 1]] <- which(wall[[2 * i - 1]])
      inds[[2 * i]] <- which(wall[[2 * i]])
    }

    index <- numeric(0)
    for (i in 1:(2 * numvars)) index <- c(index, inds[[i]])
    index <- unique(index)
    numcut <- length(ans$index) - length(index)
    ans1 <- list(z = ans$z, differencing = ans$differencing, dint = ans$dint, dseas = ans$dseas,
        period = ans$period, modes = modes, dmean = ans$dmean, itmean = ans$itmean, getHess = ans$getHess,
        p = p, q = q, pseas = pseas, qseas = qseas, lmodel = lmodel, slmodel = slmodel,
        numcut = numcut, index = index, weeded = TRUE, n = ans$n, xreg = ans$xreg, r = ans$r,
        s = ans$s, b = ans$b, call = ans$call, flag0 = ans$flag0, numeach = ans$numeach,
        strReg = ans$strReg, namexreg = ans$namexreg, regOnly = ans$regOnly, intindex = ans$intindex,
        intname = ans$intname, fixx = ans$fixx, indfixx = ans$indfixx)
    class(ans1) <- "arfima"
    ans1


}

# computes and reports the distance between all modes


#' The distance between modes of an \code{arfima} fit.
#'
#' The distance between modes of an \code{arfima} fit.
#'
#'
#' @param ans An object of class "arfima".
#' @param p The p in the p-norm to be used.
#' @param digits The number of digits to print.
#' @return A list of two data frames: one with distances in operator space, the
#' second with distances in the transformed (PACF) space.
#' @author JQ (Justin) Veensta
#' @references Veenstra, J.Q. Persistence and Antipersistence:  Theory and
#' Software (PhD Thesis)
#' @keywords ts
#' @examples
#'
#' \donttest{
#' set.seed(8564)
#' sim <- arfima.sim(1000, model = list(phi = c(0.2, 0.1), dfrac = 0.4, theta = 0.9))
#' fit <- arfima(sim, order = c(2, 0, 1), back=T)
#'
#' fit
#'
#' distance(fit)
#' }
#'
#' @export distance
distance <- function(ans, p = 2, digits = 4) {
    if (class(ans) != "arfima")
        stop("distance only defined for arfima objects")
    if (p <= 0)
        stop("invalid p")
    num <- length(ans$modes)
    distances1 <- distances2 <- matrix(0, nrow = num, ncol = num)
    if (num == 1)
        stop("only one mode")
    for (i in 1:(num)) {
        for (j in 1:num) {
            vals <- ans$modes[[i]]$pars - ans$modes[[j]]$pars
            distances1[i, j] <- (sum(abs(vals)^p))^(1/p)
            vals <- ans$modes[[i]]$parsp - ans$modes[[j]]$parsp
            distances2[i, j] <- (sum(abs(vals)^p))^(1/p)
        }
    }
    distances1 <- round(distances1, digits = digits)
    distances2 <- round(distances2, digits = digits)
    rownames(distances1) <- colnames(distances1) <- rownames(distances2) <- colnames(distances2) <- paste("Mode",
        1:num)
    distances <- vector("list", 2)
    distances[[1]] <- distances1
    distances[[2]] <- distances2
    names(distances) <- c(paste(p, "-norm distance", sep = ""), paste(p, "-norm  distance in transformed space",
        sep = ""))
    distances
}





#' Finds the best modes of an \code{arfima} fit.
#'
#' Finds the best modes of an \code{arfima} fit with respect to log-likelihood.
#'
#' This is the easiest way to remove modes with lower log-likelihoods.
#'
#' @param object An object of class "arfima".
#' @param bestn The top number of modes to keep with respect to the
#' log-likelihood.
#' @return The \code{bestn} "best" modes.
#' @author JQ (Justin) Veenstra
#' @seealso \code{\link{arfima}}
#' @keywords ts
#' @examples
#'
#' \donttest{
#' set.seed(8765)
#' sim <- arfima.sim(1000, model = list(phi = 0.4, theta = 0.9, dfrac = 0.4))
#' fit <- arfima(sim, order = c(1, 0, 1), back=T)
#' fit
#' fit <- bestModes(fit, 2)
#' fit
#' }
#'
#' @export bestModes
bestModes <- function(object, bestn) {

    if (class(object) != "arfima")
        stop("bestModes only defined for arfima objects")
    if (!object$weeded)
        stop("Please weed the object first before calling bestModes.")
    m <- length(object$modes)
    if (bestn > m)
        stop("bestn is larger than the number of modes")

    modes <- object$modes
    logls <- sapply(1:m, function(i) modes[[i]]$logl)

    ord <- order(logls, decreasing = TRUE)

    ord <- ord[1:bestn]

    modes <- modes[ord]
    object$modes <- modes
    object <- weed(object, eps2 = 0)
    object
}
