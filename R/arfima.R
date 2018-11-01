

#' Simulates, fits, and predicts persistent and anti-persistent time series.
#' arfima
#'
#' Simulates with arfima.sim, fits with arfima, and predicts with a method for
#' the generic function.  Plots predictions and the original time series. Has
#' the capability to fit regressions with ARFIMA/ARIMA-FGN/ARIMA-PLA errors, as
#' well as transfer functions/dynamic regression.
#'
#' \tabular{ll}{ Package: \tab arfima\cr Type: \tab Package\cr Version: \tab
#' 1.4-0\cr Date: \tab 2017-06-20\cr License: \tab MIT \cr }
#'
#' A list of functions:
#'
#' \code{\link{arfima.sim}} - Simulates an ARFIMA, ARIMA-FGN, or ARIMA-PLA
#' (three classes of mixed ARIMA hyperbolic decay processes) process, with
#' possible seasonal components.
#'
#' \code{\link{arfima}} - Fits an ARIMA-HD (default single-start) model to a series,
#' with options for regression with ARIMA-HD errors and dynamic regression
#' (transfer functions).  Allows for fixed parameters as well as choices for
#' the optimizer to be used.
#'
#' \code{\link{arfima0}} - Simplified version of \code{arfima}
#'
#' \code{\link{weed}} - Weeds out modes too close to each other in the same
#' fit.  The modes with the highest log-likelihoods are kept
#'
#' \code{\link{print.arfima}} - Prints the relevant output of an \code{arfima}
#' fitted object, such as parameter estimates, standard errors, etc.
#'
#' \code{\link{summary.arfima}} - A much more detailed version of
#' \code{print.arfima}
#'
#' \code{\link{coef.arfima}} - Extracts the coefficients from a \code{arfima}
#' object
#'
#' \code{\link{vcov.arfima}} - Theoretical and observed covariance matrices of
#' the coefficients
#'
#' \code{\link{residuals.arfima}} - Extracts the residuals or regression
#' residuals from a \code{arfima} object
#'
#' \code{\link{fitted.arfima}} - Extracts the fitted values from a
#' \code{arfima} object
#'
#'
#' \code{\link{tacvfARFIMA}} - Computes the theoretical autocovariance function
#' of a supplied model.  The model is checked for stationarity and
#' invertibility.
#'
#' \code{\link{iARFIMA}} - Computes the Fisher information matrix of all
#' non-FGN components of the given model.  Can be computed (almost) exactly or
#' through a psi-weights approximation.  The approximation takes more time.
#'
#' \code{\link{IdentInvertQ}} - Checks whether the model is identifiable,
#' stationary, and invertible.  Identifiability is checked through the
#' information matrix of all non-FGN components, as well as whether both types
#' of fractional noise are present, both seasonally and non-seasonally.
#'
#' \code{\link{lARFIMA}} and \code{\link{lARFIMAwTF}} - Computes the
#' log-likelihood of a given model with a given series.  The second admits
#' transfer function data.
#'
#' \code{\link{predict.arfima}} - Predicts from an \code{arfima} object.
#' Capable of exact minimum mean squared error predictions even with integer d
#' > 0 and/or integer dseas > 0. Does not include transfer function/leading
#' indicators as of yet.  Returns a \code{predarfima} object, which is composed
#' of: predictions, and standard errors (exact and, if possible, limiting).
#'
#' \code{\link{print.predarfima}} - Prints the relevant output from a
#' \code{predarfima} object: the predictions and their standard deviations.
#'
#' \code{\link{plot.predarfima}} - Plots a \code{predarfima} object.  This
#' includes the original time series, the forecasts and as default the
#'  standard 95\% prediction intervals (exact and, if available, limiting).
#'
#' \code{\link{logLik.arfima}}, \code{\link{AIC.arfima}},
#' \code{\link{BIC.arfima}} - Extracts the requested values from an
#' \code{\link{arfima}} object
#'
#' \code{\link{distance}} - Calculates the distances between the modes
#'
#' \code{\link{removeMode}} - Removes a mode from a fit
#'
#' \code{\link{tacvf}} - Calculates the theoretical autocovariance functions
#' (tacvfs) from a fitted \code{arfima} object
#'
#' \code{\link{plot.tacvf}} - Plots the tacvfs
#'
#' \code{\link{print.tacvf}} - Prints the tacvfs
#'
#' \code{\link{tacfplot}} - Plots the theoretical autocorrelation functions
#' (tacfs) of different models on the same data
#'
#' \code{\link{SeriesJ}}, \code{\link{tmpyr}} - Two datasets included with the
#' package
#'
#' @name arfima-package
#' @docType package
#' @author JQ (Justin) Veenstra, A. I. McLeod
#'
#' Maintainer: JQ (Justin) Veenstra <jqveenstra@gmail.com>
#' @references Veenstra, J.Q. Persistence and Antipersistence:  Theory and
#' Software (PhD Thesis)
#' @keywords package
#' @examples
#'
#' \donttest{
#' set.seed(8564)
#' sim <- arfima.sim(1000, model = list(phi = c(0.2, 0.1), dfrac = 0.4, theta = 0.9))
#' fit <- arfima(sim, order = c(2, 0, 1), back=TRUE)
#'
#' fit
#'
#' data(tmpyr)
#'
#' fit1 <- arfima(tmpyr, order = c(1, 0, 1), numeach = c(3, 3), dmean = FALSE)
#' fit1
#'
#' plot(tacvf(fit1), maxlag = 30, tacf = TRUE)
#'
#' fit2 <- arfima(tmpyr, order = c(1, 0, 0), numeach = c(3, 3), autoweed = FALSE,
#' dmean = FALSE)
#'
#' fit2
#'
#' fit2 <- weed(fit2)
#'
#' fit2
#'
#' tacfplot(fits = list(fit1, fit2))
#'
#' fit3 <- removeMode(fit2, 2)
#'
#' fit3
#'
#' coef(fit2)
#' vcov(fit2)
#'
#' fit1fgn <- arfima(tmpyr, order = c(1, 0, 1), numeach = c(3, 3),
#' dmean = FALSE, lmodel = "g")
#' fit1fgn
#'
#' fit1hd <- arfima(tmpyr, order = c(1, 0, 1), numeach = c(3, 3),
#' dmean = FALSE, lmodel = "h")
#' fit1hd
#'
#' data(SeriesJ)
#' attach(SeriesJ)
#'
#' fitTF <- arfima(YJ, order= c(2, 0, 0), xreg = XJ, reglist =
#' list(regpar = c(1, 2, 3)), lmodel = "n", dmean = FALSE)
#' fitTF
#'
#' detach(SeriesJ)
#'
#' set.seed(4567)
#'
#' sim <- arfima.sim(1000, model = list(phi = 0.3, dfrac = 0.4, dint = 1),
#' sigma2 = 9)
#'
#' X <- matrix(rnorm(2000), ncol = 2)
#'
#' simreg <- sim + crossprod(t(X), c(2, 3))
#'
#' fitreg <- arfima(simreg, order = c(1, 1, 0), xreg = X)
#'
#' fitreg
#'
#' plot(sim)
#'
#' lines(residuals(fitreg, reg = TRUE)[[1]], col = "blue")
#' ##pretty much a perfect match.
#' }
#'
NULL




#' Series J, Gas Furnace Data
#'
#' Gas furnace data, sampling interval 9 seconds; observations for 296 pairs of
#' data points.
#'
#' XJ is input gas rate in cubic feet per minute, YJ is percentage carbon
#' dioxide (CO2) in outlet gas.  X is the regressor.
#'
#' Box, Jenkins, and Reinsel (2008) fit an AR(2) to YJ, with transfer function
#' specifications r = 2, s = 2, and b = 3, regressing on XJ.  Our package
#' agrees with their results.
#'
#' @name SeriesJ
#' @docType data
#' @format List with ts objects XJ and YJ.
#' @references Box, G. E. P., Jenkins, G. M., and Reinsel, G. C. (2008) Time
#' Series Analysis: Forecasting and Control.  4th Edition. John Wiley and Sons,
#' Inc., New Jersey.
#'
#' Veenstra, J. and McLeod, A. I. (Working Paper). The arfima R package: Exact
#' Methods for Hyperbolic Decay Time Series
#' @source Box, Jenkins and Reinsel(2008).  Time Series Analysis: Forecasting
#' and Control.
#' @keywords datasets
#' @examples
#'
#' data(SeriesJ)
#' attach(SeriesJ)
#'
#' fitTF <- arfima(YJ, order= c(2, 0, 0), xreg = XJ, reglist =
#' list(regpar = c(2, 2, 3)), lmodel = "n")
#' fitTF ## agrees fairly closely with Box et. al.
#'
#'
#' detach(SeriesJ)
#'
NULL





#' Temperature Data
#'
#' Central England mean yearly temperatures from 1659 to 1976
#'
#' Hosking notes that while the ARFIMA(1, d, 1) has a lower AIC, it is not much
#' lower than the AIC of the ARFIMA(1, d, 0).
#'
#' Bhansali and Kobozka find: muHat = 9.14, d = 0.28, phi = -0.77, and theta =
#' -0.66 for the ARFIMA(1, d, 1), which is close to our result, although our
#' result reveals trimodality if \code{numeach} is large enough.  The third
#' mode is close to Hosking's fit of an ARMA(1, 1) to these data, while the
#' second is very antipersistent.
#'
#' Our package gives a very close result to Hosking for the ARFIMA(1, d, 0)
#' case, although there is also a second mode.  Given how close it is to the
#' boundary, it may or may not be spurious.  A check with \code{dmean = FALSE}
#' shows that it is not the optimized mean giving a spurious mode.
#'
#' If, however, we use \code{whichopt = 1}, we only have one mode.  Note that
#' Nelder-Mead sometimes does take out non-spurious modes, or add spurious
#' modes to the surface.
#'
#' @name tmpyr
#' @docType data
#' @format A ts tmpyr
#' @references Parker, D.E., Legg, T.P., and Folland, C.K. (1992).  A new daily
#' Central England Temperature Series, 1772-1991. Int. J. Clim., Vol 12, pp
#' 317-342
#'
#' Manley,G. (1974).  Central England Temperatures: monthly means 1659 to 1973.
#' Q.J.R. Meteorol. Soc., Vol 100, pp 389-405.
#'
#' Hosking, J. R. M. (1984). Modeling persistence in hydrological time series
#' using fractional differencing, Water Resour. Res., 20(12)
#'
#' Bhansali, R. J. and Koboszka, P. S. (2003) Prediction of Long-Memory Time
#' Series In Doukhan, P., Oppenheim, G. and Taqqu, M. S. (Eds) Theory and
#' Applications of Long-Range Dependence (pp355-368) Birkhauser Boston Inc.
#'
#' Veenstra, J.Q. Persistence and Antipersistence:  Theory and
#' Software (PhD Thesis)
#'
#' @source \url{http://www.metoffice.gov.uk/hadobs/hadcet/}
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(tmpyr)
#'
#' fit <- arfima(tmpyr, order = c(1, 0, 1), numeach = c(3, 3), dmean = TRUE, back=TRUE)
#' fit
#' ##suspect that fourth mode may be spurious, even though not close to a boundary
#' ##may be an induced mode from the optimization of the mean
#'
#' fit <- arfima(tmpyr, order = c(1, 0, 1), numeach = c(3, 3), dmean = FALSE, back=TRUE)
#' fit
#'
#' ##perhaps so
#'
#'
#' plot(tacvf(fit), maxlag = 30, tacf = TRUE)
#'
#' fit1 <- arfima(tmpyr, order = c(1, 0, 0), dmean = TRUE, back=TRUE)
#' fit1
#'
#' fit2 <- arfima(tmpyr, order = c(1, 0, 0), dmean = FALSE, back=TRUE)
#' fit2  ##still bimodal.  Second mode may or may not be spurious.
#'
#' fit3 <- arfima(tmpyr, order = c(1, 0, 0), dmean = FALSE, whichopt = 1, numeach = c(3, 3))
#' fit3  ##Unimodal.  So the second mode was likely spurious.
#'
#' plot(tacvf(fit2), maxlag = 30, tacf = TRUE)
#' ##maybe not spurious.  Hard to tell without visualizing the surface.
#'
#' ##compare to plotted tacf of fit1:  looks alike
#' plot(tacvf(fit1), maxlag = 30, tacf = TRUE)
#'
#' tacfplot(list(fit1, fit2))
#' }
#'
NULL
