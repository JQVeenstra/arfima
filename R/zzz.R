arfima_const_predvarnum <- 5 #7 with bootstrap

message = "Note that the arfima package has new defaults starting with
1.4-0: type arfimachanges() for a list, as well as some other notes.
NOTE: some of these are quite important!"


.onAttach <- function(libname, pkgname) {
  packageStartupMessage(message)
}


changes1.7.0 = "
Changes in arfima starting in 1.4-0:
  1. arfima() now defaults to searching for only one mode:
    previously it had, as default, for p time series parameters,
    2^p fits attempted to find multiple modes.  You can go back to
    this by setting numeach = c(2, 2) and seasonal$numeach = c(2, 2)
    should you prefer.  Generally speaking, numeach = c(r, m) starts
    r fits for each ARMA (resp. SARMA) parameters and m fits for each
    fractional parameter, for a total of (rm)^p fits (in the nonseasonal
    case).
    Note that in the next major version of the package, there will be more
    flexibility in multiple mode finding.
  2. For backwards compatibility for the above two points, set back = T
    in the call to arfima.  That is, setting back = T will set numeach = c(2, 2)
    seasonal$numeach = c(2, 2).
  3. For your convenience, you can also pass a matrix of starting values instead
    of setting numeach (or the random start option).
  4. predict.arfima() (that is, calling predict on an arfima object) now
    has parameters prop.use (default = 0.25) and modes (default = 'all').
    The latter allows you to select which modes you would like to forecast
    from.
    The former speeds up forecasting by only using the last part of the
    model, with minimum use of 1000 points.  This minimum can be overridden
    by using the min.prop argument, which defaults to NULL, and must be an
    integer >=1 if set.
  5. predict.arfima now uses the median of the bootstrapped predictions, and
    the name of the value has been changed from meanvalnp to medvalnp.
      NOTE: bootstrap predictions are not a part of the package from 1.5-0
        until further notice.  There is some fault with them that needs to
        be worked out.

Changes in arfima starting in 1.5-0:
  1.  The Boot.arfima and Boot.ARFIMA functions have been temporarily removed
    from the package.  There is an issue I need to track down.
  2.  As a consequence of the above, there are no longer bootstrap predictions;
    this is actually the major issue I need to track down above.

Changes in arfima starting in 1.6-0:
  1.  Prediction with integrated series and xreg has been fixed.
Changes in arfima starting in 1.6-2:
  1.  Prediction with integrated series without xreg has been fixed.  This bug was
    created by the previous fix in 1.6-0 and 1.6-1.

Changes in arfima starting in 1.7-0:
  1.  Added the sim_from_fixed function 

Note: on prediction intervals with xreg, the intervals are only with respect to
the time series itself.  I will update to include standard errors on beta
as soon as I am able.

Finally, please note that in the next major version of the package, there will
be multiple other changes, some breaking backwards compatibility.
"
#' Prints changes to the package since the last update.  Started in 1.4-0
arfimachanges <- function()
  cat(changes1.7.0)


