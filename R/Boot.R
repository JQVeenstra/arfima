#' Generic Bootstrap Function
#' 
#' Generic function to bootstrap a fitted model.
#' 
#' At present, the only function implemented is \code{\link{Boot.arfima}}.
#' 
#' @param obj fitted object
#' @param R number of bootstrap replicates
#' @param ... optional arguments
#' @return Parametric bootstrap simulation
#' @author A.I. McLeod and Justin Veenstra
#' @seealso \code{\link{Boot.arfima}}
#' @keywords ts
#' @export Boot
Boot <- function(obj, R = 1, ...) {
    UseMethod("Boot")
} 
