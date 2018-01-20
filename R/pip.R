#' Compute posterior inclusion probabilities (PIPs)
#'
#' Computes posterior probabilities that a feature is a true member of an assigned cluster
#'
#' @param pvalue a vector of p-values.
#' @param group a vector of group indicators (optional).
#' @param ... optional arguments to control a local FDR estimation.
#'
#' @return \code{pip} returns a vector of posterior inclusion probabilities
#'
#' @export pip
#'
#' @importFrom qvalue lfdr
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
pip <- function(pvalue, group = NULL, 
    ...) {
    prob <- 1 - lfdr(pvalue, ...)
    
    if (!is.null(group)) {
        k <- length(unique(group))
        prob <- vector("numeric", 
            length = length(pvalue))
        for (i in 1:k) {
            prob[group == i] <- 1 - 
                lfdr(pvalue[group == 
                  i], ...)
        }
    }
    
    return(prob)
}
