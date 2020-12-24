#' Compute posterior inclusion probabilities (PIPs)
#'
#' Computes posterior probabilities that a feature is a true member of an assigned cluster
#'
#' @param pvalue a vector of p-values.
#' @param group a vector of group indicators (optional).
#' @param pi0 a vector of pi0 values (optional)
#' @param ... optional arguments for \code{\link[qvalue]{lfdr}} to control a local FDR estimation.
#'
#' @return \code{pip} returns a vector of posterior inclusion probabilities
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' 
#' @export
pip <- function(pvalue, group = NULL, pi0 = NULL,
    ...) {
    if(is.null(pi0)) {
        message("Using qvalue::pi0est to estimate pi0 values.")
        if(is.null(group)) {
            prob <- 1 - qvalue::lfdr(pvalue, ...)
        } else {
            k <- length(unique(group))
            prob <- vector("numeric",
                length = length(pvalue))
            for (i in 1:k) {
                prob[group == i] <- 1 -
                    qvalue::lfdr(pvalue[group == i], ...)
            }
        }
    } else {
        message(paste0("Using pi0 values:", pi0))
        if(is.null(group)) {
            prob <- 1 - qvalue::lfdr(pvalue, ...)
        } else {
            k <- length(unique(group))
            if(length(pi0) == 1) {
                pi0 <- rep(pi0, k)
            } else if(length(pi0) != k) {
                stop("When providing pi0 values, the length of pi0 has to be either 1 or equal the number of groups.")
            }
            prob <- vector("numeric",
                           length = length(pvalue))
            for (i in 1:k) {
                prob[group == i] <- 1 -
                    qvalue::lfdr(pvalue[group == i], pi0=pi0[i], ...)
            }
        }
    }

    return(prob)
}
