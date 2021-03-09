#' Compute posterior inclusion probabilities (PIPs)
#'
#' Computes posterior probabilities that a feature is a true member of an assigned cluster
#'
#' @param pvalue a vector of p-values.
#' @param group a vector of group indicators (optional).
#' If provided, PIP analysis is stratified.
#' Assumes groups are in 1:k where k is the number of unique groups.
#' @param pi0 a vector of pi0 values (optional).
#' Its length has to be either 1 or equal the number of groups.
#' @param verbose If TRUE, reports information.
#' @param ... optional arguments for \code{\link[qvalue]{lfdr}} to control a local FDR estimation.
#'
#' @return \code{pip} returns a vector of posterior inclusion probabilities
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' 
#' @export
pip <- function(
                pvalue,
                group = NULL,
                pi0 = NULL,
                verbose = TRUE,
                ...
                ) {
    # pvalue is always mandatory
    if ( missing( pvalue ) )
        stop( '`pvalue` is required!' )

    m <- length(pvalue)
    
    if ( is.null( pi0 ) ) {
        if (verbose)
            message("Using qvalue::pi0est to estimate pi0 values.")
        if ( is.null( group ) ) {
            # pooled PIP
            prob <- 1 - qvalue::lfdr(pvalue, ...)
        } else {
            # stratified PIP
            if ( length( group ) != m )
                stop( '`group` length (', length(group), ') must match `pvalue` length (', m, ')' )
            k <- length(unique(group))
            prob <- vector("numeric", length = m )
            for (i in 1:k) {
                # skip empty groups quietly
                if ( sum(group == i) > 0 ) 
                    prob[group == i] <- 1 - qvalue::lfdr( pvalue[group == i], ... )
            }
        }
    } else {
        if (verbose)
            message("Using pi0 values:", pi0)
        if (is.null(group)) {
            # pooled PIP
            # NOTE: pi0 not used here
            prob <- 1 - qvalue::lfdr(pvalue, ...)
        } else {
            # stratified PIP
            if ( length( group ) != m )
                stop( '`group` length (', length(group), ') must match `pvalue` length (', m, ')' )
            k <- length(unique(group))
            if ( length( pi0 ) == 1 ) {
                pi0 <- rep(pi0, k)
            } else if (length(pi0) != k) {
                stop("When providing pi0 values, the length of pi0 has to be either 1 or equal the number of groups.")
            }
            prob <- vector("numeric", length = m )
            for (i in 1:k) {
                # skip empty groups quietly
                if ( sum(group == i) > 0 ) 
                    prob[group == i] <- 1 - qvalue::lfdr( pvalue[group == i], pi0 = pi0[i], ... )
            }
        }
    }

    return(prob)
}
