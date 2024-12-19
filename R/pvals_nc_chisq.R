#' Calculate high precision p-values from Jackstraw data under non-central chi squared distribution for null statistics
#'
#' This function takes a Jackstraw object for convenience, estimates the non-centrality parameter (NCP) from the null statistics (using [ncp_est()]), then uses this non-central chi squared model to calculate p-values for the observed statistics.
#' The goal is to be able to calculate small p-values with arbitrary precision, which is normally not possible with the ordinary Jackstraw functions since the minimum empirical p-value is determined by the inverse of the number of null samples (`1 /( s * B )` where `s` and `B` are Jackstraw parameters (see [jackstraw_lfa()])).
#' In contrast, large p-values values are typically similar between empirical and NCP model versions, and in fact this must hold if the null model is correctly specified.
#' This model works well empirically for [jackstraw_lfa()] and [jackstraw_alstructure()], whose statistics are deviances that would ordinarily be central chi-squared distributed were the test factors not estimated from the same data.
#' This model is not appropriate for data from [jackstraw_pca()], whose statistics have a roughly F distribution, and possibly many other functions from this package.
#'
#' @param out The object returned by [jackstraw_lfa()] and [jackstraw_alstructure()].
#' In particular, a list with at least the following two named elements: `null.stat` containing the null statistics on which the NCP null model is estimated, and `obs.stat` containing the observed or actual statistics for which p-values are calculated.
#' @param df The degrees of freedom of the test used to calculate the statistics.
#' Note that for ordinary jackstraw runs this is `r-1`, since the null model/intercept is counted among the `r` latent variables, but alternate values can be provided for more complex tests.
#' @param null.stat Alternate way to provide these values if the input is not a Jackstraw object, so `out` is not available.
#' @param obs.stat Alternate way to provide these values if the input is not a Jackstraw object, so `out` is not available.
#'
#' @return A list with two named values:
#' - `p.value`: the vector of p-values calculated for the input vector `obs.stat`
#' - `ncp`: the numeric vector of two values of NCP estimates calculated from `null.stat`, namely the maximum likelihood (MLE) and method-of-moments estimates (see [ncp_est()]).  Only MLE is used to calculate p-values.
#'
#' @examples
#' # instead of running jackstraw here, we simulate toy data
#' df <- 1
#' ncp_true <- 1
#' # null data are truly non-central chi squared distributed
#' # observed data are the same but with huge power reflected in a larger NCP
#' out <- list(
#'     null.stat = rchisq( 100, df, ncp_true ),
#'     obs.stat = rchisq( 11, df, ncp_true + 10 )
#' )
#'
#' # this calculates new p-values with much higher precision for highly significant cases
#' out2 <- pvals_nc_chisq( out, df )
#' # these are the desired p-values
#' out2$p.value
#' # and NCP estimates from two methods
#' out2$ncp
#'
#' \dontrun{
#' # This is the more typical usage in practice
#' # first run Jackstraw-LFA, which returns empirical p-values with limited precision and raw stats
#' out <- jackstraw_lfa( dat, r, FUN = function(x) lfa( x, r ) )
#' # then fit NCP model and use it to calculate high-precision p-values!
#' out2 <- pvals_nc_chisq( out, r-1 )
#'
#' # compare p-values!  In linear scale they should agree very well
#' plot( out$p.value, out2$p.value, pch = '.' )
#' # but in log scale it's clear the NCP values can take on much smaller values
#' plot( out$p.value, out2$p.value, pch = '.', log = 'xy' )
#' }
#'
#' @seealso
#' [ncp_est()] for more notes on estimating non-centrality parameters only.
#'
#' [jackstraw_lfa()], [jackstraw_alstructure()]
#'
#' @export
pvals_nc_chisq <- function( out, df, null.stat = out$null.stat, obs.stat = out$obs.stat ) {
    # determine the logic of presence of inputs, which is more complicated than usual
    if ( missing( out ) ) {
        if ( missing( null.stat ) )
            stop( '`null.stat` must be provided if `out` is missing!' )
        if ( missing( obs.stat ) )
            stop( '`obs.stat` must be provided if `out` is missing!' )
    }
    if ( missing( df ) )
        stop( '`df` is required!' )
    
    # estimate non-centrality params
    ncp_ests  <- ncp_est( null.stat, df )
    # and p-values according to this model
    pvals <- stats::pchisq( obs.stat, df, ncp_ests[1], lower.tail = FALSE )
    pvals[pvals==0] <- .Machine$double.xmin
    # return all the useful data, the first in the same style as the jackstraw_* functions
    return( list( p.value = pvals, ncp = ncp_ests ) )
}