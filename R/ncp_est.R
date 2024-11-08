#' Estimate parameter of non-central Chi squared distributed data
#'
#' This function estimates the non-centrality parameter (NCP) of data assuming a non-central chi squared distribution with known degrees of freedom.  The function first uses the method of moments (MOM) to get a rough estimate of the NCP, then uses that as a starting point for calculating the maximum likelihood estimate (MLE).
#'
#' There is one case where this function does not calculate MLEs directly, due to numerical issues, namely when the MOM estimate is negative, in which case estimates of zero for both are returned.
#' Negative MOM NCP estimates arise in two common cases.
#' First, the data has a central chi squared distribution, in which case the NCP estimate may take on small negative values, and overall these estimates have a mean of zero and variance dependent on sample size; in this case the return values of zero are appropriate.
#' Second, the data may not have a central or non-central chi squared distribution at all, in which case the estimation model is completely misspecified; we opted to return zero estimates quietly, without issuing errors or warnings, as this case overlaps with central chi-squared data and it's not trivial to tell the difference.
#' This function does not aim to detect misspecified data, it is up to the researcher to compare the model against the observed statistics for goodness of fit, for example using q-q plots.
#'
#' @param x The vector of data whose parameter we are interested in.
#' The data is assumed to follow a non-central chi squared distribution, and appropriate inputs are deviance or likelihood ratio test statistics, or other tests that ordinarily follow a central chi-squared distribution but become non-central due to "double dipping" (where fixed variables being tested were in fact estimated from the same data).
#' Any values that are NA, infinites, or negatives are ignored (removed before all calculations).
#' @param df The degrees of freedom of the test used to calculate the statistics.
#'
#' @return A numeric vector with two values, namely the MLE and MOM estimates.
#'
#' @examples
#' # true parameters of toy data
#' df <- 1
#' ncp_true <- 1
#' x <- rchisq( 100, df, ncp_true )
#'
#' # get estimates!
#' ncp_ests <- ncp_est( x, df )
#' # this is MLE
#' ncp_ests[1]
#' # this is MOM estimate
#' ncp_ests[2]
#'
#' @seealso
#' [pvals_nc_chisq()], a wrapper around this function that facilitates p-value calculation for Jackstraw objects.
#' 
#' @export
ncp_est <- function( x, df ) {
    # check inputs
    if ( missing( x ) )
        stop( '`x` is required!' )
    if ( missing( df ) )
        stop( '`df` is required!' )
    
    # remove NAs or my code complains
    if ( anyNA( x ) ) 
        x <- x[ !is.na( x ) ]
    # there can also be infinite values apparently...
    indexes <- which( !is.finite( x ) )
    if ( length( indexes ) > 0 )
        x <- x[ - indexes ]
    # lastly, negative values can occur when models are poorly fit, which also mess with the likelihood calculations, but they ought to be super rare so let's remove those too (don't set to zero, that also has undefined log-likelihood)
    indexes <- which( x < 0 )
    if ( length( indexes ) > 0 )
        x <- x[ - indexes ]
    
    # method of moments crude estimate, good starting point
    ncp_est_mom <- mean( x ) - df
    # this might look negative if data doesn't have the right distribution, but also if it is a central chisq (so this may be normal).
    # unfortunately this case behaves poorly with MLE code, better to just treat as zero
    if ( ncp_est_mom < 0 )
        return( c( 0, 0 ) )
    
    # now get MLE
    # first, negative log likelihood of data given parameters
    ll <- function( ncp )
        -sum( stats::dchisq( x, df, ncp, log = TRUE ) )

    # maximum likelihood fit of ncp
    model_fit <- stats4::mle(
                             ll,
                             start = list( ncp = ncp_est_mom ),
                             lower = list( ncp = 0 ),
                             method = "L-BFGS-B"
                         )
    ncp_est_mle <- stats4::coef( model_fit )[1]

    return( c( ncp_est_mle, ncp_est_mom ) )
}
