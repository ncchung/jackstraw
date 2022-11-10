#' Non-Parametric Jackstraw for Logistic Factor Analysis
#'
#' Test association between the observed variables and their latent variables captured by logistic factors (LFs).
#'
#' This function uses logistic factor analysis (LFA) from Hao et al. (2016).
#' Particularly, the deviance in logistic regression (the full model with \code{r} LFs vs. the intercept-only model) is used to assess significance.
#'
#' @param dat a genotype matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param r a number of significant LFs.
#' @param FUN a function to use for LFA (by default, it uses the \code{lfa} package)
#' @param r1 a numeric vector of LFs of interest (implying you are not interested in all \code{r} LFs).
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations. There will be a total of \code{s*B} null statistics.
#' @param covariate a data matrix of covariates with corresponding \code{n} observations (do not include an intercept term).
#' @param verbose a logical specifying to print the computational progress.
#'
#' @return \code{jackstraw_lfa} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their LFs}
#' \item{obs.stat}{\code{m} observed deviances}
#' \item{null.stat}{\code{s*B} null deviances}
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#'
#' @seealso  \link{jackstraw_pca} \link{jackstraw} \link{jackstraw_subspace}
#'
#' @examples
#' \dontrun{
#' ## simulate genotype data from a logistic factor model: drawing rbinom from logit(BL)
#' m <- 5000; n <- 100; pi0 <- .9
#' m0 <- round(m*pi0)
#' m1 <- m - round(m*pi0)
#' B <- matrix(0, nrow=m, ncol=1)
#' B[1:m1,] <- matrix(runif(m1*n, min=-.5, max=.5), nrow=m1, ncol=n)
#' L <- matrix(rnorm(n), nrow=1, ncol=n)
#' BL <- B %*% L
#' prob <- exp(BL)/(1+exp(BL))
#'
#' dat <- matrix(rbinom(m*n, 2, as.numeric(prob)), m, n)
#'
#' ## apply the jackstraw_lfa
#' out <- jackstraw_lfa(dat, r = 2)
#' }
#'
#' @export
jackstraw_lfa <- function(
                          dat,
                          r,
                          FUN = function(x) lfa::lfa(x, r),
                          r1 = NULL,
                          s = NULL,
                          B = NULL,
                          covariate = NULL,
                          verbose = TRUE
                          ) {
    # check mandatory data
    if ( missing( dat ) )
        stop( '`dat` is required!' )
    if ( missing( r ) )
        stop( '`r` is required!' )
    if ( !is.matrix( dat ) )
        stop( '`dat` must be a matrix!' )

    if ( !is.null( FUN ) ) {
        FUN <- match.fun( FUN )
    } else {
        stop("Please provide a function to estimate latent variables.")
    }

    # make sure `dat` is a regular matrix
    # also get correct dimensions
    if ( is.matrix( dat ) ) {
        # regular case
        m <- nrow( dat )
        n <- ncol( dat )
    } else
        stop( '`dat` must be a matrix!  Observed class: ', class( dat ) )

    # more validations of mandatory parameters
    if ( !(r > 0 && r < n) )
        stop( "`r` is not in valid range between `1` and `n-1` (`n` is number of individuals)." )

    # if there are covariates, the dimensions must agree
    # covariate can be either a vector or a matrix, test both cases
    if ( !is.null( covariate ) ) {
        if ( is.matrix( covariate ) ) {
            if ( nrow( covariate ) != n )
                stop( 'Matrix `covariate` must have `n` rows, has: ', nrow( covariate ), ', expected: ', n )
        } else {
            if ( length( covariate ) != n )
                stop( 'Vector `covariate` must have `n` elements, has: ', length( covariate ), ', expected: ', n )
        }
    }

    if (is.null(s)) {
        s <- round(m/10)
        if (verbose)
            message( "Number of null variables (s) to be permuted is not specified: s=round(0.10*m)=", s, "." )
    }
    if (is.null(B)) {
        B <- round(m * 10/s)
        if (verbose)
            message( "Number of resampling iterations (B) is not specified: B=round(m*10/s)=", B, "." )
    }

    if ( is.null( r1 ) )
        r1 <- 1:r
    if ( all( 1:r %in% r1 ) ) {
        # no adjustment LVs
        r0 <- NULL
        LFr0 <- NULL
    } else {
        r0 <- ( 1 : r )[ -r1 ]
    }

    ## LFr[,r] is an intercept term
    LFr <- FUN( dat )
    # check dimensions for this matrix
    if (ncol(LFr) != r)
        stop("The number of latent variables ", r, " is not equal to the number of column(s) provided by `FUN`")
    if (nrow(LFr) != n)
        stop("The number of individuals ", n, " is not equal to the number of rows provided by `FUN`")
    # suubset further if `r1` is non-trivial
    LFr1 <- LFr[ , r1, drop = FALSE ]
    # create null LFs if `r1` is non-trivial (otherwise LFr0 is NULL)
    if (!is.null(r0))
        LFr0 <- LFr[ , r0, drop = FALSE ]

    # both LF models get the same covariates appended (as they should be)
    # the resulting matrices are never `NULL` (that is not handled by this `gcatest` function).
    # there is an awkward ambiguity here, that LFr0 is assumed to not contain the intercept (so it is added), but regression will fail if that assumption is wrong!
    obs <- delta_deviance_lf(
                        X = dat,
                        LF0 = cbind(LFr0, matrix(1, n, 1), covariate),
                        LF1 = cbind(LFr, covariate)
                    )

    # Estimate null association
    # statistics
    null <- matrix(0, nrow = s, ncol = B)
    LFr0.js <- NULL

    if ( verbose )
        cat( paste0( "\nComputating null statistics (", B, " total iterations): " ) )
    for (i in 1:B) {
        # select the random subset of rows/loci to edit
        random_s <- sample.int( m, s )
        # extract that data and permute it, returning a matrix
        s.nulls <- dat[ random_s, , drop = FALSE ]
        s.nulls <- t( apply( s.nulls, 1, sample ) )
        # make a copy of this whole data containing the random subset
        jackstraw.dat <- dat
        jackstraw.dat[random_s, ] <- s.nulls

        # get LFs for this random data
        LFr.js <- FUN(jackstraw.dat)
        # dimensions are not rechecked here, it's just assumed that they are correct
        LFr1.js <- LFr.js[, r1, drop = FALSE]
        if (!is.null(r0))
            LFr0.js <- LFr.js[, r0, drop = FALSE]

        null[, i] <- delta_deviance_lf(
                            X = s.nulls,
                            LF0 = cbind(LFr0.js, matrix(1, n, 1), covariate),
                            LF1 = cbind(LFr.js, covariate)
                        )

        if ( verbose )
            cat(paste(i, " "))
    }

    p.value <- empPvals( obs, null )
    
    return(
        list(
            call = match.call(),
            p.value = p.value,
            obs.stat = obs,
            null.stat = null
        )
    )
}
