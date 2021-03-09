#' Non-Parametric Jackstraw for Logistic Factor Analysis
#'
#' Test association between the observed variables and their latent variables captured by logistic factors (LFs).
#'
#' This function uses logistic factor analysis (LFA) from Wei et al. (2014). Particularly, a deviation \code{dev} in logistic regression
#' (the full model with \code{r} LFs vs. the intercept-only model) is used to assess association.
#'
#' @param dat a genotype matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param r a number of significant LFs.
#' @param FUN a function to use for LFA (by default, it uses the lfa package)
#' @param devR use a R function to compute deviance. By default, FALSE (uses C++).
#' @param r1 a numeric vector of LFs of interest (implying you are not interested in all \code{r} LFs).
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations. There will be a total of \code{s*B} null statistics.
#' @param covariate a data matrix of covariates with corresponding \code{n} observations (do not include an intercept term).
#' @param verbose a logical specifying to print the computational progress.
#' @param seed a seed for the random number generator.
#'
#' @return \code{jackstraw_lfa} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their LFs}
#' \item{obs.stat}{\code{m} observed devs}
#' \item{null.stat}{\code{s*B} null devs}
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
#' B[1:m1,] <- matrix(runif(m1*n, min=-.5, max=.5), nrow=m1, ncol=1)
#' L <- matrix(rnorm(n), nrow=1, ncol=n)
#' BL <- B %*% L
#' prob <- exp(BL)/(1+exp(BL))
#'
#' dat <- matrix(rbinom(m*n, 2, as.numeric(prob)), m, n)
#'
#' ## apply the jackstraw_lfa
#' out <- jackstraw_lfa(dat, r = 2)
#'
#' ## apply the jackstraw_lfa using self-contained R functions
#' out <- jackstraw_lfa(dat, r = 2, FUN = function(x) lfa.corpcor(x, 2)[, , drop = FALSE], devR = TRUE)
#' }
#' 
#' @export
jackstraw_lfa <- function(
                          dat,
                          r,
                          FUN = function(x) lfa::lfa(x, r),
                          devR = FALSE,
                          r1 = NULL,
                          s = NULL,
                          B = NULL,
                          covariate = NULL,
                          verbose = TRUE,
                          seed = NULL
                          ) {
    # check mandatory data
    if ( missing( dat ) )
        stop( '`dat` is required!' )
    if ( missing( r ) )
        stop( '`r` is required!' )
    if ( !is.matrix( dat ) )
        stop( '`dat` must be a matrix!' )

    # more validations of mandatory parameters
    m <- nrow(dat)
    n <- ncol(dat)
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

    if (!is.null(seed))
        set.seed(seed)
    
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

    if ( !is.null( FUN ) ) {
        FUN <- match.fun( FUN )
    } else {
        stop("Please provide a function to estimate latent variables.")
    }

    ## LFr[,r] is an intercept term
    LFr <- FUN(dat)
    LFr1 <- LFr[, r1, drop = FALSE]
    if (r != ncol(LFr))
        stop("The number of latent variables ", r, "is not equal to the number of column(s) provided by `FUN`")
    if (!is.null(r0))
        LFr0 <- LFr[, r0, drop = FALSE]

    if ( devR ) {
        ## uses a deviance computation from R base
        obs <- dev.R(
            dat,
            LFr1 = cbind( LFr1, covariate ),
            LFr0 = cbind( LFr0, covariate )
        )
    } else {
        ## jackstraw:::devdiff
        obs <- devdiff(
            dat,
            LF_alt = cbind(LFr, covariate),
            LF_null = cbind(LFr0, matrix(1, n, 1), covariate)
        )
    }
    
    # Estimate null association
    # statistics
    null <- matrix(0, nrow = s, ncol = B)
    LFr0.js <- NULL

    if ( verbose )
        cat( paste0( "\nComputating null statistics (", B, " total iterations): " ) )
    for (i in 1:B) {
        # select the random subset of rows/loci to edit
        random.s <- sample( 1:m, size = s, replace = FALSE )
        # extract that data and permute it, returning a matrix
        s.nulls <- t( apply( dat[random.s, , drop = FALSE], 1, function(x) sample(x) ) )
        # make a copy of this whole data containing the random subset
        jackstraw.dat <- dat
        jackstraw.dat[random.s, ] <- s.nulls

        # get LFs for this random data
        LFr.js <- FUN(jackstraw.dat)
        LFr1.js <- LFr.js[, r1, drop = FALSE]
        if (!is.null(r0))
            LFr0.js <- LFr.js[, r0, drop = FALSE]

        if ( devR ) {
            ## deviance computation from R base
            null[, i] <- dev.R(
                s.nulls,
                LFr1 = cbind(LFr1.js, covariate),
                LFr0 = cbind(LFr0.js, covariate)
            )
        } else {
            ## jackstraw:::devdiff originally from lfa
            null[, i] <- devdiff(
                s.nulls,
                LF_alt = cbind(LFr.js, covariate),
                LF_null = cbind(LFr0.js, matrix(1, n, 1), covariate)
            )
        }

        if ( verbose )
            cat(paste(i, " "))
    }

    p.value <- qvalue::empPvals(as.vector(obs), as.vector(null))

    return(
        list(
            call = match.call(),
            p.value = p.value,
            obs.stat = obs,
            null.stat = null
        )
    )
}
