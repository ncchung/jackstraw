#' Non-Parametric Jackstraw for ALStructure
#'
#' Test association between the observed variables and population structure estimated by ALStructure.
#'
#' This function uses ALStructure from Cabreros and Storey (2019). A deviation \code{dev} in logistic regression
#' (the full model with \code{r} LFs vs. the intercept-only model) is used to assess association.
#'
#' @param dat a genotype matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param r a number of significant LFs.
#' @param FUN a function to ALStructure
#' @param r1 a numeric vector of LFs of interest (implying you are not interested in all \code{r} LFs).
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations. There will be a total of \code{s*B} null statistics.
#' @param covariate a data matrix of covariates with corresponding \code{n} observations (do not include an intercept term).
#' @param verbose a logical specifying to print the computational progress.
#'
#' @return \code{jackstraw_alstructure} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their LFs}
#' \item{obs.stat}{\code{m} observed deviances}
#' \item{null.stat}{\code{s*B} null deviances}
#'
#' @examples
#' \dontrun{
#' # load genotype data to analyze (not shown) into this variable
#' X
#' # choose the number of ancestries
#' r <- 3
#' 
#' # load alstructure package (install from https://github.com/StoreyLab/alstructure)
#' library(alstructure)
#' # define the function this way, a function of the genotype matrix only
#' FUN <- function(x) t( alstructure(x, d_hat = r)$Q_hat )
#'
#' # calculate p-values (and other statistics) for each SNP
#' out <- jackstraw_alstructure( X, r, FUN )
#' }
#'
#' @references Chung and Storey (2015) Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics, 31(4): 545-554 \url{https://academic.oup.com/bioinformatics/article/31/4/545/2748186}
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#'
#' @seealso  \link{jackstraw_pca} \link{jackstraw}
#' 
#' @export
jackstraw_alstructure <- function(
                                  dat,
                                  r,
                                  FUN,
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
    if ( missing( FUN ) )
        stop( '`FUN` is required!' )
    if ( !is.matrix( dat ) )
        stop( '`dat` must be a matrix!' )
    if ( !is.function( FUN ) )
        stop( '`FUN` must be a function!' )
    
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

    if (is.null(s)) {
        s <- round(m/10)
        if (verbose)
            message("A number of null variables (s) to be permuted is not specified: s=round(0.10*m)=", s, ".")
    }
    if (is.null(B)) {
        B <- round(m * 10/s)
        if (verbose)
            message("A number of resampling iterations (B) is not specified: B=round(m*10/s)=", B, ".")
    }

    if (is.null(r1))
        r1 <- 1:r
    if (all(seq(r) %in% r1)) {
        # no adjustment LVs
        r0 <- NULL
        LFr0 <- NULL
    } else {
        r0 <- seq(r)[-r1]
    }

    if (!is.null(FUN)) {
        FUN <- match.fun(FUN)
    } else {
        stop("Please provide a function to estimate latent variables.")
    }

    ## ALS does not have an intercept term like LFA
    LFr <- FUN(dat)
    LFr1 <- LFr[, r1, drop = FALSE]
    if (r != ncol(LFr))
        stop( "The number of latent variables ", r, "is not equal to the number of column(s) provided by `FUN`" )
    
    if (!is.null(r0))
        LFr0 <- LFr[, r0, drop = FALSE]

    # NOTE: there are some issues here, see `jacsktraw_lfa` for notes
    obs <- gcatest::delta_deviance_lf(
                        X = dat,
                        LF0 = cbind(LFr0, matrix(1, n, 1), covariate),
                        LF1 = cbind(LFr, covariate)
                    )
    
    # Estimate null association
    # statistics
    null <- matrix(0, nrow = s, ncol = B)
    LFr0.js <- NULL

    if (verbose)
        cat(paste0("\nComputating null statistics (", B, " total iterations): "))
    for (i in 1:B) {
        random.s <- sample.int( m, s )
        s.nulls <- dat[ random.s, , drop = FALSE ]
        s.nulls <- t( apply( s.nulls, 1, sample ) )
        jackstraw.dat <- dat
        jackstraw.dat[random.s, ] <- s.nulls

        LFr.js <- FUN(jackstraw.dat)
        LFr1.js <- LFr.js[, r1, drop = FALSE]
        if (!is.null(r0))
            LFr0.js <- LFr.js[, r0, drop = FALSE]

        null[, i] <- gcatest::delta_deviance_lf(
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
