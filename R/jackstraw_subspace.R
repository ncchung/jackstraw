#' Jackstraw for the User-Defined Dimension Reduction Methods
#'
#' Test association between the observed variables and their latent variables, captured by a user-defined dimension reduction method.
#'
#' This function computes \code{m} p-values of linear association between \code{m} variables and their latent variables,
#' captured by a user-defined dimension reduction method.
#' Its resampling strategy accounts for the over-fitting characteristics due to direct computation of PCs from the observed data
#' and protects against an anti-conservative bias.
#' 
#' This function allows you to specify a parametric distribution of a noise term. It is an experimental feature. Then, a small number \code{s} of observed variables
#' are replaced by synthetic null variables generated from a specified distribution.
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param r a number of significant latent variables.
#' @param FUN Provide a specific function to estimate LVs. Must output \code{r} estimated LVs in a \code{n*r} matrix.
#' @param noise specify a parametric distribution to generate a noise term. If \code{NULL}, a non-parametric jackstraw test is performed.
#' @param r1 a numeric vector of latent variables of interest.
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations.
#' @param covariate a model matrix of covariates with \code{n} observations. Must include an intercept in the first column.
#' @param verbose a logical specifying to print the computational progress.
#'
#' @return \code{jackstraw_subspace} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their principal components}
#' \item{obs.stat}{\code{m} observed statistics}
#' \item{null.stat}{\code{s*B} null statistics}
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2015) Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics, 31(4): 545-554 \url{https://academic.oup.com/bioinformatics/article/31/4/545/2748186}
#' @references Chung (2020) Statistical significance of cluster membership for unsupervised evaluation of cell identities. Bioinformatics, 36(10): 3107–3114 \url{https://academic.oup.com/bioinformatics/article/36/10/3107/5788523}
#'
#' @seealso \link{jackstraw_pca} \link{jackstraw}
#'
#' @examples
#' ## simulate data from a latent variable model: Y = BL + E
#' B = c(rep(1,50),rep(-1,50), rep(0,900))
#' L = rnorm(20)
#' E = matrix(rnorm(1000*20), nrow=1000)
#' dat = B %*% t(L) + E
#' dat = t(scale(t(dat), center=TRUE, scale=TRUE))
#'
#' ## apply the jackstraw with the svd as a function
#' out = jackstraw_subspace(dat, FUN = function(x) svd(x)$v[,1,drop=FALSE], r=1, s=100, B=50)
#' 
#' @export
jackstraw_subspace <- function(
                               dat, 
                               r, 
                               FUN,
                               r1 = NULL,
                               s = NULL,
                               B = NULL, 
                               covariate = NULL,
                               noise = NULL,
                               verbose = TRUE
                               ) {
    # check mandatory data
    if ( missing( dat ) )
        stop( '`dat` is required!' )
    if ( missing( r ) )
        stop( '`r` is required!' )
    if ( missing( FUN ) )
        stop( "`FUN`, a function to estimate latent variables, is required!" )
    if ( !is.matrix( dat ) )
        stop( '`dat` must be a matrix!' )
    if ( !is.function( FUN ) )
        stop( "`FUN` must be a function!" )
    
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
            message( "Number of null variables (s) to be permuted is not specified: s=round(0.10*m)=", s, "." )
    }
    if (is.null(B)) {
        B <- round(m * 10/s)
        if (verbose)
            message( "Number of resampling iterations (B) is not specified: B=round(m*10/s)=", B, "." )
    }

    FUN <- match.fun(FUN)
    LV <- FUN(dat)
    if (r != ncol(LV)) 
        stop( "The number of latent variables ", r, " is not equal to the number of column(s) provided by FUN: ", ncol(LV) )
    
    if (is.null(r1)) 
        r1 <- 1:r
    
    if (all(seq(r) %in% r1)) {
        # no adjustment LVs
        r0 <- NULL
        ALV <- NULL
        ALV.js <- NULL
    } else {
        # r0 adjustment LVs
        r0 <- seq(r)[-r1]
        ALV <- LV[, r0, drop = FALSE]
        LV <- LV[, r1, drop = FALSE]  ## note that LV firstly contained the r latent variables; then reduced to the r1 latent variables of interest.
    }
    
    obs <- FSTAT(
        dat = dat,
        LV = LV, 
        ALV = ALV,
        covariate = covariate
    )$fstat
    
    if (!is.null(noise)) {
        noise <- match.fun(noise)
        message("The distribution for the noise term is specified; performing the parametric jackstraw test.")
    }
    
    if ( verbose )
        cat(paste0("\nComputating null statistics (", B, " total iterations): "))

    # Estimate null association
    # statistics
    null <- matrix(0, nrow = s, ncol = B)
    for (i in 1:B) {
        if ( verbose )
            cat(paste(i, " "))

        random.s <- sample.int( m, s )
        if (!is.null(noise)) {
            s.nulls <- matrix( noise(n * s), nrow = s, ncol = n )
        } else {
            s.nulls <- dat[ random.s, , drop = FALSE ]
            s.nulls <- t( apply( s.nulls, 1, sample ) )
        }
        jackstraw.dat <- dat
        jackstraw.dat[random.s, ] <- s.nulls
        
        LV.js <- FUN(jackstraw.dat)
        if (!is.null(r0)) {
            ALV.js <- LV.js[, r0, drop = FALSE]
            LV.js <- LV.js[, r1, drop = FALSE]
        }
        null[, i] <- FSTAT(
            dat = s.nulls, 
            LV = LV.js,
            ALV = ALV.js, 
            covariate = covariate
        )$fstat
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

