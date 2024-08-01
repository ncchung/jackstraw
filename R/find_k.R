#' Find a number of clusters or principal components
#'
#' There are a wide range of algorithms and visual techniques to identify
#' a number of clusters or principal components embedded in the observed data.
#'
#' It is critical to explore the eigenvalues, cluster stability, and visualization.
#' See R packages \code{bootcluster}, \code{EMCluster}, and \code{nFactors}.
#'
#' Please see the R package \code{SC3}, which provides \code{estkTW()} function to
#' find the number of significant eigenvalues according to the Tracy-Widom test.
#'
#' \code{ADPclust} package includes \code{adpclust()} function that runs the algorithm
#' on a range of K values. It helps you to identify the most suitable number of clusters.
#'
#' This package also provides an alternative methods in \code{permutationPA}.
#' Through a resampling-based Parallel Analysis, it finds a number of significant components.
#'
#' @export
find_k <- function() {
  print("See ? find_k for helpful functions.")
  return(NULL)
}

#' Permutation Parallel Analysis
#'
#' Estimate a number of significant principal components from a permutation test.
#'
#' Adopted from \code{sva::num.sv}, and based on Buja and Eyuboglu (1992)
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param B a number (a positive integer) of resampling iterations.
#' @param threshold a numeric value between 0 and 1 to threshold p-values.
#' @param verbose a logical indicator as to whether to print the progress.
#'
#' @return \code{permutationPA} returns
#' \item{r}{an estimated number of significant principal components based on thresholding p-values at \code{threshold}}
#' \item{p}{a list of p-values for significance of principal components}
#'
#' @references Buja A and Eyuboglu N. (1992) Remarks on parallel analysis. Multivariate Behavioral Research, 27(4), 509-540
#' @export
permutationPA <- function(
                          dat,
                          B = 100,
                          threshold = 0.05,
                          verbose = TRUE
                          ) {
    # check mandatory data
    if ( missing( dat ) )
        stop( '`dat` is required!' )
    if ( !is.matrix( dat ) )
        stop( '`dat` must be a matrix!' )
    
    n <- ncol(dat)
    m <- nrow(dat)

    uu <- corpcor::fast.svd(dat, tol = 0)
    ndf <- n - 1
    # looks like a sort of variance explained
    dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
    
    # draw random samples, to get null statistics from?
    dstat0 <- matrix( 0, nrow = B, ncol = ndf )
    if ( verbose )
        message("Estimating a number of significant principal component: ")
    for (i in 1:B) {
        if ( verbose )
            cat(paste(i, " "))
        dat0 <- t( apply( dat, 1, sample ) )
        uu0 <- corpcor::fast.svd(dat0, tol = 0)
        dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
    }
    # calculate p-values
    p <- rep(1, n)
    for ( i in 1:ndf ) {
        p[i] <- mean( dstat0[, i] >= dstat[i] )
    }
    # monotonize p-values
    for ( i in 2:ndf ) {
        p[ i ] <- max( p[ i - 1 ], p[ i ] )
    }
    r <- sum( p <= threshold )
    return( list( r = r, p = p ) )
}

