#' Non-Parametric Jackstraw for Partitioning Around Medoids (PAM)
#'
#' Test the cluster membership for Partitioning Around Medoids (PAM)
#'
#' PAM assigns \code{m} rows into \code{K} clusters. This function enable statistical
#' evaluation if the cluster membership is correctly assigned. Each of \code{m} p-values refers to
#' the statistical test of that row with regard to its assigned cluster.
#' Its resampling strategy accounts for the over-fitting characteristics due to direct computation of clusters from the observed data
#' and protects against an anti-conservative bias.
#'
#' For a large dataset, PAM could be too slow. Consider using \code{cluster::clara} and \code{jackstraw::jackstraw_clara}.
#'
#' The input data (\code{dat}) must be of a class `matrix`.
#'
#' @param dat a matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param pam.dat an output from applying \code{cluster::pam()} on \code{dat}.
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations.
#' @param center a logical specifying to center the rows. By default, \code{TRUE}.
#' @param covariate a model matrix of covariates with \code{n} observations. Must include an intercept in the first column.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param pool a logical specifying to pool the null statistics across all clusters. By default, \code{TRUE}.
#' @param ... optional arguments to control the k-means clustering algorithm (refers to \code{kmeans}).
#'
#' @return \code{jackstraw_pam} returns a list consisting of
#' \item{F.obs}{\code{m} observed F statistics between variables and cluster medoids.}
#' \item{F.null}{F null statistics between null variables and cluster medoids, from the jackstraw method.}
#' \item{p.F}{\code{m} p-values of membership.}
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung (2020) Statistical significance of cluster membership for unsupervised evaluation of cell identities. Bioinformatics, 36(10): 3107–3114 \url{https://academic.oup.com/bioinformatics/article/36/10/3107/5788523}
#' @examples
#' \dontrun{
#' library(cluster)
#' dat = t(scale(t(Jurkat293T), center=TRUE, scale=FALSE))
#' pam.dat <- pam(dat, k=2)
#' jackstraw.out <- jackstraw_pam(dat, pam.dat = pam.dat)
#' }
#' 
#' @export
jackstraw_pam <- function(
                          dat,
                          pam.dat,
                          s = NULL,
                          B = NULL,
                          center = TRUE,
                          covariate = NULL,
                          verbose = FALSE,
                          pool = TRUE,
                          ...
                          ) {
    # check mandatory data
    if ( missing( dat ) )
        stop( '`dat` is required!' )
    if ( missing( pam.dat ) )
        stop( '`pam.dat` is required!' )
    if ( !is.matrix( dat ) )
        stop( '`dat` must be a matrix!' )
    if ( !methods::is( pam.dat, "pam" ) )
        stop( "`pam.dat` must be an object of class `pam`. See ?pam.object." )
    
    m <- nrow(dat)
    n <- ncol(dat)

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
            message( "A number of null variables (s) to be permuted is not specified: s=round(0.10*m)=", s, "." )
    }
    if (is.null(B)) {
        B <- round(m * 10/s)
        if (verbose)
            message("A number of resampling iterations (B) is not specified: B=round(m*10/s)=", B, "." )
    }

    k <- nrow(pam.dat$medoids)

    if (verbose)
        cat(paste0("\nComputating null statistics (", B, " total iterations): "))
    
    # compute the observed
    # statistics between rows and
    # cluster medoids
    F.obs <- vector("numeric", m)
    for (i in 1:k) {
        F.obs[pam.dat$clustering == i] <- FSTAT(
            dat[pam.dat$clustering == i, , drop = FALSE],
            LV = t(pam.dat$medoids[i, , drop = FALSE]),
            covariate = covariate
        )$fstat
    }

    # set-up empty matrices for
    # null statistics
    F.null <- vector("list", length = k)
    for (j in 1:B) {
        if (verbose)
            cat(paste(j, " "))
        
        jackstraw.dat <- dat
        # randomly choose s variables
        # to permute
        ind <- sample.int( m, s )
        jackstraw.dat[ind, ] <- apply(
            dat[ind, , drop = FALSE],
            1,
            function(x) sample(x, replace = TRUE)
        )
        if (center)
            jackstraw.dat[ind, ] <- t(scale(
                t( jackstraw.dat[ ind, , drop = FALSE ] ),
                center = TRUE,
                scale = FALSE
            ))
        
        # re-cluster the jackstraw data
        pam.null <- cluster::pam(
                                 jackstraw.dat,
                                 k = k,
                                 medoids = pam.dat$id.med,
                                 ...
                             )
        
        for (i in 1:k) {
            ind.i <- intersect(
                ind,
                which( pam.null$clustering == i )
            )
            if ( length(ind.i) > 0 ) {
                F.null[[i]] <- c(
                    F.null[[i]],
                    FSTAT(
                        dat = jackstraw.dat[ind.i, , drop = FALSE],
                        LV = t( pam.null$medoids[i, , drop = FALSE] ), # colinearity here (in toy cases)???
                        covariate = covariate
                    )$fstat
                )
            }
        }
    }

    # compute p-values
    p.F <- vector("numeric", m)
    if (pool) {
        p.F <- empPvals( F.obs, unlist( F.null ) )
    } else {
        for (i in 1:k) {
            # warn about a relatively low
            # number of null statistics
            if (length(F.null[[i]]) < (B * s/k * 0.1))
                warning( "The number of empirical null statistics for the cluster [", i, "] is [", length(F.null[[i]]), "].")
            p.F[pam.dat$clustering == i] <- empPvals( F.obs[pam.dat$clustering == i], F.null[[i]] )
        }
    }
    
    return(
        list(
            call = match.call(),
            F.obs = F.obs,
            F.null = F.null,
            p.F = p.F
        )
    )
}
