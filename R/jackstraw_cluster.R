#' Jackstraw for the User-Defined Clustering Algorithm
#'
#' Test the cluster membership using a user-defined clustering algorithm
#'
#' The clustering algorithms assign \code{m} rows into \code{K} clusters. This function enable statistical
#' evaluation if the cluster membership is correctly assigned. Each of \code{m} p-values refers to
#' the statistical test of that row with regard to its assigned cluster.
#' Its resampling strategy accounts for the over-fitting characteristics due to direct computation of clusters from the observed data
#' and protects against an anti-conservative bias.
#'
#' The user is expected to explore the data with a given clustering algorithm and
#' determine the number of clusters \code{k}.
#' Furthermore, provide \code{cluster} and \code{centers} as given by applying \code{algorithm} onto \code{dat}.
#' The rows of \code{centers} correspond to \code{k} clusters, as well as available levels in \code{cluster}.
#' This function allows you to specify a parametric distribution of a noise term. It is an experimental feature.
#' 
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param k a number of clusters.
#' @param cluster a vector of cluster assignments.
#' @param centers a matrix of all cluster centers.
#' @param algorithm a clustering algorithm to use, where an output must include `cluster` and `centers`. For exact specification, see \code{\link[stats]{kmeans}}.
#' @param noise specify a parametric distribution to generate a noise term. If \code{NULL}, a non-parametric jackstraw test is performed.
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations.
#' @param center a logical specifying to center the rows. By default, \code{TRUE}.
#' @param covariate a model matrix of covariates with \code{n} observations. Must include an intercept in the first column.
#' @param pool a logical specifying to pool the null statistics across all clusters. By default, \code{TRUE}.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param ... additional, optional arguments to `algorithm`.
#'
#' @return \code{jackstraw_cluster} returns a list consisting of
#' \item{F.obs}{\code{m} observed F statistics between variables and cluster centers.}
#' \item{F.null}{F null statistics between null variables and cluster centers, from the jackstraw method.}
#' \item{p.F}{\code{m} p-values of membership.}
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung (2020) Statistical significance of cluster membership for unsupervised evaluation of cell identities. Bioinformatics, 36(10): 3107–3114 \url{https://academic.oup.com/bioinformatics/article/36/10/3107/5788523}
#' 
#' @export
jackstraw_cluster <- function(
                              dat, 
                              k,
                              cluster,
                              centers,
                              algorithm = function(x, centers, ...) stats::kmeans(x, centers, ...),
                              s = 1,
                              B = 1000,
                              center = TRUE,
                              noise = NULL,
                              covariate = NULL,
                              pool = TRUE,
                              verbose = FALSE,
                              ...
                              ) {
    # check mandatory data
    if ( missing( dat ) )
        stop( '`dat` is required!' )
    if ( missing( k ) )
        stop( '`k` is required!' )
    if ( missing( cluster ) )
        stop( '`cluster` is required!' )
    if ( missing( centers ) )
        stop( '`centers` is required!' )
    if ( !is.matrix( dat ) )
        stop( '`dat` must be a matrix!' )

    m <- nrow(dat)
    n <- ncol(dat)
    
    # check additional dimensions
    if ( length(cluster) != m )
        stop( 'Length of `cluster` (', length(cluster), ') does not equal numbber of rows of data (', m , ').')
    if ( length(unique(cluster)) != k )
        stop( 'Number of clusters in `cluster` (', length(unique(cluster)), ') does not equal `k` (, k, )' )
    if ( nrow(centers) != k )
        stop( 'Number of rows in `centers` (', nrow(centers), ') does not equal `k` (', k, ')' )
    if ( ncol(centers) != n )
        stop( 'Number of columns in `centers` (', ncol(centers), ') does not equal number of columns in `dat` (', n, ')' )
    
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

    algorithm <- match.fun(algorithm)
    
    # compute the observed
    # statistics between observed
    # variables and cluster centers
    F.obs <- vector("numeric", m)
    for (i in 1:k) {
        F.obs[cluster == i] <- FSTAT(
            dat[cluster == i, , drop = FALSE], 
            LV = t( centers[i, , drop = FALSE] ), 
            covariate = covariate
        )$fstat
    }
    
    if (!is.null(noise)) {
        noise <- match.fun(noise)
        if (verbose)
            message("The distribution for the noise term is specified; performing the parametric jackstraw test.")
    }

    if (verbose)
        cat(paste0("\nComputating null statistics (", B, " total iterations): "))

    # set-up empty matrices for
    # null statistics
    F.null <- vector("list", length = k)
    for (j in 1:B) {
        if ( verbose )
            cat(paste(j, " "))

        jackstraw.dat <- dat
        # randomly choose s variables
        # to permute
        ind <- sample.int( m, s )
        if (!is.null(noise)) {
            jackstraw.dat[ind, ] <- matrix(noise(n * s), nrow = s, ncol = n)
        } else {
            jackstraw.dat[ind, ] <- apply(
                dat[ind, , drop = FALSE],
                1, 
                function(x) sample(x, replace = TRUE)
            )
        }
        if (center) {
            jackstraw.dat[ind, ] <- t(scale(
                t( jackstraw.dat[ ind, , drop = FALSE ] ),
                center = TRUE,
                scale = FALSE
            ))
        }
        
        # re-cluster the jackstraw data
        recluster <- algorithm(
            jackstraw.dat, 
            centers = centers,
            ...
        )
        
        for (i in 1:k) {
            ind.i <- intersect( ind, which(recluster$cluster == i) )
            if (length(ind.i) > 0) {
                
                F.null[[i]] <- c(
                    F.null[[i]], 
                    as.vector(
                        FSTAT(
                            dat = jackstraw.dat[ind.i, , drop = FALSE], 
                            LV = t(recluster$centers[i, , drop = FALSE]), 
                            covariate = covariate
                        )$fstat
                    )
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
            if (length(F.null[[i]]) == 0)
                stop( "There are no null statistics for the cluster ", i, ". Check if B is large enough, clusters stable, and k selected appropriately. Also consider running the algorithm with the pool option." )

            if (length(F.null[[i]]) < (B * s/k * 0.1))
                warning( "The number of empirical null statistics for the cluster [", i, "] is [", length(F.null[[i]]), "]." )
            
            p.F[cluster == i] <- empPvals( F.obs[ cluster == i ], F.null[[ i ]] )
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
