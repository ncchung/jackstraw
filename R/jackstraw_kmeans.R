#' Non-Parametric Jackstraw for K-means Clustering
#'
#' Test the cluster membership for K-means clustering
#'
#' K-means clustering assign \code{m} rows into \code{K} clusters. This function enable statistical
#' evaluation if the cluster membership is correctly assigned. Each of \code{m} p-values refers to
#' the statistical test of that row with regard to its assigned cluster.
#' Its resampling strategy accounts for the over-fitting characteristics due to direct computation of clusters from the observed data
#' and protects against an anti-conservative bias.
#'
#' The input data (\code{dat}) must be of a class `matrix`.
#'
#' @param dat a matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param kmeans.dat an output from applying \code{kmeans()} onto \code{dat}.
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations.
#' @param center a logical specifying to center the rows of the null samples. By default, \code{TRUE}.
#' @param covariate a model matrix of covariates with \code{n} observations. Must include an intercept in the first column.
#' @param match a logical specifying to match the observed clusters and jackstraw clusters using minimum Euclidean distances.
#' @param pool a logical specifying to pool the null statistics across all clusters. By default, \code{TRUE}.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param ... optional arguments to control the k-means clustering algorithm (refers to \code{kmeans}).
#'
#' @return \code{jackstraw_kmeans} returns a list consisting of
#' \item{F.obs}{\code{m} observed F statistics between variables and cluster centers.}
#' \item{F.null}{F null statistics between null variables and cluster centers, from the jackstraw method.}
#' \item{p.F}{\code{m} p-values of membership.}
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung (2020) Statistical significance of cluster membership for unsupervised evaluation of cell identities. Bioinformatics, 36(10): 3107–3114 \url{https://academic.oup.com/bioinformatics/article/36/10/3107/5788523}
#' @examples
#' \dontrun{
#' dat = t(scale(t(Jurkat293T), center=TRUE, scale=FALSE))
#' kmeans.dat <- kmeans(dat, centers=2, nstart = 10, iter.max = 100)
#' jackstraw.out <- jackstraw_kmeans(dat, kmeans.dat)
#' }
#' 
#' @export
jackstraw_kmeans <- function(
                             dat,
                             kmeans.dat,
                             s = NULL,
                             B = NULL,
                             center = FALSE,
                             covariate = NULL,
                             match = TRUE,
                             pool = TRUE,
                             verbose = FALSE,
                             ...
                             ) {
    # check mandatory data
    if ( missing( dat ) )
        stop( '`dat` is required!' )
    if ( missing( kmeans.dat ) )
        stop( '`kmeans.dat` is required!' )
    if ( !is.matrix( dat ) )
        stop( '`dat` must be a matrix!' )
    if ( !methods::is( kmeans.dat, "kmeans" ) )
        stop( "`kmeans.dat` must be an object of class `kmeans`. See ?kmeans." )

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
            message( "A number of resampling iterations (B) is not specified: B=round(m*10/s)=", B, "." )
    }
    
    k <- nrow(kmeans.dat$centers)
    obs.centers <- kmeans.dat$centers
    rownames(obs.centers) <- paste0("obs",1:k)

    if (verbose)
        cat(paste0("\nComputating null statistics (", B, " total iterations): "))

    # compute the observed
    # statistics between rows and
    # cluster centers
    F.obs <- vector("numeric", m)
    for (i in 1:k) {
        F.obs[kmeans.dat$cluster == i] <- FSTAT(
            dat[kmeans.dat$cluster == i, , drop = FALSE],
            LV = t(kmeans.dat$centers[i, , drop = FALSE]),
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
        kmeans.null <- stats::kmeans(
                                  jackstraw.dat,
                                  centers = kmeans.dat$centers,
                                  ...
                              )
        
        # with stable clusters, numeric identities of clusters are typically matched, after resampling s
        if (match) {
            jck.centers <- kmeans.null$centers
            rownames(jck.centers) <- paste0("jck",1:k)
            # min euclidean dist to match jck.centers with obs.centers
            jck.clmatch <- apply(
                jck.centers,
                1,
                function(y) which.min(apply(obs.centers, 1, function(x,y) stats::dist(rbind(x,y)),y))
            )
            if (verbose)
                message("Numeric identities of clusters are matched according to min euclidean distances.")

            for (i in 1:k) {
                ind.i <- intersect(
                    ind,
                    which(kmeans.null$cluster == which(jck.clmatch == i) )
                )

                if (verbose)
                    message( "Iteration ", j, ": jackstraw cluster ", which(jck.clmatch == i), " matched with obs cluster ", i,"; # intersections ",length(ind.i))
                if (length(ind.i) > 0) {
                    F.null[[i]] <- c(
                        F.null[[i]],
                        as.vector(
                            FSTAT(
                                dat = jackstraw.dat[ind.i, , drop = FALSE],
                                LV = t(kmeans.null$centers[which(jck.clmatch == i), , drop = FALSE]),
                                covariate = covariate
                            )$fstat
                        )
                    )
                }
            }
        } else {
            # no matching required
            for (i in 1:k) {
                ind.i <- intersect(ind, which(kmeans.null$cluster == i))
                if (length(ind.i) > 0) {
                    F.null[[i]] <- c(
                        F.null[[i]],
                        as.vector(
                            FSTAT(
                                dat = jackstraw.dat[ind.i, , drop = FALSE],
                                LV = t(kmeans.null$centers[i, , drop = FALSE]),
                                covariate = covariate
                            )$fstat
                        )
                    )
                }
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

            p.F[kmeans.dat$cluster == i] <- empPvals( F.obs[kmeans.dat$cluster == i], F.null[[i]] )
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
