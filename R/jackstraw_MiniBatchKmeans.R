#' Non-Parametric Jackstraw for Mini Batch K-means Clustering
#'
#' Test the cluster membership for K-means clustering
#'
#' K-means clustering assign \code{m} rows into \code{K} clusters. This function enable statistical
#' evaluation if the cluster membership is correctly assigned. Each of \code{m} p-values refers to
#' the statistical test of that row with regard to its assigned cluster.
#' Its resampling strategy accounts for the over-fitting characteristics due to direct computation of clusters from the observed data
#' and protects against an anti-conservative bias.
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param MiniBatchKmeans.output an output from applying \code{ClusterR::MiniBatchKmeans()} onto \code{dat}. This provides more controls over the algorithm and subsequently the initial centroids used.
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations.
#' @param center a logical specifying to center the rows. By default, \code{TRUE}.
#' @param covariate a model matrix of covariates with \code{n} observations. Must include an intercept in the first column.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param batch_size the size of the mini batches.
#' @param initializer the method of initialization. By default, \code{kmeans++}.
#' @param pool a logical specifying to pool the null statistics across all clusters. By default, \code{TRUE}.
#' @param ... optional arguments to control the Mini Batch K-means clustering algorithm (refers to \code{ClusterR::MiniBatchKmeans}).
#'
#' @return \code{jackstraw_MiniBatchKmeans} returns a list consisting of
#' \item{F.obs}{\code{m} observed F statistics between variables and cluster centers.}
#' \item{F.null}{F null statistics between null variables and cluster centers, from the jackstraw method.}
#' \item{p.F}{\code{m} p-values of membership.}
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung (2020) Statistical significance of cluster membership for unsupervised evaluation of cell identities. Bioinformatics, 36(10): 3107–3114 \url{https://academic.oup.com/bioinformatics/article/36/10/3107/5788523}
#' @examples
#' \dontrun{
#' library(ClusterR)
#' dat = t(scale(t(Jurkat293T), center=TRUE, scale=FALSE))
#' MiniBatchKmeans.output <- MiniBatchKmeans(data=dat, clusters = 2, batch_size = 300,
#' initializer = "kmeans++")
#' jackstraw.output <- jackstraw_MiniBatchKmeans(dat,
#' MiniBatchKmeans.output = MiniBatchKmeans.output)
#' }
#' 
#' @export
jackstraw_MiniBatchKmeans <- function(
                                      dat,
                                      MiniBatchKmeans.output = NULL,
                                      s = NULL,
                                      B = NULL,
                                      center = TRUE,
                                      covariate = NULL,
                                      verbose = FALSE,
                                      batch_size = floor(nrow(dat)/100),
                                      initializer = 'kmeans++',
                                      pool = TRUE,
                                      ...
                                      ) {
    # check mandatory data
    if ( missing( dat ) )
        stop( '`dat` is required!' )
    if ( missing( MiniBatchKmeans.output ) )
        stop( '`MiniBatchKmeans.output` is required!' )
    if ( !is.matrix( dat ) )
        stop( '`dat` must be a matrix!' )
    if ( !methods::is( MiniBatchKmeans.output, "k-means clustering" ) )
        stop("`MiniBatchKmeans.output` must be an object of class `k-means clustering` as a result from applying ClusterR::MiniBatchKmeans. See ?ClusterR::MiniBatchKmeans.")
    
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

    # silence this warning, which is not a real deprecation (function is not getting replaced, and we're not using the troublesome case fuzzy=TRUE anyway (default is FALSE)).  Warning Message:
    # - `predict_MBatchKMeans()` was deprecated in ClusterR 1.3.0.
    # - i Beginning from version 1.4.0, if the fuzzy parameter is TRUE the function 'predict_MBatchKMeans' will return only the probabilities, whereas currently it also returns the hard clusters
    # From NEWS: I added a deprecation warning in the ‘predict_MBatchKMeans()’ function because starting from version 1.4.0, if the ‘fuzzy’ parameter is TRUE then the function will return only the probabilities, whereas currently it also returns the hard clusters. Moreover, I added the ‘updated_output’ parameter which shows the new output format when set to TRUE.
    suppressWarnings(
        MiniBatchKmeans.output$cluster <- ClusterR::predict_MBatchKMeans( dat, MiniBatchKmeans.output$centroids )
    )
    k <- clusters <- nrow( MiniBatchKmeans.output$centroids )
    
    if ( verbose )
        cat(paste0("\nComputating null statistics (", B, " total iterations): "))

    # compute the observed
    # statistics between rows and
    # cluster centers
    F.obs <- vector("numeric", m)
    for (i in 1:k) {
        F.obs[MiniBatchKmeans.output$cluster == i] <- FSTAT(
            dat[MiniBatchKmeans.output$cluster == i, , drop = FALSE],
            LV = t(MiniBatchKmeans.output$centroids[i, , drop = FALSE]),
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
        if (center) {
            jackstraw.dat[ind, ] <- t(scale(
                t(jackstraw.dat[ind, , drop = FALSE ]),
                center = TRUE,
                scale = FALSE
            ))
        }

        # re-cluster the jackstraw data
        jackstraw.MiniBatchKmeans <- ClusterR::MiniBatchKmeans(
                                                   data=jackstraw.dat,
                                                   CENTROIDS = MiniBatchKmeans.output$centroids,
                                                   clusters = clusters,
                                                   batch_size = batch_size,
                                                   initializer = initializer,
                                                   ...
                                               )
        # suppress warning as done earlier (explanation above)
        suppressWarnings(
            jackstraw.MiniBatchKmeans$cluster <- ClusterR::predict_MBatchKMeans(jackstraw.dat, jackstraw.MiniBatchKmeans$centroids)
        )

        for (i in 1:k) {
            ind.i <- intersect(ind, which(jackstraw.MiniBatchKmeans$cluster == i))
            if (length(ind.i) > 0) {
                F.null[[i]] <- c(
                    F.null[[i]],
                    as.vector(
                        FSTAT(
                            dat = jackstraw.dat[ind.i, , drop = FALSE],
                            LV = t( jackstraw.MiniBatchKmeans$centroids[i, , drop = FALSE] ),
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
            if (length(F.null[[i]]) < (B * s/k * 0.1))
                warning( "The number of empirical null statistics for the cluster [", i, "] is [", length(F.null[[i]]), "].")
            
            p.F[MiniBatchKmeans.output$cluster == i] <- empPvals( F.obs[MiniBatchKmeans.output$cluster == i], F.null[[i]] )
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
