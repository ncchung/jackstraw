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
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param kmeans.dat an output from applying \code{kmeans()} onto \code{dat}.
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations.
#' @param covariate a model matrix of covariates with \code{n} observations. Must include an intercept in the first column.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param seed a seed for the random number generator.
#' @param ... optional arguments to control the k-means clustering algorithm (refers to \code{kmeans}).
#'
#' @return \code{jackstraw_kmeans} returns a list consisting of
#' \item{F.obs}{\code{m} observed F statistics between variables and cluster centers.}
#' \item{F.null}{F null statistics between null variables and cluster centers, from the jackstraw method.}
#' \item{p.F}{\code{m} p-values of membership.}
#'
#' @export jackstraw_kmeans
#' @importFrom qvalue empPvals
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung (2018) Statistical significance for cluster membership. biorxiv, doi:10.1101/248633 \url{https://www.biorxiv.org/content/early/2018/01/16/248633}
jackstraw_kmeans <- function(dat, 
    kmeans.dat, s = 1, B = 1000, 
    covariate = NULL, verbose = FALSE, 
    seed = NULL, ...) {
    if (is.null(seed)) 
        set.seed(seed)
    m <- nrow(dat)
    n <- ncol(dat)
    
    ## sanity check
    if (class(kmeans.dat) != "kmeans") {
        stop("The class of `kmeans.dat` must `kmeans'.")
    }
    k <- nrow(kmeans.dat$centers)
    
    if (verbose == TRUE) {
        cat(paste0("\nComputating null statistics (", 
            B, " total iterations): "))
    }

    # compute the observed
    # statistics between rows and
    # cluster centers
    F.obs <- vector("numeric", 
        m)
    for (i in 1:k) {
        F.obs[kmeans.dat$cluster == 
            i] <- FSTAT(dat[kmeans.dat$cluster == 
            i, , drop = FALSE], 
            LV = t(kmeans.dat$centers[i, 
                , drop = FALSE]), 
            covariate = covariate)$fstat
    }
    
    # set-up empty matrices for
    # null statistics
    F.null <- vector("list", length = k)
    for (j in 1:B) {
        if (verbose == TRUE) {
            cat(paste(j, " "))
        }

        jackstraw.dat <- dat
        # randomly choose s variables
        # to permute
        ind <- sample(seq(m), s)
        jackstraw.dat[ind, ] <- apply(dat[ind, 
            , drop = FALSE], 1, 
            function(x) sample(x, 
                replace = TRUE))
        jackstraw.dat[ind, ] <- t(scale(t(jackstraw.dat[ind, 
            ]), center = TRUE, 
            scale = FALSE))
        
        # re-cluster the jackstraw data
        kmeans.null <- kmeans(jackstraw.dat, 
            centers = kmeans.dat$centers, 
            ...)
        
        for (i in 1:k) {
            ind.i <- intersect(ind, 
                which(kmeans.null$cluster == 
                  i))
            if (length(ind.i) > 
                0) {

                F.null[[i]] <- c(F.null[[i]], 
                  as.vector(FSTAT(dat = jackstraw.dat[ind.i, 
                    , drop = FALSE], 
                    LV = t(kmeans.null$centers[i, 
                      , drop = FALSE]), 
                    , covariate = covariate)$fstat))
            }
        }
    }
    
    # compute p-values
    p.F <- vector("numeric", m)
    for (i in 1:k) {
        # warn about a relatively low
        # number of null statistics
        if (length(F.null[[i]]) < 
            (B * s/k * 0.1)) {
            warning(paste0("The number of empirical null statistics for the cluster [", 
                i, "] is [", length(F.null[[i]]), 
                "]."))
        }
        p.F[kmeans.dat$cluster == 
            i] <- empPvals(F.obs[kmeans.dat$cluster == 
            i], F.null[[i]])
    }
    
    return(list(call = match.call(), 
        F.obs = F.obs, F.null = F.null, 
        p.F = p.F))
}
