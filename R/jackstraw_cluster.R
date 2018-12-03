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
#' @param algorithm a clustering algorithm to use, where an output must include `cluster` and `centers`. For exact specification, see \code{kmeans}.
#' @param noise specify a parametric distribution to generate a noise term. If \code{NULL}, a non-parametric jackstraw test is performed.
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations.
#' @param center a logical specifying to center the rows. By default, \code{TRUE}.
#' @param covariate a model matrix of covariates with \code{n} observations. Must include an intercept in the first column.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param seed a seed for the random number generator.
#' @param ... optional arguments to control the clustering algorithm.
#'
#' @return \code{jackstraw_cluster} returns a list consisting of
#' \item{F.obs}{\code{m} observed F statistics between variables and cluster centers.}
#' \item{F.null}{F null statistics between null variables and cluster centers, from the jackstraw method.}
#' \item{p.F}{\code{m} p-values of membership.}
#'
#' @export jackstraw_cluster
#' @importFrom qvalue empPvals
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung (2018) Statistical significance for cluster membership. biorxiv, doi:10.1101/248633 \url{https://www.biorxiv.org/content/early/2018/01/16/248633}
jackstraw_cluster <- function(dat, 
    k, cluster = NULL, centers = NULL,
    algorithm = function(x, centers) kmeans(x, centers, ...),
    s = 1, B = 1000, center = TRUE, noise = NULL, covariate = NULL,
    verbose = FALSE, seed = NULL, 
    ...) {
    if (is.null(seed)) 
        set.seed(seed)
    m <- nrow(dat)
    n <- ncol(dat)
    
    if (is.null(cluster) | is.null(centers)) {
        stop("Supply the original cluster assignments and estimated centers.")
    }
    algorithm <- match.fun(algorithm)
    
    ## sanity check
    if (k != length(unique(cluster))) {
        stop("The input k must equal the number of clusters available.")
    }
    if (k != nrow(centers)) {
        stop("The input k must equal the number of available centers, that are rows of `center`.")
    }
    if (length(unique(cluster)) != 
        nrow(centers)) {
        stop("The number of clusters must equal the number of centers, that are rows of `center`.")
    }
    
    # compute the observed
    # statistics between observed
    # variables and cluster centers
    F.obs <- vector("numeric", m)
    for (i in 1:k) {
        F.obs[cluster == i] <- FSTAT(dat[cluster == 
            i, , drop = FALSE], 
            LV = t(centers[i, , 
                drop = FALSE]), 
            covariate = covariate)$fstat
    }
    
    if (!is.null(noise)) {
        noise <- match.fun(noise)
        message("The distribution for the noise term is specified; performing the parametric jackstraw test.")
    }

    if (verbose == TRUE) {
        cat(paste0("\nComputating null statistics (", 
            B, " total iterations): "))
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
        if (!is.null(noise)) {
            jackstraw.dat[ind, ] <- matrix(noise(n * 
                s), nrow = s, ncol = n)
        } else {
            jackstraw.dat[ind, ] <- apply(dat[ind, 
                , drop = FALSE], 1, 
                function(x) sample(x, 
                    replace = TRUE))
        }
        if(center == TRUE) {
            jackstraw.dat[ind, ] <- t(scale(t(jackstraw.dat[ind, 
                ]), center = TRUE, scale = FALSE))
        }
        
        # re-cluster the jackstraw data
        recluster <- algorithm(jackstraw.dat, 
            centers = centers, 
            ...)
        
        for (i in 1:k) {
            ind.i <- intersect(ind, 
                which(recluster$cluster == 
                  i))
            if (length(ind.i) > 
                0) {
                
                F.null[[i]] <- c(F.null[[i]], 
                  as.vector(FSTAT(dat = jackstraw.dat[ind.i, 
                    , drop = FALSE], 
                    LV = t(recluster$centers[i, 
                      , drop = FALSE]), 
                    covariate = covariate)$fstat))
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
        p.F[cluster == i] <- empPvals(F.obs[cluster == 
            i], F.null[[i]])
    }
    
    return(list(call = match.call(), 
        F.obs = F.obs, F.null = F.null, 
        p.F = p.F))
}
