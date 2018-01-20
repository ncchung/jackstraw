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
#' @param FUN optionally, provide a specfic function to estimate LVs. Must output \code{r} estimated LVs in a \code{n*r} matrix.
#' @param noise specify a parametric distribution to generate a noise term. If \code{NULL}, a non-parametric jackstraw test is performed.
#' @param r a number of significant latent variables.
#' @param r1 a numeric vector of latent variables of interest.
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations.
#' @param covariate a model matrix of covariates with \code{n} observations. Must include an intercept in the first column.
#' @param verbose a logical specifying to print the computational progress.
#' @param seed a seed for the random number generator.
#'
#' @return \code{jackstraw_subspace} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their principal components}
#' \item{obs.stat}{\code{m} observed statistics}
#' \item{null.stat}{\code{s*B} null statistics}
#'
#' @importFrom corpcor fast.svd
#' @export jackstraw_subspace
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2015) Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics, 31(4): 545-554 \url{http://bioinformatics.oxfordjournals.org/content/31/4/545}
#' @references Chung (2018) Statistical significance for cluster membership. biorxiv, doi:10.1101/248633 \url{https://www.biorxiv.org/content/early/2018/01/16/248633}
#'
#' @seealso \link{jackstraw_pca} \link{jackstraw}
#'
#' @examples
#' set.seed(1234)
#' ## simulate data from a latent variable model: Y = BL + E
#' B = c(rep(1,50),rep(-1,50), rep(0,900))
#' L = rnorm(20)
#' E = matrix(rnorm(1000*20), nrow=1000)
#' dat = B %*% t(L) + E
#' dat = t(scale(t(dat), center=TRUE, scale=TRUE))
#'
#' ## apply the jackstraw with the svd as a function
#' out = jackstraw_subspace(dat, FUN = function(x) svd(x)$v[,1,drop=FALSE], r=1, s=100, B=50)
jackstraw_subspace <- function(dat, 
    FUN, r = NULL, 
    r1 = NULL, s = NULL, B = NULL, 
    covariate = NULL, noise = NULL, verbose = TRUE, 
    seed = NULL) {
    if (!is.null(seed)) 
        set.seed(seed)
    if (is.null(r)) 
        stop("Must provide a number of latent variables, r.")
    
    m <- dim(dat)[1]
    n <- dim(dat)[2]
    if (is.null(s)) {
        s <- round(m/10)
        message(paste0("A number of null variables (s) to be permuted is not specified: s=round(0.10*m)=", 
            s, "."))
    }
    if (is.null(B)) {
        B <- round(m * 10/s)
        message(paste0("A number of resampling iterations (B) is not specified: B=round(m*10/s)=", 
            B, "."))
    }
    if (is.null(FUN)) {
        stop("Please provide a function to estimate latent variables.")
    }
    
    FUN <- match.fun(FUN)
    LV <- FUN(dat)
    if (r != ncol(LV)) 
        stop(paste0("The number of latent variables ", 
            r, "is not equal to the number of column(s) provided by ", 
            FUN))
    
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
    obs <- FSTAT(dat = dat, LV = LV, 
        ALV = ALV, covariate = covariate)$fstat
    
    if (!is.null(noise)) {
        noise <- match.fun(noise)
        message("The distribution for the noise term is specified; performing the parametric jackstraw test.")
    }
    
    if (verbose == TRUE) {
        cat(paste0("\nComputating null statistics (", 
            B, " total iterations): "))
    }

    # Estimate null association
    # statistics
    null <- matrix(0, nrow = s, ncol = B)
    for (i in 1:B) {
        if (verbose == TRUE) {
            cat(paste(i, " "))
        }

        random.s <- sample(1:m, 
            size = s, replace = FALSE)
        if (!is.null(noise)) {
            s.nulls <- matrix(noise(n * 
                s), nrow = s, ncol = n)
        } else {
            s.nulls <- t(apply(dat[random.s, 
                , drop = FALSE], 
                1, function(x) sample(x)))
        }
        jackstraw.dat <- dat
        jackstraw.dat[random.s, ] <- s.nulls
        
        LV.js <- FUN(jackstraw.dat)
        if (!is.null(r0)) {
            ALV.js <- LV.js[, r0, 
                drop = FALSE]
            LV.js <- LV.js[, r1, 
                drop = FALSE]
        }
        null[, i] <- FSTAT(dat = s.nulls, 
            LV = LV.js, ALV = ALV.js, 
            covariate = covariate)$fstat
    }
    
    p.value <- empPvals(as.vector(obs), as.vector(null))
    
    return(list(call = match.call(), 
        p.value = p.value, obs.stat = obs, 
        null.stat = null))
}

