#' Non-Parametric Jackstraw for Principal Component Analysis (PCA) using the augmented implicitly restarted Lanczos bidiagonalization algorithm (IRLBA)
#'
#' Test association between the observed variables and their latent variables captured by principal components (PCs). PCs are computed using the augmented implicitly restarted Lanczos bidiagonalization algorithm (IRLBA; see \code{\link[irlba]{irlba}}).
#'
#' This function computes \code{m} p-values of linear association between \code{m} variables and their PCs.
#' Its resampling strategy accounts for the over-fitting characteristics due to direct computation of PCs from the observed data
#' and protects against an anti-conservative bias.
#'
#' Provide the data matrix, with \code{m} variables as rows and \code{n} observations as columns.
#' Given that there are \code{r} significant PCs, this function tests for linear association between
#' \code{m} varibles and their \code{r} PCs.
#'
#' You could specify a subset of significant PCs that you are interested in (\code{PC}). If \code{PC} is given,
#' then this function computes statistical significance of association between \code{m} variables and \code{PC},
#' while adjusting for other PCs (i.e., significant PCs that are not your interest).
#' For example, if you want to identify variables associated with 1st and 2nd PCs,
#' when your data contains three significant PCs, set \code{r=3} and \code{PC=c(1,2)}.
#'
#' Please take a careful look at your data and use appropriate graphical and statistical criteria
#' to determine a number of significant PCs, \code{r}. The number of significant PCs depends on the data structure and the context.
#' In a case when you fail to specify \code{r}, it will be estimated from a permutation test (Buja and Eyuboglu, 1992)
#' using a function \link{permutationPA}.
#'
#' If \code{s} is not supplied, \code{s} is set to about 10\% of \code{m} variables.
#' If \code{B} is not supplied, \code{B} is set to \code{m*10/s}.
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param r1 a numeric vector of principal components of interest. Choose a subset of \code{r} significant PCs to be used.
#' @param r a number (a positive integer) of significant principal components. See \link{permutationPA} and other methods.
#' @param s a number (a positive integer) of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number (a positive integer) of resampling iterations. There will be a total of \code{s*B} null statistics.
#' @param covariate a data matrix of covariates with corresponding \code{n} observations (do not include an intercept term).
#' @param verbose a logical specifying to print the computational progress.
#' @param seed a numeric seed for the random number generator.
#' @param ... additional arguments to \code{\link[irlba]{irlba}}.
#'
#' @return \code{jackstraw_irlba} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their principal components}
#' \item{obs.stat}{\code{m} observed F-test statistics}
#' \item{null.stat}{\code{s*B} null F-test statistics}
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2015) Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics, 31(4): 545-554 \url{http://bioinformatics.oxfordjournals.org/content/31/4/545}
#'
#' @seealso \link{jackstraw} \link{jackstraw_subspace} \link{permutationPA}
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
#' ## apply the jackstraw
#' out = jackstraw_irlba(dat, r=1)
#'
#' ## Use optional arguments
#' ## For example, set s and B for a balance between speed of the algorithm and accuracy of p-values
#' \dontrun{
#' ## out = jackstraw_irlba(dat, r=1, s=10, B=1000, seed=5678)
#' }
#' 
#' @export
jackstraw_irlba <- jackstraw_tpca <- function(dat,
    r1 = NULL, r = NULL, s = NULL,
    B = NULL, covariate = NULL,
    verbose = TRUE, seed = NULL, ...) {
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
    if (!is.null(seed))
        set.seed(seed)
    if (is.null(r)) {
        warning("The number of significant PCs (r) is missing; this is strongly advised to determine r using appropriate statistical and graphical criteria.")
        r <- permutationPA(dat = dat,
            threshold = 0.05, verbose = verbose)$r
        message(paste0("Permutation Parallel Analysis, with a threshold of 0.05, estimated r = ",
            r, "."))
    }
    if (!(r > 0 && r < n)) {
        stop("r is not in valid range between 1 and n-1.")
    }
    if (is.null(r1))
        r1 <- 1:r
    if (all(seq(r) %in% r1)) {
        # no adjustment LVs
        r0 <- NULL
        ALV <- NULL
    } else {
        r0 <- seq(r)[-r1]
    }

    # Calculate observed
    # association statistics
    svd.dat <- irlba::irlba(dat, nv=r, ...)
    LV <- svd.dat$v[, r1, drop = FALSE]
    if (!is.null(r0))
        ALV <- svd.dat$v[, r0,
            drop = FALSE]

    obs <- FSTAT(dat = dat, LV = LV,
        ALV = ALV, covariate = covariate)$fstat

    # Estimate null association
    # statistics
    null <- matrix(0, nrow = s,
        ncol = B)
    ALV.js <- NULL

    if (verbose == TRUE)
        cat(paste0("\nComputating null statistics (",
            B, " total iterations): "))
    for (i in 1:B) {
        random.s <- sample(1:m,
            size = s, replace = FALSE)
        s.nulls <- t(apply(dat[random.s,
            , drop = FALSE], 1,
            function(x) sample(x)))
        jackstraw.dat <- dat
        jackstraw.dat[random.s, ] <- s.nulls

        svd.jackstraw.dat <- irlba::irlba(jackstraw.dat, nv=r, ...)
        LV.js <- svd.jackstraw.dat$v[,
            r1, drop = FALSE]
        if (!is.null(r0))
            ALV.js <- svd.jackstraw.dat$v[,
                r0, drop = FALSE]
        null[, i] <- FSTAT(dat = s.nulls,
            LV = LV.js, ALV = ALV.js,
            covariate = covariate)$fstat

        if (verbose == TRUE)
            cat(paste(i, " "))
    }

    p.value <- qvalue::empPvals(as.vector(obs), as.vector(null))

    return(list(call = match.call(),
        p.value = p.value, obs.stat = obs,
        null.stat = null))
}
