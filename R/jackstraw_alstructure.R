#' Non-Parametric Jackstraw for ALStructure
#'
#' Test association between the observed variables and population structure estimated by ALStructure.
#'
#' This function uses ALStructure from Cabreros and Storey (2019). A deviation \code{dev} in logistic regression
#' (the full model with \code{r} LFs vs. the intercept-only model) is used to assess association.
#'
#' @param dat a genotype matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param FUN a function to ALStructure
#' @param devR use a R function to compute deviance. By default, FALSE (uses C++).
#' @param r a number of significant LFs.
#' @param r1 a numeric vector of LFs of interest (implying you are not interested in all \code{r} LFs).
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations. There will be a total of \code{s*B} null statistics.
#' @param covariate a data matrix of covariates with corresponding \code{n} observations (do not include an intercept term).
#' @param verbose a logical specifying to print the computational progress.
#' @param seed a seed for the random number generator.
#'
#' @return \code{jackstraw_alstructure} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their LFs}
#' \item{obs.stat}{\code{m} observed devs}
#' \item{null.stat}{\code{s*B} null devs}
#'
#' @importFrom corpcor fast.svd
#' @importFrom qvalue empPvals
#' @importFrom alstructure alstructure
#' @export jackstraw_alstructure
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#'
#' @seealso  \link{jackstraw_pca} \link{jackstraw}
jackstraw_alstructure <- function(dat,
    FUN = function(x) t(alstructure(x, d_hat = r, svd_method = "truncated_svd", tol = 0.001, max_iters = 1000)$Q_hat[, , drop = FALSE]),
    devR = FALSE, r = NULL, r1 = NULL, s = NULL,
    B = NULL, covariate = NULL,
    verbose = TRUE, seed = NULL) {
    if (!is.null(seed))
        set.seed(seed)

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

    if (!(r > 0 && r < n)) {
        stop("r is not in valid range between 1 and n-1.")
    }
    if (is.null(r1))
        r1 <- 1:r
    if (all(seq(r) %in% r1)) {
        # no adjustment LVs
        r0 <- NULL
        LFr0 <- NULL
    } else {
        r0 <- seq(r)[-r1]
    }

    if (!is.null(FUN)) {
        FUN <- match.fun(FUN)
    } else {
        stop("Please provide a function to estimate latent variables.")
    }

    ## ALS does not have an intercept term like LFA
    LFr <- FUN(dat)
    LFr1 <- LFr[, r1, drop = FALSE]
    if (r != ncol(LFr))
        stop(paste0("The number of latent variables ",
            r, "is not equal to the number of column(s) provided by ",
            FUN))
    if (!is.null(r0))
        LFr0 <- LFr[, r0, drop = FALSE]

    if (devR == FALSE) {
        ## see jackstraw:::devdiff
        obs <- devdiff(dat,
                       LF_alt = cbind(LFr, covariate),
                       LF_null = cbind(LFr0, matrix(1, n, 1), covariate))
    } else {
        ## compute a deviance from R base
        obs <- dev.R(dat,
                     LFr1 = cbind(LFr1, covariate),
                     LFr0 = cbind(LFr0, covariate))
    }

    # Estimate null association
    # statistics
    null <- matrix(0, nrow = s,
        ncol = B)
    LFr0.js <- NULL

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

        LFr.js <- FUN(jackstraw.dat)
        LFr1.js <- LFr.js[, r1,
            drop = FALSE]
        if (!is.null(r0))
            LFr0.js <- LFr.js[,
                r0, drop = FALSE]

        if (devR == FALSE) {
            ## see jackstraw:::devdiff
            null[, i] <- devdiff(s.nulls,
                LF_alt = cbind(LFr.js, covariate),
                LF_null = cbind(LFr0.js, matrix(1, n, 1), covariate))
        } else {
            ## uses a deviance computation
            ## function from R base
            null[, i] <- dev.R(s.nulls,
                LFr1 = cbind(LFr1.js, covariate),
                LFr0 = cbind(LFr0.js, covariate))
        }

        if (verbose == TRUE)
            cat(paste(i, " "))
    }

    p.value <- empPvals(as.vector(obs), as.vector(null))

    return(list(call = match.call(),
        p.value = p.value, obs.stat = obs,
        null.stat = null))
}
