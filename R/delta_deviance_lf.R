# NOTE: copy of gcatest devel version (not yet on bioconductor)
delta_deviance_lf <- function(X, LF0, LF1) {
    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (missing(LF0))
        stop("`LF0` matrix is required!")
    if (missing(LF1))
        stop("`LF1` matrix is required!")
    # nothing can be null (`jackstraw`'s default LF0 is NULL i.e. intercept
    # only, disagreeing with `gcatest`, best to be explicit)
    if (is.null(X))
        stop("`X` cannot be null!")
    if (is.null(LF0))
        stop("`LF0` cannot be null!")
    if (is.null(LF1))
        stop("`LF1` cannot be null!")
    if (!is.matrix(X)) # check class
        stop("`X` must be a matrix!")
    n <- ncol(X) # m not used in this case
    if (nrow(LF0) != n) # check dimensions
        stop("Number of individuals in `X` and `LF0` disagrees!")
    if (nrow(LF1) != n)
        stop("Number of individuals in `X` and `LF1` disagrees!")
    if (anyNA(LF0)) # check LFs for missing values
        stop("`LF0` must not have missing values!")
    if (anyNA(LF1))
        stop("`LF1` must not have missing values!")
    # start actual processing
    devdiff <- apply(X, 1, .delta_deviance_snp_lf, LF0, LF1)
    return(devdiff)
}

