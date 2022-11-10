# NOTE: copy of gcatest devel version (not yet on bioconductor)
.delta_deviance_snp_lf <- function(xi, LF0, LF1) {
    if (missing(xi))
        stop("Genotype vector `xi` is required!")
    if (missing(LF0))
        stop("`LF0` matrix is required!")
    if (missing(LF1))
        stop("`LF1` matrix is required!")
    # nothing can be null (`jackstraw`'s default LF0 is NULL i.e. intercept
    # only, disagreeing with `gcatest`, best to be explicit)
    if (is.null(xi))
        stop("`xi` cannot be null!")
    if (is.null(LF0))
        stop("`LF0` cannot be null!")
    if (is.null(LF1))
        stop("`LF1` cannot be null!")
    # check dimensions
    n <- length(xi)
    if (nrow(LF0) != n)
        stop("Number of individuals in `xi` and `LF0` disagree!")
    if (nrow(LF1) != n)
        stop("Number of individuals in `xi` and `LF1` disagree!")
    # check LFs for missing values
    if (anyNA(LF0))
        stop("`LF0` must not have missing values!")
    if (anyNA(LF1))
        stop("`LF1` must not have missing values!")
    # remove individuals with NA genotypes (do not contribute to deviance)
    indexes_keep <- !is.na(xi)
    if (!any(indexes_keep))
        stop("All individuals at one locus were missing (unusual)!")
    if (any(!indexes_keep)) {
        xi <- xi[indexes_keep]
        LF0 <- LF0[indexes_keep, , drop = FALSE]
        LF1 <- LF1[indexes_keep, , drop = FALSE]
    }
    # perform two logistic regressions, under the null and alternative,
    # resulting in estimated allele frequencies under each model
    p0 <- lfa::af_snp(xi, LF0)
    p1 <- lfa::af_snp(xi, LF1)
    # now compute delta deviance uses a special numerically stable algorithm
    devdiff <- .delta_deviance_snp(xi, p0, p1)
    return(devdiff)
}
