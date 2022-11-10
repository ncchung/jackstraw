# NOTE: copy of gcatest devel version (not yet on bioconductor)
# combined delta deviances for numerical accuracy (more stable than subtracting
# separate deviances). Handles poor model fit cases (avoid errors/NAs when
# reasonable)
.delta_deviance_snp <- function(xi, pi0, pi1) {
    if (missing(xi))
        stop("Genotype vector `xi` is required!")
    if (missing(pi0))
        stop("Individual-specific allele frequency vector `pi0` is required!")
    if (missing(pi1))
        stop("Individual-specific allele frequency vector `pi1` is required!")
    if (anyNA(xi))
        stop("Genotype vector `xi` must not have NA values!")
    # try to understand NAs in pi0, pi1, and how to deal with them
    if (anyNA(pi0))
        return(NA)  # LFA/glm.fit can just fail sometimes (not an error)
    if (anyNA(pi1))
        return(NA)  # LFA/glm.fit can just fail sometimes (not an error)
    # If LFA returns all 0s or 1s for AFs, return valid numbers here (including
    # Inf or -Inf, no NAs/errors)!
    p1_all_0 <- all(pi1 == 0)
    if (all(pi0 == 0))
        return(if (p1_all_0) 0 else Inf)
    if (p1_all_0)
        return(-Inf)  # if we're here then !all(pi0 == 0)
    p1_all_1 <- all(pi1 == 1)
    if (all(pi0 == 1))
        return(if (p1_all_1) 0 else Inf)
    if (p1_all_1)
        return(-Inf)  # if we're here then !all(pi0 == 1)
    # Don't expect `p0_all_0 && p1_all_1` or other way around, those are NaN.

    # Calculate `dd` sum in parts, add factor of 2 in the end.
    if (storage.mode(xi) != "integer")
        xi <- as.integer(xi)  # for checks below, make sure xi are integers
    ind <- xi != 0L
    if (any(pi0[ind] == 0) || any(pi1[ind] == 0))
        return(NA)  # AFs==0 or 1 are bad fits and result in NA (not error).
    # Note `0 * log( 0 ) = 0` is correct limit, naive calc gives NA.
    dd <- sum(xi[ind] * log(pi1[ind]/pi0[ind]))
    # now add: ( 2 - xi ) * log( ( 1 - pi1 ) / ( 1 - pi0 ) with same rules.
    ind <- xi != 2L
    if (any(pi0[ind] == 1) || any(pi1[ind] == 1))
        return(NA)
    dd <- dd + sum((2L - xi[ind]) * log((1 - pi1[ind])/(1 - pi0[ind])))
    return(2 * dd)  # include the final factor of 2
}
