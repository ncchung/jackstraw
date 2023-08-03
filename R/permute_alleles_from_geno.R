# correctly "permutes" alleles in genotype data, to result in null data satisfying HWE (Binomial).
# - assumes input is one vector of genotypes (one SNP, all values in 0,1,2)
# - decomposes genotypes into alleles (twice as many entries, values in 0,1), permutes alleles, then regroups into genotypes
# - if there are NAs, they do not contribute to permutation, but their number is preserved and location permuted
# in other words, genotype counts are not preserved (unlike regular permutation), but allele and NA counts are preserved
permute_alleles_from_geno <- function( x ) {
    # let's not validate data (i.e. that it indeed contains exclusively data in c(0,1,2,NA) ), this should be fast and mostly used internally

    # get counts of interest
    # number of individuals (total number of alleles of both types (0 and 1) is 2*n when there's no missingness)
    n <- length( x )
    # number of individuals that are NA (there is no allele-level number of NAs to model, NA is at genotype level, that is, individuals cannot be missing just one of their alleles, its both or neither)
    nn <- sum( is.na( x ) )
    # number of complete (non-missing) individuals
    nc <- n - nn

    # handle the weird edge case that entire vector is missing (better not to die no matter what)
    # return the same vector (there is nothing to permute, it is as desired already)
    if ( nc == 0L )
        return( x )
    # now we may assume there is non-missing data
    
    # number of counted (1) alleles: heterozygotes (x=1) contribute one, homozygotes for counted allele (x=2) contribute 2
    n1 <- sum( x, na.rm = TRUE )
    # number of non-counted (0) alleles
    n0 <- 2L * nc - n1
    # here is one simple validation we should perform, if data contains values higher than 2, it could result in a negative infered value of n0, let's make sure that didn't happen or more obscure errors will just occur later
    if ( n0 < 0L )
        stop( 'The number of non-counted alleles was infered to be negative.  This happens if genotypes values sum to more than twice the number of individuals!' )
    
    # construct new genotypes from permuted approach
    # first create allele vector, ordered for convenience
    y <- c( rep.int( 1L, n1 ), rep.int( 0L, n0 ) )
    # permute this data now (it is permuting at the allele level as desired!)
    y <- sample( y )
    # by construction the length of `y` is `2 * nc`, split into two vectors each of length `nc`, then add them up, these are the new genotype values for the complete individuals!
    # Overwrite input `x` here (original is no longer needed)
    indexes <- 1L : nc
    x <- y[ indexes ] + y[ -indexes ]

    # if there's no NAs left, we're done, return what we have
    if ( nn == 0L )
        return( x )
    # else we need to add NAs to this data
    # lazy way is to just append them
    x <- c( x, rep.int( NA, nn ) )
    # then permute again!
    x <- sample( x )
    # now we are finally done, return this!
    return( x )
}
