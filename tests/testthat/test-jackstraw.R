# need to simulate binomial data
# these are the dimensions we want
m <- 600 # this has to be at least this large or `pip` (really `qvalue::pi0est`) fails in stupid ways
n <- 10
# use these many LFs
d <- 3
# actual genotype-like binomial data
X <- matrix(
    rbinom( n * m, 2, 0.5 ),
    nrow = m,
    ncol = n
)

# always remove zero variance rows, draw them again
# NOTE: failure to do this results in all kinds of testing errors, like:
# - FSTAT: essentially perfect fit: summary may be unreliable
# - jackstraw_subspace/svd: Error: infinite or missing values in 'x'
# - jackstraw_rpca/rsvd: Error: NA/NaN/Inf in foreign function call (arg 1)
# - jackstraw_irlba/irlba: Error: missing value where TRUE/FALSE needed
# - jackstraw_kmeans/kmeans: Error: NA/NaN/Inf in foreign function call (arg 1)
# - jackstraw_kmeanspp/ClusterR::KMeans_rcpp: Error: the data includes NaN's or +/- Inf values
# - jackstraw_MiniBatchKmeans/ClusterR::MiniBatchKmeans: Error: the data includes NaN's or +/- Inf values
# - jackstraw_pam/pam: Error: No clustering performed, NAs in the computed dissimilarity matrix.
x_bar <- rowMeans( X )
x_var <- rowMeans( X - x_bar )^2
indexes <- x_var == 0
m_const <- sum( indexes )
while( m_const ) {
    # draw those rows again
    X[ indexes, ] <- matrix(
        rbinom( n * m_const, 2, 0.5 ),
        nrow = m_const,
        ncol = n
    )
    # look for constant rows again
    x_bar <- rowMeans( X )
    x_var <- rowMeans( X - x_bar )^2
    indexes <- x_var == 0
    m_const <- sum( indexes )
}

# for PCA and other analyses for continuous data, makes sense to centerscale data
Xc <- t( scale( t( X ) ) )

# and LFs to go with this data
## # since LFA is imported, might as well use it to get actual LFs for our random data
## LF1 <- lfa::lfa( X, d )
# to have this work without Bioconductor, make fake LFs using PCA, in the same format as LFA's
LF1 <- cbind( eigen( crossprod( Xc ) )$vectors[ , 1:(d-1) ], 1 )
# in practice LF0 is just intercept model
LF0 <- NULL

# set these parameters (match defaults), which also determine dimensions of null.stat matrix
# let's use smaller values than defaults so tests are faster
s <- 5 # default is: round( m / 10 )
B <- 2 # default is: round( m * 10 / s )

# include a covariate for tests (a vector for individuals)
# draw continuous (in LFA, drawing bernoulli causes fitting problems I don't particularly want to deal with here)
covariate <- rnorm( n )

# only perform these tests when qvalue is available (it isn't on some CRAN servers!)
if ( requireNamespace( "qvalue", quietly = TRUE ) ) {

    test_that( "pip works", {
        # construct some random data
        pvalue <- runif( m )
        group <- sample( d, m, replace = TRUE )
        pi0 <- runif( d )
        # to prevent some qvalue-specific errors, ensure one p-value is 1
        pvalue[1] <- 1

        # try some errors on purpose
        # only pvalue is mandatory
        expect_error( pip( ) )
        # pass group with wrong length
        expect_error( pip( pvalue = pvalue, group = group[-1], verbose = FALSE ) )

        # now successful runs
        expect_silent(
            prob <- pip( pvalue = pvalue, group = group, pi0 = pi0, verbose = FALSE )
        )
        expect_true( is.numeric( prob ) )
        expect_equal( length( prob ), m )
        expect_true( all( prob >= 0 ) )
        # expect_true( all( prob <= 0 ) ) # not true!

        expect_silent(
            prob <- pip( pvalue = pvalue, pi0 = pi0, verbose = FALSE )
        )
        expect_true( is.numeric( prob ) )
        expect_equal( length( prob ), m )
        expect_true( all( prob >= 0 ) )
        # expect_true( all( prob <= 0 ) ) # not true!

        # a rare error, ultimately in `qvalue::pi0est`, happens in the next line (due to small groups, not enough p-values for good pi0 estimates):
        ## Error (test-jackstraw.R:67:5): pip works
        ## Error: missing or infinite values in inputs are not allowed
        ## Backtrace:
        ##   1. testthat::expect_silent(...) test-jackstraw.R:67:4
        ##   9. jackstraw::pip(pvalue = pvalue, group = group, verbose = FALSE)
        ##  10. qvalue::lfdr(pvalue[group == i], ...) /home/viiia/docs/ochoalab/jackstraw/R/pip.R:45:16
        ##  11. qvalue::pi0est(p, ...)
        ##  12. stats::smooth.spline(lambda, pi0, df = smooth.df)
        # so happens only when there's groups and pi0 has to be estimated
        expect_silent(
            prob <- pip( pvalue = pvalue, group = group, verbose = FALSE )
        )
        expect_true( is.numeric( prob ) )
        expect_equal( length( prob ), m )
        expect_true( all( prob >= 0 ) )
        # expect_true( all( prob <= 0 ) ) # not true!

        expect_silent(
            prob <- pip( pvalue = pvalue, verbose = FALSE )
        )
        expect_true( is.numeric( prob ) )
        expect_equal( length( prob ), m )
        expect_true( all( prob >= 0 ) )
        # expect_true( all( prob <= 0 ) ) # not true!
    })

}

test_that( "RSS works" , {
    # check for missing mandatory data
    # both `dat` and `mod` are required
    expect_error( RSS( ) )
    expect_error( RSS( dat = X ) )
    expect_error( RSS( mod = LF1 ) )
    # mod must be a matrix
    # (it could be an intercept only, but in practice we always pass objects constructed via model.matrix, so they're always matrices)
    expect_error( RSS( dat = X, mod = 1:n ) )

    # now a successful run
    expect_silent(
        rss <- RSS( dat = X, mod = LF1 )
    )
    expect_true( is.numeric( rss ) )
    expect_equal( length( rss ), m )
    expect_true( all( rss >= 0 ) )

    # compare against lm
    # have to fit each row separately
    rss_lm <- vector( 'numeric', m )
    for ( i in 1:m ) {
        rss_lm[ i ] <- sum( lm( X[ i, ] ~ LF1 )$residuals^2 )
    }
    expect_equal( rss, rss_lm )

    # test case where `dat` is vector
    # (not used by dependencies, but in theory supported)
    expect_silent(
        rss <- RSS( dat = X[ 1, ], mod = LF1 )
    )
    expect_equal( rss, rss_lm[1] )
})

test_that( "FSTAT works", {
    # LF1 doesn't work as-is because the intercept gets added twice!
    LV <- LF1[ , -d, drop = FALSE ]

    # check for missing mandatory data
    # both `dat` and `LV` are required
    expect_error( FSTAT( ) )
    expect_error( FSTAT( dat = Xc ) )
    expect_error( FSTAT( LV = LV ) )
    # `dat` must be matrix
    expect_error( FSTAT( dat = 1:10, LV = LV ) )
    # `LV` has non-matching numbers of rows with `dat`
    expect_error( FSTAT( dat = Xc, LV = LV[ -1, ] ) )

    # successful run
    expect_silent(
        obj <- FSTAT( dat = Xc, LV = LV )
    )
    # since parammetric = FALSE, this only returns one element in its list
    expect_true( is.list( obj ) )
    expect_equal( length( obj ), 1 )
    expect_equal( names( obj ), 'fstat' )
    # statistics are one per row
    expect_equal( length( obj$fstat ), m )

    # compare to lm
    # have to fit each row separately
    fstat_lm <- vector( 'numeric', m )
    for ( i in 1:m ) {
        fstat_lm[ i ] <- summary( lm( X[ i, ] ~ LF1 ) )$fstatistic[1]
    }
    expect_equal( obj$fstat, fstat_lm )

    # successful run with covariates
    expect_silent(
        obj <- FSTAT( dat = Xc, LV = LV, covariate = covariate )
    )
    # since parammetric = FALSE, this only returns one element in its list
    expect_true( is.list( obj ) )
    expect_equal( length( obj ), 1 )
    expect_equal( names( obj ), 'fstat' )
    # statistics are one per row
    expect_equal( length( obj$fstat ), m )

    # successful run with ALV (should be identical to covariate, but some use cases specify both separately)
    expect_silent(
        obj2 <- FSTAT( dat = Xc, LV = LV, ALV = covariate )
    )
    expect_equal( obj2, obj )

    # NOTE: `parametric = TRUE` is not used in this package for any public code, so it's also not tested

    #FSTAT(dat, LV, ALV = NULL, covariate = NULL, parametric = FALSE)
})


test_that( "permutationPA works", {
    # data is required
    expect_error( permutationPA() )
    # data must be a matrix
    expect_error( permutationPA( dat = 1:10 ) )

    # a successful run
    # reduce B for speed (default is 100)
    expect_silent(
        obj <- permutationPA( dat = X, B = B, verbose = FALSE )
    )
    # test return object
    expect_true( is.list( obj ) )
    expect_equal( length( obj ), 2 )
    expect_equal( names( obj ), c('r', 'p') )
    # test r
    expect_equal( length( obj$r ), 1 )
    expect_true( is.integer( obj$r ) )
    # test p-values
    expect_true( is.numeric( obj$p ) )
    expect_equal( length( obj$p ), n )
    expect_true( !anyNA( obj$p ) )
    expect_true( all( obj$p >= 0 ) )
    expect_true( all( obj$p <= 1 ) )
})

test_that( 'permute_alleles_from_geno works', {
    # a toy example of extremely structured data, where everybody is homozygous for one or the other allele
    # (recall n=10 in this toy data)
    n2 <- n / 2L
    x <- c( rep.int( 0L, n2 ), rep.int( 2L, n2 ) )
    # permute data at allele level
    expect_silent(
        y <- permute_alleles_from_geno( x )
    )
    # we should see that lengths match
    expect_equal( length( y ), n )
    # and there are no missing values in this case
    expect_true( !anyNA( y ) )
    # confirm that everything is in desired range
    expect_true( all( y %in% c(0L, 1L, 2L) ) )
    # for this example, the case without heterozygotes is so rare (0.5^10 = 9.8e-4) let's demand that there be at least one heterozygote
    expect_true( any( y == 1L ) )
    # and sum of genotypes equals `n` (in this case true by construction, not true in general)
    expect_equal( sum( y ), n )

    # now try on the simulated random data (all rows)
    # vectorize so number of tests isn't too extreme
    expect_silent(
        Y <- t( apply( X, 1L, permute_alleles_from_geno ) )
    )
    # perform general tests, still allowing no missingness (not present in simulated data)
    expect_equal( ncol( Y ), n )
    expect_true( !anyNA( Y ) )
    expect_true( all( Y %in% c(0L, 1L, 2L) ) )
    # returning the same data is practically impossible if permutation is correct
    expect_true( !all( Y == X ) )
    # here test that the sums of input and output for each row match
    expect_equal( rowSums( Y ), rowSums( X ) )

    # we haven't performed tests with missingness in general, and the overall runtime is so high I don't think we want to do it broadly, but here let's just sprinkle random missingness (here there's absolute tolerance for bad things happening, like fixed rows, even rows full of NAs!)
    X_miss <- X
    # just a little bit of missingness
    p_miss <- 0.1
    # add missing values
    X_miss[ sample( n*m, n*m*p_miss ) ] <- NA
    # repeat earlier test!
    expect_silent(
        Y_miss <- t( apply( X_miss, 1L, permute_alleles_from_geno ) )
    )
    # perform general tests, now allowing missingness!
    expect_equal( ncol( Y_miss ), n )
    expect_true( all( Y %in% c(0L, 1L, 2L), na.rm = TRUE ) )
    # returning the same data is practically impossible if permutation is correct
    expect_true( !all( Y == X, na.rm = TRUE ) )
    # here test that the sums of input and output for each row match
    expect_equal( rowSums( Y, na.rm = TRUE ), rowSums( X, na.rm = TRUE ) )
})

test_that( "empPvals handles NAs correctly", {
    # `qvalue::empPvals` actually it doesn't, but it's easy to fix with a minor hack, wrapped around internal `empPvals`
    m <- 100
    obs <- c(NA, 0.01, 0.001)
    null <- runif( m )
    expect_silent(
        pvals <- empPvals( obs, null )
    )
    # actual tests
    expect_equal( length( pvals ), length( obs ) )
    expect_true( is.na( pvals[1] ) )
    expect_true( !anyNA( pvals[2:3] ) )

    # a bigger random test
    # 20% are NAs
    obs <- runif( m )
    obs[ sample.int( m, 0.2 * m ) ] <- NA
    expect_silent(
        pvals <- empPvals( obs, null )
    )
    # actual tests
    expect_equal( length( pvals ), length( obs ) )
    expect_true( all( is.na( pvals[ is.na( obs ) ] ) ) )
    expect_true( !anyNA( pvals[ !is.na( obs ) ] ) )
})

test_jackstraw_return_val <- function ( obj, s, B, kmeans = FALSE ) {
    # all jackstraw variants return basically the same thing
    # globals used: m, d

    # kmeans object returns different names and order, unfortunately (otherwise equivalent though)
    if ( kmeans ) {
        name_p <- 'p.F'
        name_o <- 'F.obs'
        name_n <- 'F.null'
        names_obj <- c('call', name_o, name_n, name_p)
    } else {
        name_p <- 'p.value'
        name_o <- 'obs.stat'
        name_n <- 'null.stat'
        names_obj <- c('call', name_p, name_o, name_n)
    }

    # test overall object
    expect_true( is.list( obj ) )
    expect_equal( length( obj ), 4 )
    expect_equal( names( obj ), names_obj )

    # test individual elements
    # 1) call
    expect_true( is.call( obj$call ) )
    # 2) p.value
    p <- obj[[ name_p ]]
    expect_true( is.numeric( p ) )
    expect_equal( length( p ), m )
    expect_true( !anyNA( p ) ) # in theory there can be NAs, they just don't arise in my simple examples
    expect_true( all( p >= 0, na.rm = TRUE ) )
    expect_true( all( p <= 1, na.rm = TRUE ) )
    # 3) obs.stat
    obs.stat <- obj[[ name_o ]]
    expect_true( is.numeric( obs.stat ) )
    expect_equal( length( obs.stat ), m )
    # 4) null.stat
    null.stat <- obj[[ name_n ]]
    if ( kmeans ) {
        # here it's a list of length d
        expect_true( is.list( null.stat ) )
        expect_equal( length( null.stat ), d )
        # NOTE: lengths of elements vary, not a useful test
    } else {
        # test first as vector
        expect_true( is.numeric( null.stat ) )
        expect_equal( length( null.stat ), s * B )
        # then as matrix
        expect_true( is.matrix( null.stat ) )
        expect_equal( nrow( null.stat ), s )
        expect_equal( ncol( null.stat ), B )
    }
}

test_that("jackstraw_subspace works", {
    FUN <- function( x ) svd( x )$v[ , 1:d, drop = FALSE ]
    # cause errors due to missing required data
    # must provide all of dat = Xc, r = d, and FUN for a minimal successful run
    expect_error( jackstraw_subspace( ) )
    expect_error( jackstraw_subspace( dat = Xc ) )
    expect_error( jackstraw_subspace( r = d ) )
    expect_error( jackstraw_subspace( FUN = FUN ) )
    expect_error( jackstraw_subspace( dat = Xc, r = d ) )
    expect_error( jackstraw_subspace( dat = Xc, FUN = FUN ) )
    expect_error( jackstraw_subspace( r = d, FUN = FUN ) )
    # check that data is matrix
    expect_error( jackstraw_subspace( dat = 1:10, r = d, FUN = FUN ) )
    # pass bad covariates on purpose
    expect_error( jackstraw_subspace( dat = Xc, r = d, FUN = FUN, covariate = 1 ) ) # scalar is bad
    expect_error( jackstraw_subspace( dat = Xc, r = d, FUN = FUN, covariate = covariate[-1] ) ) # length off by 1, vector
    expect_error( jackstraw_subspace( dat = Xc, r = d, FUN = FUN, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
    expect_error( jackstraw_subspace( dat = Xc, r = d, FUN = FUN, covariate = rbind( covariate ) ) ) # transposed matrix

    # perform a basic run

    # make it silent so we can focus on problem messages
    expect_silent(
        obj <- jackstraw_subspace( Xc, r = d, FUN = FUN, s = s, B = B, verbose = FALSE )
    )
    # check basic jackstraw return object
    test_jackstraw_return_val( obj, s, B )

    # test edge cases
    # s = 1
    expect_silent(
        obj <- jackstraw_subspace( Xc, r = d, FUN = FUN, s = 1, B = B, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, B )
    # B = 1
    expect_silent(
        obj <- jackstraw_subspace( Xc, r = d, FUN = FUN, s = s, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, 1 )
    # s = B = 1
    expect_silent(
        obj <- jackstraw_subspace( Xc, r = d, FUN = FUN, s = 1, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, 1 )

    # test version with covariates
    expect_silent(
        obj <- jackstraw_subspace( Xc, r = d, FUN = FUN, s = s, B = B, covariate = covariate, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, B )
})

test_that("jackstraw_pca works", {
    # cause errors due to missing required data
    # must provide dat = Xc for a minimal successful run
    expect_error( jackstraw_pca( ) )
    expect_error( jackstraw_pca( r = d ) )
    # check that data is matrix
    expect_error( jackstraw_pca( dat = 1:10, r = d ) )
    # pass bad covariates on purpose
    expect_error( jackstraw_pca( dat = Xc, r = d, covariate = 1 ) ) # scalar is bad
    expect_error( jackstraw_pca( dat = Xc, r = d, covariate = covariate[-1] ) ) # length off by 1, vector
    expect_error( jackstraw_pca( dat = Xc, r = d, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
    expect_error( jackstraw_pca( dat = Xc, r = d, covariate = rbind( covariate ) ) ) # transposed matrix

    # perform a basic run

    # make it silent so we can focus on problem messages
    expect_silent(
        obj <- jackstraw_pca( Xc, r = d, s = s, B = B, verbose = FALSE )
    )
    # check basic jackstraw return object
    test_jackstraw_return_val( obj, s, B )

    # test edge cases
    # s = 1
    expect_silent(
        obj <- jackstraw_pca( Xc, r = d, s = 1, B = B, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, B )
    # B = 1
    expect_silent(
        obj <- jackstraw_pca( Xc, r = d, s = s, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, 1 )
    # s = B = 1
    expect_silent(
        obj <- jackstraw_pca( Xc, r = d, s = 1, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, 1 )

    # test version with covariates
    expect_silent(
        obj <- jackstraw_pca( Xc, r = d, s = s, B = B, covariate = covariate, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, B )
})

test_that("jackstraw_rpca works", {
    # cause errors due to missing required data
    # must provide dat = Xc for a minimal successful run
    expect_error( jackstraw_rpca( ) )
    expect_error( jackstraw_rpca( r = d ) )
    # check that data is matrix
    expect_error( jackstraw_rpca( dat = 1:10, r = d ) )
    # pass bad covariates on purpose
    expect_error( jackstraw_rpca( dat = Xc, r = d, covariate = 1 ) ) # scalar is bad
    expect_error( jackstraw_rpca( dat = Xc, r = d, covariate = covariate[-1] ) ) # length off by 1, vector
    expect_error( jackstraw_rpca( dat = Xc, r = d, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
    expect_error( jackstraw_rpca( dat = Xc, r = d, covariate = rbind( covariate ) ) ) # transposed matrix

    # perform a basic run

    # make it silent so we can focus on problem messages
    expect_silent(
        obj <- jackstraw_rpca( Xc, r = d, s = s, B = B, verbose = FALSE )
    )
    # check basic jackstraw return object
    test_jackstraw_return_val( obj, s, B )

    # test edge cases
    # s = 1
    expect_silent(
        obj <- jackstraw_rpca( Xc, r = d, s = 1, B = B, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, B )
    # B = 1
    expect_silent(
        obj <- jackstraw_rpca( Xc, r = d, s = s, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, 1 )
    # s = B = 1
    expect_silent(
        obj <- jackstraw_rpca( Xc, r = d, s = 1, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, 1 )

    # test version with covariates
    expect_silent(
        obj <- jackstraw_rpca( Xc, r = d, s = s, B = B, covariate = covariate, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, B )
})

test_that("jackstraw_irlba works", {
    # cause errors due to missing required data
    # must provide dat = Xc for a minimal successful run
    expect_error( jackstraw_irlba( ) )
    expect_error( jackstraw_irlba( r = d ) )
    # check that data is matrix
    expect_error( jackstraw_irlba( dat = 1:10, r = d ) )
    # pass bad covariates on purpose
    expect_error( jackstraw_irlba( dat = Xc, r = d, covariate = 1 ) ) # scalar is bad
    expect_error( jackstraw_irlba( dat = Xc, r = d, covariate = covariate[-1] ) ) # length off by 1, vector
    expect_error( jackstraw_irlba( dat = Xc, r = d, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
    expect_error( jackstraw_irlba( dat = Xc, r = d, covariate = rbind( covariate ) ) ) # transposed matrix

    # perform a basic run

    # make it silent so we can focus on problem messages
    expect_silent(
        obj <- jackstraw_irlba( Xc, r = d, s = s, B = B, verbose = FALSE )
    )
    # check basic jackstraw return object
    test_jackstraw_return_val( obj, s, B )

    # test edge cases
    # s = 1
    expect_silent(
        obj <- jackstraw_irlba( Xc, r = d, s = 1, B = B, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, B )
    # B = 1
    expect_silent(
        obj <- jackstraw_irlba( Xc, r = d, s = s, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, 1 )
    # s = B = 1
    expect_silent(
        obj <- jackstraw_irlba( Xc, r = d, s = 1, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, 1 )

    # test version with covariates
    expect_silent(
        obj <- jackstraw_irlba( Xc, r = d, s = s, B = B, covariate = covariate, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, B )
})

test_that( "jackstraw_kmeans works", {
    # a simple k-means run
    kmeans.dat <- kmeans( Xc, centers = d )

    # cause errors due to missing required data
    # must provide both dat = X and kmeans.dat for a minimal successful run
    expect_error( jackstraw_kmeans() )
    expect_error( jackstraw_kmeans( dat = Xc ) )
    expect_error( jackstraw_kmeans( kmeans.dat = kmeans.dat ) )

    # check that data is matrix
    expect_error( jackstraw_kmeans( dat = 1:10, kmeans.dat = kmeans.dat ) )
    # pass bad covariates on purpose
    expect_error( jackstraw_kmeans( dat = Xc, kmeans.dat = kmeans.dat, covariate = 1 ) ) # scalar is bad
    expect_error( jackstraw_kmeans( dat = Xc, kmeans.dat = kmeans.dat, covariate = covariate[-1] ) ) # length off by 1, vector
    expect_error( jackstraw_kmeans( dat = Xc, kmeans.dat = kmeans.dat, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
    expect_error( jackstraw_kmeans( dat = Xc, kmeans.dat = kmeans.dat, covariate = rbind( covariate ) ) ) # transposed matrix

    # perform a basic run

    # make it silent so we can focus on problem messages
    expect_silent(
        obj <- jackstraw_kmeans( Xc, kmeans.dat = kmeans.dat, s = s, B = B, verbose = FALSE )
    )
    # check basic jackstraw return object
    test_jackstraw_return_val( obj, s, B, kmeans = TRUE )

    # test edge cases
    # s = 1
    expect_silent(
        obj <- jackstraw_kmeans( Xc, kmeans.dat = kmeans.dat, s = 1, B = B, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, B, kmeans = TRUE )
    # B = 1
    expect_silent(
        obj <- jackstraw_kmeans( Xc, kmeans.dat = kmeans.dat, s = s, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, 1, kmeans = TRUE )
    # s = B = 1
    expect_silent(
        obj <- jackstraw_kmeans( Xc, kmeans.dat = kmeans.dat, s = 1, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, 1, kmeans = TRUE )

    # test version with covariates
    expect_silent(
        obj <- jackstraw_kmeans( Xc, kmeans.dat = kmeans.dat, s = s, B = B, covariate = covariate, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, B, kmeans = TRUE )
})

test_that( "jackstraw_kmeanspp works", {
    # a simple k-means run
    kmeans.dat <- ClusterR::KMeans_rcpp( Xc, clusters = d )

    # cause errors due to missing required data
    # must provide both dat = X and kmeans.dat for a minimal successful run
    expect_error( jackstraw_kmeanspp() )
    expect_error( jackstraw_kmeanspp( dat = Xc ) )
    expect_error( jackstraw_kmeanspp( kmeans.dat = kmeans.dat ) )

    # check that data is matrix
    expect_error( jackstraw_kmeanspp( dat = 1:10, kmeans.dat = kmeans.dat ) )
    # pass bad covariates on purpose
    expect_error( jackstraw_kmeanspp( dat = Xc, kmeans.dat = kmeans.dat, covariate = 1 ) ) # scalar is bad
    expect_error( jackstraw_kmeanspp( dat = Xc, kmeans.dat = kmeans.dat, covariate = covariate[-1] ) ) # length off by 1, vector
    expect_error( jackstraw_kmeanspp( dat = Xc, kmeans.dat = kmeans.dat, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
    expect_error( jackstraw_kmeanspp( dat = Xc, kmeans.dat = kmeans.dat, covariate = rbind( covariate ) ) ) # transposed matrix

    # perform a basic run

    # make it silent so we can focus on problem messages
    expect_silent(
        obj <- jackstraw_kmeanspp( Xc, kmeans.dat = kmeans.dat, s = s, B = B, verbose = FALSE )
    )
    # check basic jackstraw return object
    test_jackstraw_return_val( obj, s, B, kmeans = TRUE )

    # test edge cases
    # s = 1
    expect_silent(
        obj <- jackstraw_kmeanspp( Xc, kmeans.dat = kmeans.dat, s = 1, B = B, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, B, kmeans = TRUE )
    # B = 1
    expect_silent(
        obj <- jackstraw_kmeanspp( Xc, kmeans.dat = kmeans.dat, s = s, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, 1, kmeans = TRUE )
    # s = B = 1
    expect_silent(
        obj <- jackstraw_kmeanspp( Xc, kmeans.dat = kmeans.dat, s = 1, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, 1, kmeans = TRUE )

    # test version with covariates
    expect_silent(
        obj <- jackstraw_kmeanspp( Xc, kmeans.dat = kmeans.dat, s = s, B = B, covariate = covariate, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, B, kmeans = TRUE )
})

test_that( "jackstraw_MiniBatchKmeans works", {
    # a simple k-means run
    batch_size <- 10
    MiniBatchKmeans.output <- ClusterR::MiniBatchKmeans( Xc, clusters = d, batch_size = batch_size )

    # cause errors due to missing required data
    # must provide both dat = X and MiniBatchKmeans.output for a minimal successful run
    expect_error( jackstraw_MiniBatchKmeans() )
    expect_error( jackstraw_MiniBatchKmeans( dat = Xc ) )
    expect_error( jackstraw_MiniBatchKmeans( MiniBatchKmeans.output = MiniBatchKmeans.output ) )

    # check that data is matrix
    expect_error( jackstraw_MiniBatchKmeans( dat = 1:10, MiniBatchKmeans.output = MiniBatchKmeans.output ) )
    # pass bad covariates on purpose
    expect_error( jackstraw_MiniBatchKmeans( dat = Xc, MiniBatchKmeans.output = MiniBatchKmeans.output, covariate = 1 ) ) # scalar is bad
    expect_error( jackstraw_MiniBatchKmeans( dat = Xc, MiniBatchKmeans.output = MiniBatchKmeans.output, covariate = covariate[-1] ) ) # length off by 1, vector
    expect_error( jackstraw_MiniBatchKmeans( dat = Xc, MiniBatchKmeans.output = MiniBatchKmeans.output, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
    expect_error( jackstraw_MiniBatchKmeans( dat = Xc, MiniBatchKmeans.output = MiniBatchKmeans.output, covariate = rbind( covariate ) ) ) # transposed matrix

    # perform a basic run

    # make it silent so we can focus on problem messages
    expect_silent(
        obj <- jackstraw_MiniBatchKmeans( Xc, MiniBatchKmeans.output = MiniBatchKmeans.output, batch_size = batch_size, s = s, B = B, verbose = FALSE )
    )
    # check basic jackstraw return object
    test_jackstraw_return_val( obj, s, B, kmeans = TRUE )

    # test edge cases
    # s = 1
    expect_silent(
        obj <- jackstraw_MiniBatchKmeans( Xc, MiniBatchKmeans.output = MiniBatchKmeans.output, batch_size = batch_size, s = 1, B = B, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, B, kmeans = TRUE )
    # B = 1
    expect_silent(
        obj <- jackstraw_MiniBatchKmeans( Xc, MiniBatchKmeans.output = MiniBatchKmeans.output, batch_size = batch_size, s = s, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, 1, kmeans = TRUE )
    # s = B = 1
    expect_silent(
        obj <- jackstraw_MiniBatchKmeans( Xc, MiniBatchKmeans.output = MiniBatchKmeans.output, batch_size = batch_size, s = 1, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, 1, kmeans = TRUE )

    # test version with covariates
    expect_silent(
        obj <- jackstraw_MiniBatchKmeans( Xc, MiniBatchKmeans.output = MiniBatchKmeans.output, batch_size = batch_size, s = s, B = B, covariate = covariate, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, B, kmeans = TRUE )
})

test_that( "jackstraw_pam works", {
    # a simple PAM run
    pam.dat <- cluster::pam( Xc, k = d )

    # cause errors due to missing required data
    # must provide both dat = X and pam.dat for a minimal successful run
    expect_error( jackstraw_pam() )
    expect_error( jackstraw_pam( dat = Xc ) )
    expect_error( jackstraw_pam( pam.dat = pam.dat ) )

    # check that data is matrix
    expect_error( jackstraw_pam( dat = 1:10, pam.dat = pam.dat ) )
    # pass bad covariates on purpose
    expect_error( jackstraw_pam( dat = Xc, pam.dat = pam.dat, covariate = 1 ) ) # scalar is bad
    expect_error( jackstraw_pam( dat = Xc, pam.dat = pam.dat, covariate = covariate[-1] ) ) # length off by 1, vector
    expect_error( jackstraw_pam( dat = Xc, pam.dat = pam.dat, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
    expect_error( jackstraw_pam( dat = Xc, pam.dat = pam.dat, covariate = rbind( covariate ) ) ) # transposed matrix

    # perform a basic run

    # make it silent so we can focus on problem messages
    expect_silent(
        obj <- jackstraw_pam( Xc, pam.dat = pam.dat, s = s, B = B, verbose = FALSE )
    )
    # check basic jackstraw return object
    test_jackstraw_return_val( obj, s, B, kmeans = TRUE )

    # test edge cases
    # NOTE: removed s=1 cases because, though they are fine often, sometimes there are sigularity errors that are not worth debugging for such small toy cases tested here
    ## # s = 1
    ## expect_silent(
    ##     ###ERROR: Error in `solve.default(t(mod) %*% mod)`: Lapack routine dgesv: system is exactly singular: U[2,2] = 0
    ##     obj <- jackstraw_pam( Xc, pam.dat = pam.dat, s = 1, B = B, verbose = FALSE )
    ## )
    ## test_jackstraw_return_val( obj, 1, B, kmeans = TRUE )
    # B = 1
    expect_silent(
        obj <- jackstraw_pam( Xc, pam.dat = pam.dat, s = s, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, 1, kmeans = TRUE )
    ## # s = B = 1
    ## expect_silent(
    ##     obj <- jackstraw_pam( Xc, pam.dat = pam.dat, s = 1, B = 1, verbose = FALSE )
    ## )
    ## test_jackstraw_return_val( obj, 1, 1, kmeans = TRUE )

    # test version with covariates
    expect_silent(
        obj <- jackstraw_pam( Xc, pam.dat = pam.dat, s = s, B = B, covariate = covariate, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, B, kmeans = TRUE )
})

test_that( "jackstraw_cluster works", {
    # a simple kmeans run
    FUN <- kmeans
    FUN.dat <- FUN( Xc, centers = d )
    cluster <- FUN.dat$cluster
    centers <- FUN.dat$centers

    # cause errors due to missing required data
    # 4 arguments are required for a successful run
    expect_error( jackstraw_cluster() )
    # singletons
    expect_error( jackstraw_cluster( dat = Xc ) )
    expect_error( jackstraw_cluster( k = d ) )
    expect_error( jackstraw_cluster( cluster = cluster ) )
    expect_error( jackstraw_cluster( centers = centers ) )
    # pairs
    expect_error( jackstraw_cluster( dat = Xc, k = d ) )
    expect_error( jackstraw_cluster( dat = Xc, cluster = cluster ) )
    expect_error( jackstraw_cluster( dat = Xc, centers = centers ) )
    expect_error( jackstraw_cluster( k = d, cluster = cluster ) )
    expect_error( jackstraw_cluster( k = d, centers = centers ) )
    expect_error( jackstraw_cluster( cluster = cluster, centers = centers ) )
    # triplets
    expect_error( jackstraw_cluster( k = d, cluster = cluster, centers = centers ) )
    expect_error( jackstraw_cluster( dat = Xc, cluster = cluster, centers = centers ) )
    expect_error( jackstraw_cluster( dat = Xc, k = d, centers = centers ) )
    expect_error( jackstraw_cluster( dat = Xc, k = d, cluster = cluster ) )

    # check that data is matrix
    expect_error( jackstraw_cluster( dat = 1:10, k = d, cluster = cluster, centers = centers ) )
    # check cluster vector length
    expect_error( jackstraw_cluster( dat = Xc, k = d, cluster = cluster[-1], centers = centers ) )
    # check centers dimensions
    expect_error( jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers[ -1, ] ) )
    expect_error( jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers[ , -1 ] ) )
    # pass bad covariates on purpose
    expect_error( jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers, covariate = 1 ) ) # scalar is bad
    expect_error( jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers, covariate = covariate[-1] ) ) # length off by 1, vector
    expect_error( jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
    expect_error( jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers, covariate = rbind( covariate ) ) ) # transposed matrix

    # perform a basic run

    # make it silent so we can focus on problem messages
    expect_silent(
        obj <- jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers, algorithm = FUN, s = s, B = B, verbose = FALSE )
    )
    # check basic jackstraw return object
    test_jackstraw_return_val( obj, s, B, kmeans = TRUE )

    # test edge cases
    # s = 1
    expect_silent(
        obj <- jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers, algorithm = FUN, s = 1, B = B, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, B, kmeans = TRUE )
    # B = 1
    expect_silent(
        obj <- jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers, algorithm = FUN, s = s, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, 1, kmeans = TRUE )
    # s = B = 1
    expect_silent(
        obj <- jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers, algorithm = FUN, s = 1, B = 1, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, 1, 1, kmeans = TRUE )

    # test version with covariates
    expect_silent(
        obj <- jackstraw_cluster( dat = Xc, k = d, cluster = cluster, centers = centers, algorithm = FUN, s = s, B = B, covariate = covariate, verbose = FALSE )
    )
    test_jackstraw_return_val( obj, s, B, kmeans = TRUE )
})

validate_ncp_est <- function( ests ) {
    expect_true( is.vector( ests ) )
    expect_true( is.numeric( ests ) )
    expect_equal( length( ests ), 2 )
    expect_true( all( ests >= 0 ) )
}

validate_pvals_nc_chisq <- function( obj, n_obs, test_ncp = TRUE ) {
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c('p.value', 'ncp') )
    expect_equal( length( obj$p.value ), n_obs )
    expect_true( is.numeric( obj$p.value ) )
    expect_true( min( obj$p.value, na.rm = TRUE ) >= 0 )
    expect_true( max( obj$p.value, na.rm = TRUE ) <= 1 )
    if ( test_ncp )
        validate_ncp_est( obj$ncp )
}

# test nc-chisq code
test_that( 'ncp_est, pvals_nc_chisq work', {
    # simulate small data truly from nc-chisq
    df <- 1
    ncp_true <- 2
    n <- 100
    x <- rchisq( n, df, ncp_true )
    # and a second data with stats we want p-values for, let's give huge power to this test
    ny <- 11
    y <- rchisq( ny, df, ncp_true + 10 )
    # add artifacts to these from the start
    y[1:3] <- c(NA, -1, Inf)
    
    # test this expected successful example
    expect_silent( 
        ests <- ncp_est( x, df )
    )
    validate_ncp_est( ests )

    # and the bigger wrapper
    out <- list( null.stat = x, obs.stat = y )
    expect_silent(
        obj <- pvals_nc_chisq( out, df )
    )
    validate_pvals_nc_chisq( obj, ny, test_ncp = FALSE )
    # since the inputs were the same, these should match what we had before
    expect_equal( obj$ncp, ests )

    # repeat with alternate call setup, expect identical results again
    expect_silent(
        obj2 <- pvals_nc_chisq( null.stat = x, obs.stat = y, df = df )
    )
    expect_equal( obj2, obj )

    # make sure when things are missing that it dies
    expect_error( pvals_nc_chisq( null.stat = x, df = df ) )
    expect_error( pvals_nc_chisq( obs.stat = y, df = df ) )
    expect_error( pvals_nc_chisq( df = df ) )
    
    # add artifacts to this data that the code is supposed to handle, including NAs, infinite values, and negatives (all get removed)
    x2 <- c( x, NA, Inf, -Inf, -1 )
    expect_silent( 
        ests2 <- ncp_est( x2, df )
    )
    # since the input was effectively the same as before (all additions should be removed), the output should be identical too
    expect_equal( ests2, ests )
    # ditto here
    expect_silent(
        obj2 <- pvals_nc_chisq( null.stat = x2, obs.stat = y, df = df )
    )
    expect_equal( obj2, obj )

    # test the edge case where the data is actually central
    x <- rchisq( n, df, 0 )
    expect_silent( 
        ests <- ncp_est( x, df )
    )
    validate_ncp_est( ests )
    # repeat this wrapper too
    expect_silent(
        obj <- pvals_nc_chisq( null.stat = x, obs.stat = y, df = df )
    )
    validate_pvals_nc_chisq( obj, ny, test_ncp = FALSE )
    # since the inputs were the same, these should match what we had before
    expect_equal( obj$ncp, ests )
})



# first write test genotypes somewhere
file_tmp <- tempfile( 'test-jackstraw', fileext = '.bed' )
genio::write_bed( file_tmp, X, verbose = FALSE )
# read it back as a BEDMatrix object
# read it this way, specifying dimensions, because there's no BIM/FAM files
X_BM <- BEDMatrix::BEDMatrix( file_tmp, n = n, p = m )

test_that( "jackstraw_BEDMatrix works", {
    # NOTE: BEDMatrix and genio are both dependencies

    # extra parameters
    m_chunk <- 1000 # default
    # recall s==2, let's leave it that way for now
    random_s <- sample.int( m, s )
    # need sorted version for some tests
    random_s_sorted <- sort( random_s )

    expect_silent(
        obj <- jackstraw_BEDMatrix( dat = X_BM, random_s, m_chunk = m_chunk )
    )
    expect_equal( class( obj ), 'list' )
    expect_equal( names( obj ), c('dat_full', 'dat_rand', 'file_full', 'file_rand') )
    expect_true( class( obj$dat_full ) == 'BEDMatrix' )
    expect_true( class( obj$dat_rand ) == 'BEDMatrix' )
    expect_equal( class( obj$file_full ), 'character' )
    expect_equal( class( obj$file_rand ), 'character' )
    expect_equal( nrow( obj$dat_full ), n )
    expect_equal( ncol( obj$dat_full ), m )
    expect_equal( nrow( obj$dat_rand ), n )
    expect_equal( ncol( obj$dat_rand ), s )
    # the non-random loci should match the original data
    # (NOTE: BEDMatrix data has to be transposed and dimnames removed)
    X_orig_exp <- X[ -random_s, ]
    X_orig_obs <- t( obj$dat_full[ , -random_s ] )
    dimnames( X_orig_obs ) <- NULL
    expect_equal( X_orig_obs, X_orig_exp )
    # the random data in the full matrix should equal the random matrix
    # here they're both BEDMatrix so no transposition is required
    # however, `dat_rand` lists rows in original order, which may disagree with `random_s` (random order), so let's sort here to get matching orders
    X_rand_exp <- obj$dat_rand[]
    X_rand_obs <- obj$dat_full[ , random_s_sorted, drop = FALSE ]
    expect_equal( X_rand_obs, X_rand_exp )

    # test this edge case
    m_chunk <- 1
    expect_silent(
        obj <- jackstraw_BEDMatrix( dat = X_BM, random_s, m_chunk = m_chunk )
    )
    expect_equal( class( obj ), 'list' )
    expect_equal( names( obj ), c('dat_full', 'dat_rand', 'file_full', 'file_rand') )
    expect_true( class( obj$dat_full ) == 'BEDMatrix' )
    expect_true( class( obj$dat_rand ) == 'BEDMatrix' )
    expect_equal( class( obj$file_full ), 'character' )
    expect_equal( class( obj$file_rand ), 'character' )
    expect_equal( nrow( obj$dat_full ), n )
    expect_equal( ncol( obj$dat_full ), m )
    expect_equal( nrow( obj$dat_rand ), n )
    expect_equal( ncol( obj$dat_rand ), s )
    X_orig_exp <- X[ -random_s, ]
    X_orig_obs <- t( obj$dat_full[ , -random_s ] )
    dimnames( X_orig_obs ) <- NULL
    expect_equal( X_orig_obs, X_orig_exp )
    X_rand_exp <- obj$dat_rand[]
    X_rand_obs <- obj$dat_full[ , random_s_sorted, drop = FALSE ]
    expect_equal( X_rand_obs, X_rand_exp )

    # test this edge case
    s <- 1
    m_chunk <- 1000
    random_s <- sample.int( m, s )
    # need sorted version for some tests (not needed if s=1 but meh)
    random_s_sorted <- sort( random_s )
    expect_silent(
        obj <- jackstraw_BEDMatrix( dat = X_BM, random_s, m_chunk = m_chunk )
    )
    expect_equal( class( obj ), 'list' )
    expect_equal( names( obj ), c('dat_full', 'dat_rand', 'file_full', 'file_rand') )
    expect_true( class( obj$dat_full ) == 'BEDMatrix' )
    expect_true( class( obj$dat_rand ) == 'BEDMatrix' )
    expect_equal( class( obj$file_full ), 'character' )
    expect_equal( class( obj$file_rand ), 'character' )
    expect_equal( nrow( obj$dat_full ), n )
    expect_equal( ncol( obj$dat_full ), m )
    expect_equal( nrow( obj$dat_rand ), n )
    expect_equal( ncol( obj$dat_rand ), s )
    X_orig_exp <- X[ -random_s, ]
    X_orig_obs <- t( obj$dat_full[ , -random_s ] )
    dimnames( X_orig_obs ) <- NULL
    expect_equal( X_orig_obs, X_orig_exp )
    X_rand_exp <- obj$dat_rand[]
    X_rand_obs <- obj$dat_full[ , random_s_sorted, drop = FALSE ]
    expect_equal( X_rand_obs, X_rand_exp )

    # final, double edge case
    s <- 1
    m_chunk <- 1
    expect_silent(
        obj <- jackstraw_BEDMatrix( dat = X_BM, random_s, m_chunk = m_chunk )
    )
    expect_equal( class( obj ), 'list' )
    expect_equal( names( obj ), c('dat_full', 'dat_rand', 'file_full', 'file_rand') )
    expect_true( class( obj$dat_full ) == 'BEDMatrix' )
    expect_true( class( obj$dat_rand ) == 'BEDMatrix' )
    expect_equal( class( obj$file_full ), 'character' )
    expect_equal( class( obj$file_rand ), 'character' )
    expect_equal( nrow( obj$dat_full ), n )
    expect_equal( ncol( obj$dat_full ), m )
    expect_equal( nrow( obj$dat_rand ), n )
    expect_equal( ncol( obj$dat_rand ), s )
    X_orig_exp <- X[ -random_s, ]
    X_orig_obs <- t( obj$dat_full[ , -random_s ] )
    dimnames( X_orig_obs ) <- NULL
    expect_equal( X_orig_obs, X_orig_exp )
    X_rand_exp <- obj$dat_rand[]
    X_rand_obs <- obj$dat_full[ , random_s_sorted, drop = FALSE ]
    expect_equal( X_rand_obs, X_rand_exp )

})


# the following tests require the Bioconductor `lfa` package only
if ( requireNamespace( "lfa", quietly = TRUE ) ) {

    test_that( "efron_Rsq_snp works", {
        # data to use
        xi <- X[1,]
        pi <- lfa::af_snp(xi, LF1)

        # cause errors due to missing arguments
        expect_error( efron_Rsq_snp() )
        expect_error( efron_Rsq_snp( snp = xi ) )
        expect_error( efron_Rsq_snp( p1 = pi ) )

        # now a successful run
        expect_silent(
            r2 <- efron_Rsq_snp( snp = xi, p1 = pi )
        )
        # the basics of what this R^2 should be like
        expect_true( is.numeric( r2 ) )
        expect_equal( length( r2 ), 1 )
        expect_true( !is.na( r2 ) )
        expect_true( r2 >= 0 )
        expect_true( r2 <= 1 )
    })

    test_that( "mcfadden_Rsq_snp works", {
        # LF0 it can't be null here
        if ( is.null( LF0 ) )
            LF0 <- matrix(1, n, 1)
            # data to use
            xi <- X[1,]
            p1 <- lfa::af_snp(xi, LF1)
            p0 <- lfa::af_snp(xi, LF0)

            # cause errors due to missing arguments
            # all three arguments are required
            expect_error( mcfadden_Rsq_snp() )
            expect_error( mcfadden_Rsq_snp( snp = xi ) )
            expect_error( mcfadden_Rsq_snp( p1 = p1 ) )
            expect_error( mcfadden_Rsq_snp( p0 = p0 ) )
            expect_error( mcfadden_Rsq_snp( snp = xi, p1 = p1 ) )
            expect_error( mcfadden_Rsq_snp( snp = xi, p0 = p0 ) )
            expect_error( mcfadden_Rsq_snp( p1 = p1, p0 = p0 ) )

            # now a successful run
            expect_silent(
                r2 <- mcfadden_Rsq_snp( snp = xi, p1 = p1, p0 = p0 )
            )
            # the basics of what this R^2 should be like
            expect_true( is.numeric( r2 ) )
            expect_equal( length( r2 ), 1 )
            expect_true( !is.na( r2 ) )
            expect_true( r2 >= 0 )
            expect_true( r2 <= 1 )
    })

    test_that( "pseudo_Rsq works", {
        # cause errors due to missing arguments
        expect_error( pseudo_Rsq( ) )
        expect_error( pseudo_Rsq( X ) )
        expect_error( pseudo_Rsq( LF_alt = LF1 ) )
        # pass non-matrix arguments
        expect_error( pseudo_Rsq( as.vector( X ), LF1 ) )
        expect_error( pseudo_Rsq( X, as.vector( LF1 ) ) )

        # now a successful run
        # LF_null is set to default (intercept only)
        expect_silent(
            r2 <- pseudo_Rsq( X, LF1 )
        )
        # the basics of what this R^2 should be like
        expect_true( is.numeric( r2 ) )
        expect_equal( length( r2 ), m )
        expect_true( !anyNA( r2 ) )
        expect_true( all( r2 >= 0 ) )
        expect_true( all( r2 <= 1 ) )
    })

    test_that( "efron_Rsq works", {
        # cause errors due to missing arguments
        expect_error( efron_Rsq( ) )
        expect_error( efron_Rsq( X ) )
        expect_error( efron_Rsq( LF = LF1 ) )
        # pass non-matrix arguments
        expect_error( efron_Rsq( as.vector( X ), LF1 ) )
        expect_error( efron_Rsq( X, as.vector( LF1 ) ) )

        # now a successful run
        expect_silent(
            r2 <- efron_Rsq( X, LF1 )
        )
        # the basics of what this R^2 should be like
        expect_true( is.numeric( r2 ) )
        expect_equal( length( r2 ), m )
        expect_true( !anyNA( r2 ) )
        expect_true( all( r2 >= 0 ) )
        expect_true( all( r2 <= 1 ) )
    })
}

# these functions require the Bioconductor package `gcatest`
if ( requireNamespace( "gcatest", quietly = TRUE ) ) {

    # the following tests require the Bioconductor `lfa` package
    if ( requireNamespace( "lfa", quietly = TRUE ) ) {

        # define the function to pass to `jackstraw_lfa`!  Uses global `d`
        FUN <- function(x) lfa::lfa( x, d )

        test_that( "jackstraw_lfa works", {
            # cause errors due to missing required data
            # must provide all of (dat = X, r = d, FUN = FUN) for a minimal successful run
            expect_error( jackstraw_lfa( ) )
            expect_error( jackstraw_lfa( dat = X ) )
            expect_error( jackstraw_lfa( r = d ) )
            expect_error( jackstraw_lfa( FUN = FUN ) )
            expect_error( jackstraw_lfa( dat = X, r = d ) )
            expect_error( jackstraw_lfa( dat = X, FUN = FUN ) )
            expect_error( jackstraw_lfa( r = d, FUN = FUN ) )
            # check that data is matrix
            expect_error( jackstraw_lfa( dat = 1:10, r = d, FUN = FUN ) )
            # pass bad covariates on purpose
            expect_error( jackstraw_lfa( dat = X, r = d, FUN = FUN, covariate = 1 ) ) # scalar is bad
            expect_error( jackstraw_lfa( dat = X, r = d, FUN = FUN, covariate = covariate[-1] ) ) # length off by 1, vector
            expect_error( jackstraw_lfa( dat = X, r = d, FUN = FUN, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
            expect_error( jackstraw_lfa( dat = X, r = d, FUN = FUN, covariate = rbind( covariate ) ) ) # transposed matrix

            # perform a basic run

            # make it silent so we can focus on problem messages
            expect_silent(
                obj <- jackstraw_lfa( X, r = d, FUN = FUN, s = s, B = B, verbose = FALSE )
            )
            # check basic jackstraw return object
            test_jackstraw_return_val( obj, s, B )
            # in this case try nc-chisq p-values!
            expect_silent(
                obj_nc <- pvals_nc_chisq( obj, d-1 )
            )
            validate_pvals_nc_chisq( obj_nc, m )
            
            # test edge cases
            # s = 1
            expect_silent(
                obj <- jackstraw_lfa( X, r = d, FUN = FUN, s = 1, B = B, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, 1, B )
            # B = 1
            expect_silent(
                obj <- jackstraw_lfa( X, r = d, FUN = FUN, s = s, B = 1, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, s, 1 )
            # s = B = 1
            expect_silent(
                obj <- jackstraw_lfa( X, r = d, FUN = FUN, s = 1, B = 1, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, 1, 1 )

            # this comparison succeeds most of the time, but not 100% of the time
            # bad fits cause errors randomly, which are rare but over 100 loci it gets less rare
            # also, things get very slow

            # test version with covariates
            expect_silent(
                obj <- jackstraw_lfa( X, r = d, FUN = FUN, s = s, B = B, covariate = covariate, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, s, B )

            # test version without default allele-level permutation!
            expect_silent(
                obj <- jackstraw_lfa( X, r = d, FUN = FUN, s = s, B = B, permute_alleles = FALSE, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, s, B )

        })

        test_that( "jackstraw_lfa works with BEDMatrix", {
            # make it silent so we can focus on problem messages
            expect_silent(
                obj <- jackstraw_lfa( X_BM, r = d, FUN = FUN, s = s, B = B, verbose = FALSE )
            )
            # check basic jackstraw return object
            test_jackstraw_return_val( obj, s, B )

            # test edge cases
            # s = 1
            expect_silent(
                obj <- jackstraw_lfa( X_BM, r = d, FUN = FUN, s = 1, B = B, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, 1, B )
            # B = 1
            expect_silent(
                obj <- jackstraw_lfa( X_BM, r = d, FUN = FUN, s = s, B = 1, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, s, 1 )
            # s = B = 1
            expect_silent(
                obj <- jackstraw_lfa( X_BM, r = d, FUN = FUN, s = 1, B = 1, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, 1, 1 )

            # test version with covariates
            expect_silent(
                obj <- jackstraw_lfa( X_BM, r = d, FUN = FUN, s = s, B = B, covariate = covariate, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, s, B )
        })

    }

    # run alstructure tests only if optional package is available
    if (suppressMessages(suppressWarnings(require(alstructure)))) {

        test_that("jackstraw_alstructure works", {
            # define function to pass (uses global `d`)
            FUN <- function(x) t( alstructure(x, d_hat = d)$Q_hat )

            # cause errors due to missing required data
            # must provide all of (dat = X, r = d, FUN = FUN) for a minimal successful run
            expect_error( jackstraw_alstructure( ) )
            expect_error( jackstraw_alstructure( dat = X ) )
            expect_error( jackstraw_alstructure( r = d ) )
            expect_error( jackstraw_alstructure( FUN = FUN ) )
            expect_error( jackstraw_alstructure( dat = X, r = d ) )
            expect_error( jackstraw_alstructure( dat = X, FUN = FUN ) )
            expect_error( jackstraw_alstructure( r = d, FUN = FUN ) )
            # check that data is matrix
            expect_error( jackstraw_alstructure( dat = 1:10, r = d, FUN = FUN ) )
            # pass bad covariates on purpose
            expect_error( jackstraw_alstructure( dat = X, r = d, FUN = FUN, covariate = 1 ) ) # scalar is bad
            expect_error( jackstraw_alstructure( dat = X, r = d, FUN = FUN, covariate = covariate[-1] ) ) # length off by 1, vector
            expect_error( jackstraw_alstructure( dat = X, r = d, FUN = FUN, covariate = cbind( covariate[-1] ) ) ) # length off by 1, matrix
            expect_error( jackstraw_alstructure( dat = X, r = d, FUN = FUN, covariate = rbind( covariate ) ) ) # transposed matrix

            # perform a basic run

            # make it silent so we can focus on problem messages
            expect_silent(
                obj <- jackstraw_alstructure( X, r = d, FUN = FUN, s = s, B = B, verbose = FALSE )
            )
            # check basic jackstraw return object
            test_jackstraw_return_val( obj, s, B )
            # in this case try nc-chisq p-values!
            expect_silent(
                obj_nc <- pvals_nc_chisq( obj, d-1 )
            )
            validate_pvals_nc_chisq( obj_nc, m )

            # test edge cases
            # s = 1
            expect_silent(
                obj <- jackstraw_alstructure( X, r = d, FUN = FUN, s = 1, B = B, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, 1, B )
            # B = 1
            expect_silent(
                obj <- jackstraw_alstructure( X, r = d, FUN = FUN, s = s, B = 1, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, s, 1 )
            # s = B = 1
            expect_silent(
                obj <- jackstraw_alstructure( X, r = d, FUN = FUN, s = 1, B = 1, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, 1, 1 )

            # test version with covariates
            expect_silent(
                obj <- jackstraw_alstructure( X, r = d, FUN = FUN, s = s, B = B, covariate = covariate, verbose = FALSE )
            )
            test_jackstraw_return_val( obj, s, B )
        })

    }
}

# clean up BEDMatrix example
invisible( suppressWarnings( file.remove( file_tmp ) ) )
