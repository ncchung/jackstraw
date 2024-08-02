#' Non-Parametric Jackstraw for Logistic Factor Analysis
#'
#' Test association between the observed variables and their latent variables captured by logistic factors (LFs).
#'
#' This function uses logistic factor analysis (LFA) from Hao et al. (2016).
#' Particularly, the deviance in logistic regression (the full model with \code{r} LFs vs. the intercept-only model) is used to assess significance.
#'
#' The random outputs of the regular matrix versus the \code{BEDMatrix} versions are equal in distribution.
#' However, fixing a seed and providing the same data to both versions does not result in the same exact outputs.
#' This is because the \code{BEDMatrix} version permutes loci in a different order by necessity.
#'
#' @param dat either a genotype matrix with \code{m} rows as variables and \code{n} columns as observations, or a \code{BEDMatrix} object (see package \code{BEDMatrix}, these objects are transposed compared to the above but this works fine as-is, see example, no need to modify a \code{BEDMatrix} input).
#' A \code{BEDMatrix} input triggers a low-memory mode where permuted data is also written and processed from disk, whereas a regular matrix input stores permutations in memory.
#' The tradeoff is \code{BEDMatrix} version typically runs considerably slower, but enables analysis of very large data that is otherwise impossible.
#' @param r a number of significant LFs.
#' @param FUN a function to use for LFA (by default, it uses the \code{lfa} package)
#' @param r1 a numeric vector of LFs of interest (implying you are not interested in all \code{r} LFs).
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations. There will be a total of \code{s*B} null statistics.
#' @param covariate a data matrix of covariates with corresponding \code{n} observations (do not include an intercept term).
#' @param permute_alleles If TRUE (default), alleles (rather than genotypes) are permuted, which results in a more Binomial synthetic null when data is highly structured.
#' Changing to FALSE is not recommended, except for research purposes to confirm that it performs worse than the default.
#' @param verbose a logical specifying to print the computational progress.
#'
#' @return \code{jackstraw_lfa} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their LFs}
#' \item{obs.stat}{\code{m} observed deviances}
#' \item{null.stat}{\code{s*B} null deviances}
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2015) Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics, 31(4): 545-554 \url{https://academic.oup.com/bioinformatics/article/31/4/545/2748186}
#' @seealso  \link{jackstraw_pca} \link{jackstraw} \link{jackstraw_subspace}
#'
#' @examples
#' \dontrun{
#' ## simulate genotype data from a logistic factor model: drawing rbinom from logit(BL)
#' m <- 5000; n <- 100; pi0 <- .9
#' m0 <- round(m*pi0)
#' m1 <- m - round(m*pi0)
#' B <- matrix(0, nrow=m, ncol=1)
#' B[1:m1,] <- matrix(runif(m1*n, min=-.5, max=.5), nrow=m1, ncol=n)
#' L <- matrix(rnorm(n), nrow=1, ncol=n)
#' BL <- B %*% L
#' prob <- exp(BL)/(1+exp(BL))
#'
#' dat <- matrix(rbinom(m*n, 2, as.numeric(prob)), m, n)
#'
#' ## apply the jackstraw_lfa
#' out <- jackstraw_lfa(dat, r = 2)
#'
#' # if you had very large genotype data in plink BED/BIM/FAM files,
#' # use BEDMatrix and save memory by reading from disk (at the expense of speed)
#' library(BEDMatrix)
#' dat_BM <- BEDMatrix( 'filepath' ) # assumes filepath.bed, .bim and .fam exist
#' # run jackstraw!
#' out <- jackstraw_lfa(dat_BM, r = 2)
#' }
#'
#' @export
jackstraw_lfa <- function(
                          dat,
                          r,
                          FUN = function(x) lfa::lfa(x, r),
                          r1 = NULL,
                          s = NULL,
                          B = NULL,
                          covariate = NULL,
                          permute_alleles = TRUE,
                          verbose = TRUE
                          ) {
    # check mandatory data
    if ( missing( dat ) )
        stop( '`dat` is required!' )
    if ( missing( r ) )
        stop( '`r` is required!' )
    if ( !is.matrix( dat ) )
        stop( '`dat` must be a matrix!' )

    if ( !is.null( FUN ) ) {
        FUN <- match.fun( FUN )
    } else {
        stop("Please provide a function to estimate latent variables.")
    }

    # make sure `dat` is either a regular matrix or BEDMatrix
    # also get correct dimensions
    is_BEDMatrix <- FALSE
    if ( 'BEDMatrix' %in% class( dat ) ) {
        is_BEDMatrix <- TRUE
        # dimensions are transposed for this object
        n <- nrow( dat )
        m <- ncol( dat )
    } else if ( is.matrix( dat ) ) {
        # regular case
        m <- nrow( dat )
        n <- ncol( dat )
    } else
        stop( '`dat` must be a matrix or BEDMatrix object!  Observed class: ', class( dat ) )

    # more validations of mandatory parameters
    if ( !(r > 0 && r < n) )
        stop( "`r` is not in valid range between `1` and `n-1` (`n` is number of individuals)." )

    # if there are covariates, the dimensions must agree
    # covariate can be either a vector or a matrix, test both cases
    if ( !is.null( covariate ) ) {
        if ( is.matrix( covariate ) ) {
            if ( nrow( covariate ) != n )
                stop( 'Matrix `covariate` must have `n` rows, has: ', nrow( covariate ), ', expected: ', n )
        } else {
            if ( length( covariate ) != n )
                stop( 'Vector `covariate` must have `n` elements, has: ', length( covariate ), ', expected: ', n )
        }
    }

    if (is.null(s)) {
        s <- round(m/10)
        if (verbose)
            message( "Number of null variables (s) to be permuted is not specified: s=round(0.10*m)=", s, "." )
    }
    if (is.null(B)) {
        B <- round(m * 10/s)
        if (verbose)
            message( "Number of resampling iterations (B) is not specified: B=round(m*10/s)=", B, "." )
    }

    if ( is.null( r1 ) )
        r1 <- 1:r
    if ( all( 1:r %in% r1 ) ) {
        # no adjustment LVs
        r0 <- NULL
        LFr0 <- NULL
    } else {
        r0 <- ( 1 : r )[ -r1 ]
    }

    ## LFr[,r] is an intercept term
    # NOTE: `lfa::lfa` supports BEDMatrix as input `dat`, returns regular matrix `LF`
    LFr <- FUN( dat )
    # check dimensions for this matrix
    if (ncol(LFr) != r)
        stop("The number of latent variables ", r, " is not equal to the number of column(s) provided by `FUN`")
    if (nrow(LFr) != n)
        stop("The number of individuals ", n, " is not equal to the number of rows provided by `FUN`")
    # suubset further if `r1` is non-trivial
    LFr1 <- LFr[ , r1, drop = FALSE ]
    # create null LFs if `r1` is non-trivial (otherwise LFr0 is NULL)
    if (!is.null(r0))
        LFr0 <- LFr[ , r0, drop = FALSE ]

    # NOTE: this `gcatest` function supports BEDMatrix as input `X`/`dat`!
    # both LF models get the same covariates appended (as they should be)
    # the resulting matrices are never `NULL` (that is not handled by this `gcatest` function).
    # there is an awkward ambiguity here, that LFr0 is assumed to not contain the intercept (so it is added), but regression will fail if that assumption is wrong!
    obs <- gcatest::delta_deviance_lf(
                        X = dat,
                        LF0 = cbind(LFr0, matrix(1, n, 1), covariate),
                        LF1 = cbind(LFr, covariate)
                    )

    # Estimate null association
    # statistics
    null <- matrix(0, nrow = s, ncol = B)
    LFr0.js <- NULL

    if ( verbose )
        cat( paste0( "\nComputating null statistics (", B, " total iterations): " ) )
    for (i in 1:B) {
        # select the random subset of rows/loci to edit
        random_s <- sample.int( m, s )
        if ( is_BEDMatrix ) {
            objBM <- jackstraw_BEDMatrix( dat, random_s, permute_alleles = permute_alleles )
            # both of these are BEDMatrix objects
            jackstraw.dat <- objBM$dat_full
            # NOTE: the order of loci in this s.nulls is different than in the non-BEDMatrix version below (here `sort( random_s )` rather than `random_s` below), which doesn't affect the output p-values, but since the random steps are in a different order, even fixing a seed results in different random steps in BEDMatrix vs non-BEDMatrix versions.
            s.nulls <- objBM$dat_rand
        } else {
            # extract that data and permute it, returning a matrix
            s.nulls <- dat[ random_s, , drop = FALSE ]
            s.nulls <- t( apply( s.nulls, 1L, if ( permute_alleles ) permute_alleles_from_geno else sample ) )
            # make a copy of this whole data containing the random subset
            jackstraw.dat <- dat
            jackstraw.dat[random_s, ] <- s.nulls
        }

        # get LFs for this random data
        LFr.js <- FUN(jackstraw.dat)
        # dimensions are not rechecked here, it's just assumed that they are correct
        LFr1.js <- LFr.js[, r1, drop = FALSE]
        if (!is.null(r0))
            LFr0.js <- LFr.js[, r0, drop = FALSE]

        null[, i] <- gcatest::delta_deviance_lf(
                            X = s.nulls,
                            LF0 = cbind(LFr0.js, matrix(1, n, 1), covariate),
                            LF1 = cbind(LFr.js, covariate)
                        )

        # now we're done with these temporary files
        # NOTE: on Windows there's a peculiar issue, that these temporary files cannot be removed because BEDMatrix left them "open", silence those warnings!
        if ( is_BEDMatrix ) {
            # try to delete, ignore warnings if it failed
            invisible( suppressWarnings( file.remove( objBM$file_full ) ) )
            invisible( suppressWarnings( file.remove( objBM$file_rand ) ) )
        }

        if ( verbose )
            cat(paste(i, " "))
    }

    p.value <- empPvals( obs, null )
    
    return(
        list(
            call = match.call(),
            p.value = p.value,
            obs.stat = obs,
            null.stat = null
        )
    )
}


# internal function for creating jackstraw random samples on disk, via BEDMatrix
jackstraw_BEDMatrix <- function(
                                dat,
                                random_s,
                                permute_alleles = FALSE,
                                m_chunk = 1000
                                ) {
    if ( missing( dat ) )
        stop( '`dat` is required!' )
    if ( missing( random_s ) )
        stop( '`random_s` is required!' )

    # check class, for now makes sense to only support BEDMatrix
    if ( !( 'BEDMatrix' %in% class( dat ) ) )
        stop( '`dat` must be class "BEDMatrix"!' )

    # dimensions are transposed for this object
    n_ind <- nrow( dat )
    m_loci <- ncol( dat )
    s <- length( random_s ) # needed too

    # the order of `random_s` doesn't really matter outside, but here we can speed things up if indexes are ordered
    # ordering also helps make sense of some internal validations that are not possible otherwise
    # this list will be shortened as we advance in the loop, for efficiency in big `s` cases
    random_s <- sort( random_s )

    # create two temporary files to write modified genotype data to
    # both paths have BED extension only (will save additional time by skipping BIM and FAM, which are completely useless in this setting)
    file_full <- tempfile( 'jackstraw_BEDMatrix_full', fileext = '.bed' )
    file_rand <- tempfile( 'jackstraw_BEDMatrix_rand', fileext = '.bed' )

    # start writing in chunks for balance of memory usage and speed
    i_chunk <- 1
    while ( i_chunk <= m_loci ) {
        # calculate end of current chunk, careful not to exceed total number of loci
        i_chunk_max <- i_chunk + m_chunk - 1
        if ( i_chunk_max > m_loci )
            i_chunk_max <- m_loci
        # range of data to process
        indexes <- i_chunk : i_chunk_max

        # read chunk out of BEDMatrix object, orient immediately as needed for `genio`
        dat_chunk <- t( dat[ , indexes, drop = FALSE ] )

        # shuffle loci as needed (this is the jackstraw-specific bit!)
        # this is the list of indexes that fall in the current chunk
        #random_s_chunk <- random_s[ random_s %in% indexes ] # explicit matching
        random_s_chunk <- random_s[ random_s <= i_chunk_max ] # same thing, relies on continuity of indexes and removal of random_s cases from previous chunks
        # nothing to edit if no indexes were in this range
        if ( length( random_s_chunk ) > 0 ) {
            # map to indexes in the current chunk
            # this simple subtraction gives the correct indexes
            random_s_chunk <- random_s_chunk - i_chunk + 1L

            # extract that data
            dat_rand <- dat_chunk[ random_s_chunk, , drop = FALSE ]
            # now permute it
            dat_rand <- t( apply( dat_rand, 1L, if ( permute_alleles ) permute_alleles_from_geno else sample ) )
            # overwrite those cases immediately
            dat_chunk[ random_s_chunk, ] <- dat_rand

            # write this random data to plink BED
            genio::write_bed(
                       file_rand,
                       dat_rand,
                       verbose = FALSE,
                       append = TRUE
                   )

            # in this case we have to clean up `random_s` a bit
            # remove cases that fell in the chunk that was just processed
            random_s <- random_s[ random_s > i_chunk_max ]
        }

        # write full data to plink BED
        genio::write_bed(
                   file_full,
                   dat_chunk,
                   verbose = FALSE,
                   append = TRUE
               )

        # increment for next round
        i_chunk <- i_chunk_max + 1
    }

    # now that file writing is done, read both back with BEDMatrix, to return
    # provide dimensions for speed and also since BIM and FAM files are absent!
    dat_full <- BEDMatrix::BEDMatrix( file_full, n = n_ind, p = m_loci )
    dat_rand <- BEDMatrix::BEDMatrix( file_rand, n = n_ind, p = s )
    
    # done, return necessary data!
    return(
        list(
            dat_full = dat_full,
            # NOTE: loci in `dat_rand` are in order of appearance, matches `sort( random_s )` compared to `dat_full`!
            # but in non-BEDMatrix jackstraw, order is `random_s`, so fixing a seed won't lead to the same random results!
            dat_rand = dat_rand,
            file_full = file_full,
            file_rand = file_rand
        )
    )
}
