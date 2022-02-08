# wrapper around `qvalue::empPvals` that handles NAs in input correctly
empPvals <- function( obs, null ) {
    p.value <- qvalue::empPvals( obs, null )
    # fix NA situation (hack until ideally it is fixed in `qvalue::empPvals`)
    if ( anyNA( obs ) )
        p.value[ is.na( obs ) ] <- NA
    return( p.value )
}
