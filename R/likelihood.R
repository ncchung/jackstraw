#' Mcfadden's Pseudo R-squared
#'
#' @param X a data matrix.
#' @param LF_alt Observed logistic factors.
#' @param LF_null Null logistic factors.
#' @author Wei Hao
#'
#' @keywords internal
pseudo_Rsq <- function(X, LF_alt, LF_null = NULL) {
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( LF_alt ) )
        stop( '`LF_alt` is required!' )
    if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )
    if ( !is.matrix( LF_alt ) )
        stop( '`LF_alt` must be a matrix!' )
    
    if (is.null(LF_null))
        LF_null <- matrix(1, ncol(X), 1)
    
    m <- nrow(X)
    
    F_alt <- lfa::af(X, LF_alt)
    F_null <- lfa::af(X, LF_null)
    
    sapply(1:m, function(i) {
        mcfadden_Rsq_snp(
            X[i, ],
            F_alt[i, ],
            F_null[i, ]
        )
    })
}

#' @keywords internal
mcfadden_Rsq_snp <- function(snp, p1, p0) {
    if ( missing( snp ) )
        stop( '`snp` is required!' )
    if ( missing( p1 ) )
        stop( '`p1` is required!' )
    if ( missing( p0 ) )
        stop( '`p0` is required!' )
    
    # check for p's = 0 or 1
    IND <- (p0 != 0) & (p0 != 1) & (p1 != 0) & (p1 != 1)
    p1 <- p1[IND]
    p0 <- p0[IND]
    snp <- snp[IND]

    llalt <- sum(snp * log(p1) + (2 - snp) * log(1 - p1))
    llnull <- sum(snp * log(p0) + (2 - snp) * log(1 - p0))
    
    1 - (llalt/llnull)
}

#' Efron's Pseudo R-squared
#'
#' @param X a data matrix.
#' @param LF Observed logistic factors.
#' @author Wei Hao
#'
#' @keywords internal
efron_Rsq <- function(X, LF) {
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( LF ) )
        stop( '`LF` is required!' )
    if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )
    if ( !is.matrix( LF ) )
        stop( '`LF` must be a matrix!' )
    
    m <- nrow(X)
    
    F <- lfa::af(X, LF)

    sapply(1:m, function(i) {
        efron_Rsq_snp(X[i, ], F[i, ])
    })
}

#' @keywords internal
efron_Rsq_snp <- function(snp, p1) {
    if ( missing( snp ) )
        stop( '`snp` is required!' )
    if ( missing( p1 ) )
        stop( '`p1` is required!' )
    
    IND <- (p1 != 0) & (p1 != 1)
    p1 <- p1[IND]
    snp <- snp[IND]

    y <- as.numeric(c(snp > 0, snp == 2))
    p <- c(p1, p1)
    ybar <- mean(y)

    1 - sum((y - p)^2)/sum((y - ybar)^2)
}
