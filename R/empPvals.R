# copy of `qvalue::empPvals` modified to handle NAs in input correctly, and with hardcoded `pool=TRUE` for simplicity
# secondary reason to copy is to limit dependency on Bioconductor packages (including qvalue), and this function in particular is widely used!
empPvals <- function (stat, stat0) {
    m <- length(stat)
    m0 <- length(stat0)
    if (is.matrix(stat0))
        stat0 <- as.vector(stat0)
    v <- c(rep(TRUE, m), rep(FALSE, m0))
    v <- v[order(c(stat, stat0), decreasing = TRUE)]
    u <- 1:length(v)
    w <- 1:m
    p <- (u[v == TRUE] - w)/m0
    p <- p[rank(-stat)]
    p <- pmax(p, 1/m0)
    # fix NA situation
    if ( anyNA( stat ) )
        p[ is.na( stat ) ] <- NA
    return(p)
}
