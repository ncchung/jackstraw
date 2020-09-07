#' Pi0 estimation using a non-analytical bootstrap algorithm
#'
#' This is a function to estimate pi0 using the bootstrap, directly dervied from the qvalue version 1.43 (Bioconductor 3.0).
#' Try \code{qvalue::pi0est} or other pi0 estimation methods. This function is preserved because in the newer qvalue package (version > 2.0), the bootstrap method is replaced by a closed form solution.
#' In practice, the closed form solution may fail due to a very small pi0 (e.g., generally small p-values) or others.
#'
#' @param pvalue a vector of p-values.
#' @param lambda the value of the tuning parameter to estimate pi0 (optional)
#' @param robust an indicator of whether it is desired to make the estimate more robust for small p-values and a direct finite sample estimate of pFDR (optional)
#' @param verbose print a progress
#' @param ... optional arguments to control a local FDR estimation.
#'
#' @return \code{pi0est_bootstrap} returns a pi0 estimate
#'
#' @export pi0est_bootstrap
#'
#' @importFrom qvalue lfdr
#' @author John D. Storey, Alan Dabney
pi0est_bootstrap <- function(p=NULL, lambda=seq(0,0.90,0.05), robust=FALSE) {
    if(min(p)<0 || max(p)>1) {
        print("ERROR: p-values not in valid range.")
      return(0)
    }
    if(length(lambda)>1 && length(lambda)<4) {
        print("ERROR: If length of lambda greater than 1, you need at least 4 values.")
      return(0)
    }
    if(length(lambda)>1 && (min(lambda) < 0 || max(lambda) >= 1)) { ## change by Alan:  check for valid range for lambda
        print("ERROR: Lambda must be within [0, 1).")
      return(0)
    }
    m <- length(p)

    if(length(lambda)==1) {
        if(lambda<0 || lambda>=1) { ## change by Alan:  check for valid range for lambda
            print("ERROR: Lambda must be within [0, 1).")
          return(0)
        }

        message(paste("Using a single lambda value of ",lambda,"."))
        pi0 <- mean(p >= lambda)/(1-lambda)
        pi0 <- min(pi0,1)
    } else {
        message(paste("Using",length(lambda), "lambda values."))

        pi0 <- rep(0,length(lambda))
        for(i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])
        }

        minpi0 <- min(pi0)
        mse <- rep(0,length(lambda))
        pi0.boot <- rep(0,length(lambda))
        for(i in 1:100) {
            p.boot <- sample(p,size=m,replace=TRUE)
            for(i in 1:length(lambda)) {
                pi0.boot[i] <- mean(p.boot>lambda[i])/(1-lambda[i])
            }
            mse <- mse + (pi0.boot-minpi0)^2
        }
        pi0 <- min(pi0[mse==min(mse)])
        pi0 <- min(pi0,1)
    }
    if(pi0 <= 0) {
        print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
      return(0)
    }
    return(pi0)
}
