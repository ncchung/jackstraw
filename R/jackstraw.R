#' @details
#' The jackstraw package provides a resampling strategy and testing scheme to estimate statistical significance
#' of association between the observed data and their latent variables. Depending on the data type and the analysis aim,
#' the latent variables may be estimated by principal component analysis (PCA), K-means clustering, and related algorithms.
#' The jackstraw methods learn over-fitting characteristics inherent in this circular analysis, where the observed data are used
#' to estimate the latent variables and used again to test against those estimated latent variables.
#'
#' The jackstraw tests enable us to identify the data features (i.e., variables or observations) that are
#' driving systematic variation, in a completely unsupervised manner. Using \link{jackstraw_pca}, we can
#' find statistically significant features with regard to the top \code{r} principal components.
#' Alternatively, \link{jackstraw_kmeans} can identify the data features that are statistically significant
#' members of the data-dependent clusters. Furthermore, this package includes more general algorithms such as
#' \link{jackstraw_subspace} for the dimension reduction techniques and \link{jackstraw_cluster} for the clustering algorithms.
#'
#' Overall, it computes \code{m} p-values of association between the \code{m} data features and their corresponding latent variables.
#' From \code{m} p-values, \link{pip} computes posterior inclusion probabilities, that are useful for feature selection and visualization.
#'
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2015) Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics, 31(4): 545-554 \url{https://academic.oup.com/bioinformatics/article/31/4/545/2748186}
#' @references Chung (2020) Statistical significance of cluster membership for unsupervised evaluation of cell identities. Bioinformatics, 36(10): 3107–3114 \url{https://academic.oup.com/bioinformatics/article/36/10/3107/5788523}
#'
#' @seealso \link{jackstraw_pca} \link{jackstraw_subspace} \link{jackstraw_kmeans} \link{jackstraw_cluster}
#' @docType package
#' @name jackstraw
"_PACKAGE"
