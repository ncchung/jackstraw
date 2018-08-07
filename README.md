# jackstraw: Statistical Inference for Unsupervised Learning

This **R** package performs association tests between the observed data and their systematic patterns of variation. Systematic variation can be modeled by latent variables, that are likely arising from biological processes, experimental conditions, and environmental factors. We are often interested in estimating these patterns using principal component analysis (PCA), factor analysis (FA), K-means clustering, partition around medoids (PAM), and related methods. The jackstraw methods learn over-fitting characteristics inherent in unsupervised learning, where the observed data are used to estimate the systematic patterns and to be tested again.

Using a variety of unsupervised learning techniques, the jackstraw provides a resampling strategy and testing scheme to estimate statistical significance of association between the observed data and their systematic patterns of variation. For example, the cell cycle in microarray data may be estimated by principal components (PCs); then, we can use the jackstraw for PCA to identify genes that are significantly associated with these PCs. On the other hand, cell identities in single cell RNA-seq data are identified by K-means clustering; then, the jackstraw for clustering can evaluate reliability of computationally determined cell identities.

The jackstraw tests enable us to identify the data features (i.e., variables or observations) that are driving systematic variation, in an unsupervised manner. Using **jackstraw_pca**, we can find statistically significant features with regard to the top r principal components. Alternatively, **jackstraw_kmeans** can identify the data features that are statistically significant members of the data-dependent clusters. Furthermore, this package includes more general algorithms such as **jackstraw_subspace** for the dimension reduction techniques and **jackstraw_cluster** for the clustering algorithms.

*Chung, N.C.* and *Storey, J.D.* (2015) Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics, 31(4): 545-554
http://bioinformatics.oxfordjournals.org/content/31/4/545

# Installation

To use a stable version from CRAN:
```R
install.packages("jackstraw")
```

To use a development version from GitHub:
```R
install.packages("devtools")
library("devtools")
install_github("ncchung/jackstraw")
```
