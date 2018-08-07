# jackstraw: Statistical Inference for Unsupervised Learning

This **R** package performs association tests between the observed data and their latent variables. Latent variables may be estimated by principal component analysis, factor analysis, K-means clustering, and related methods.

The jackstraw package provides a resampling strategy and testing scheme to estimate statistical significance of association between the observed data and their latent variables. Depending on the data type and the analysis aim, the latent variables may be estimated by principal component analysis, K-means clustering, and related algorithms. The jackstraw methods learn over-fitting characteristics inherent in this circular analysis, where the observed data are used to estimate the latent variables and to again test against the estimated latent variables.

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
