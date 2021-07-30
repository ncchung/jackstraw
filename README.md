# jackstraw: Statistical Inference for Unsupervised Learning

This **R** package performs association tests between the observed data and their systematic patterns of variation. Systematic variation can be modeled by latent variables, that are likely arising from biological processes, experimental conditions, and environmental factors. We are often interested in estimating these patterns using principal component analysis (PCA), factor analysis (FA), K-means clustering, partition around medoids (PAM), and related methods. The jackstraw methods learn over-fitting characteristics inherent in unsupervised learning, where the observed data are used to estimate the systematic patterns and to be tested again.

Using a variety of unsupervised learning techniques, the jackstraw provides a resampling strategy and testing scheme to estimate statistical significance of association between the observed data and their systematic patterns of variation. For example, the cell cycle in microarray data may be estimated by principal components (PCs); then, we can use the jackstraw for PCA to identify genes that are significantly associated with these PCs. On the other hand, cell identities in single cell RNA-seq data are identified by K-means clustering; then, the jackstraw for clustering can evaluate reliability of computationally determined cell identities.

The jackstraw tests enable us to identify the variables (or observations) that are driving systematic variation, in an unsupervised manner. Using **jackstraw_pca**, we can find statistically significant variables with regard to the top r principal components. Alternatively, **jackstraw_kmeans** can identify the variables that are statistically significant members of clusters. There are many functions to support statistical inference for unsupervised learning, such as finding a number of PCs or clusters and estimating posterior probabilities from jackstraw p-values. Furthermore, this package includes more general and experimental algorithms such as **jackstraw_subspace** for the dimension reduction techniques and **jackstraw_cluster** for the clustering algorithms.

*Chung, N.C.* (2020) Statistical significance of cluster membership for unsupervised evaluation of cell identities. Bioinformatics, 36(10): 3107–3114
https://academic.oup.com/bioinformatics/article/36/10/3107/5788523

*Chung, N.C.* and *Storey, J.D.* (2015) Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics, 31(4): 545-554
https://academic.oup.com/bioinformatics/article/31/4/545/2748186

# Stable Version on CRAN

To use a stable version from CRAN:
```R
install.packages("jackstraw")
```

#### Troubleshooting

Bioconductor dependencies may fail to automatically install, namely:

- [lfa](https://bioconductor.org/packages/release/bioc/html/lfa.html)
- [gcatest](https://bioconductor.org/packages/release/bioc/html/gcatest.html)
- [qvalue](https://bioconductor.org/packages/release/bioc/html/qvalue.html)

This would result in a [warning](https://github.com/ncchung/jackstraw/issues/2).:
```R
Error: package or namespace load failed for ‘jackstraw’ in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]):
 there is no package called ‘lfa’
```

To solve this problem, please install these two packages manually using the following command:
```R
# install qvalue from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite('qvalue')
```

# Development Version on GitHub

This package is in active development. 

To install the jackstraw from GitHub:
```R
install.packages("devtools")
library("devtools")
install_github("ncchung/jackstraw")
```

The current GitHub version of `jackstraw` depends on updates for `lfa`, `gcatest`, and `genio` present only on these GitHub repositories:
```R
library(devtools)
install_github("StoreyLab/lfa")
install_github("alexviiia/gcatest")
install_github("OchoaLab/genio")
```
Eventually, the Bioconductor versions of `lfa` and `gcatest` and CRAN version of `genio` will have these updates; sorry for the temporary inconvenience.