# jackstraw: Statistical Inference for Unsupervised Learning

This **R** package performs association tests between the observed data and their systematic patterns of variation. Systematic variation can be modeled by latent variables, that can arise from biological processes, experimental conditions, environmental factors, and others. We often estimate these patterns using principal component analysis (PCA), factor analysis (FA), logistic factor analysis (LFA), K-means clustering, partition around medoids (PAM), and related methods. The jackstraw methods learn over-fitting characteristics inherent in unsupervised learning, where the observed data are used to estimate the systematic patterns and to be tested again (see [circular analysis](https://en.wikipedia.org/wiki/Circular_analysis)).

Using a variety of unsupervised learning techniques, the jackstraw provides a resampling strategy and testing scheme to estimate statistical significance of association between the observed data and their systematic patterns of variation. For example, the cell cycle in microarray data may be estimated by principal components (PCs). Then, we can use the jackstraw for PCA to identify genes that are significantly associated with these PCs. On the other hand, cell identities in single cell RNA-seq (scRNA-seq) data are often determined by K-means clustering or other unsupervised clustering algorithms. Then, the jackstraw for clustering can identify single cells that are significant members of a given cluster.

## Use cases

Using `jackstraw_pca`, we can find statistically significant variables with regard to the top `r` principal components (PCs). Alternatively, we could test association with respect to a subset of `r` PCs, which are called PCs of interest `r1`. The package also supports truncated PCA, using augmented implicitly restarted Lanczos bidiagonalization algorithm (IRLBA; `jackstraw_irlba`) or randomized Singular Value Decomposition (RSVD; `jackstraw_rpca`). 

[Logistic factor analysis (LFA)](https://academic.oup.com/bioinformatics/article/32/5/713/1744055) and [ALStructure](https://academic.oup.com/genetics/article/212/4/1009/5931257?login=false) estimate population structure from genetic data ([single-nucleotide polymorphisms](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism); SNPs). `jackstraw_lfa` and `jackstraw_alstructure` provides corresponding association tests between SNPs and population structure, as estimated by the aforementioned methods. Generally, one could directly specify an estimation method for latent variables in `jackstraw_subspace`.

Instead of continuous latent variables that are estimated by PCA, LFA, or others, one may be interested in estimating discrete clusters from a high dimensional data. For K-means clustering, `jackstraw_kmeans` evaluates whether data points are significant members of a given cluster, by testing association between observed data and cluster centers. This can help select data points that are reliable members of clusters and further improve the cluster membership. 

Related algorithms, such as [Partitioning Around Medoids (PAM) or k-medoids](https://en.wikipedia.org/wiki/K-medoids) and [Mini Batch K-means](https://dl.acm.org/doi/10.1145/1772690.1772862) algorithms, are supported by `jackstraw_pam` and `jackstraw_MiniBatchKmeans`, respectively. Generally, `jackstraw_cluster` can be used for other clustering algorithms.

There are few additional functions to support statistical inference for unsupervised learning, such as finding a number of PCs or clusters. Based on p-values, we could estimate posterior inclusion probabilities (PIPs) using `pip`.

# References

*Chung, N.C.* (2020) Statistical significance of cluster membership for unsupervised evaluation of cell identities. Bioinformatics, 36(10): 3107–3114
https://academic.oup.com/bioinformatics/article/36/10/3107/5788523

*Chung, N.C.* and *Storey, J.D.* (2015) Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics, 31(4): 545-554
https://academic.oup.com/bioinformatics/article/31/4/545/2748186

# Short Tutorials

[Association Test with Principal Components with a Gentle Introduction to Latent Variable Models](https://cbml.science/post/association-test-with-principal-components/)

[Statistical Test of Cluster Memberships with a Toy Data Set (`mtcars`)](https://cbml.science/post/test-of-cluster-memberships/)

[Unsupervised Evaluation of Cell Identities in Single Cell Genomics using the 10X Genomics Data](https://cbml.science/post/unsupervised-evaluation-of-cell-identities/)

# Installation

## Bioconductor dependencies

Bioconductor dependencies may fail to automatically install, e.g., [`lfa`](https://bioconductor.org/packages/release/bioc/html/lfa.html), [`gcatest`](https://bioconductor.org/packages/release/bioc/html/gcatest.html), [`qvalue`](https://bioconductor.org/packages/release/bioc/html/qvalue.html). This would result in a [warning](https://github.com/ncchung/jackstraw/issues/2).

To solve this problem, please install Bioconductor dependencies manually first:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c('qvalue', 'lfa', 'gcatest'))
```

## Development Version on GitHub

This package is in active development. 
Install jackstraw from GitHub:
```R
install.packages("devtools")
library("devtools")
install_github("ncchung/jackstraw")
```

To use `jackstraw_alstructure`, install the optional `alstructure` package from GitHub: 
```R
library(devtools)
install_github("StoreyLab/alstructure")
```

## Stable Version on CRAN

The stable version **jackstraw v1.3.14** is on CRAN. To install a stable version from CRAN:
```R
install.packages("jackstraw")
```

# Implementations and Extensions

Here are some implementations of the jackstraw in different contexts and application domains.

#### Implementation of the jackstraw in Python is available:

[jackstraw (Python) by Iain Carmichael](https://github.com/idc9/jackstraw)

#### Extension of [Jackstraw Inference for AJIVE Data Integration](https://arxiv.org/abs/2109.12272):

[Jackstraw significance testing for JIVE in Python](https://github.com/thomaskeefe/jive_jackstraw)

#### The jackstraw used in [Seurat](https://satijalab.org/seurat/), R toolkit for single cell genomics:

[Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)

[Determine statistical significance of PCA scores](https://satijalab.org/seurat/reference/jackstraw)

[Seurat Wizard (GUI Web App)](http://nasqar2.abudhabi.nyu.edu/SeuratV3Wizard/)

