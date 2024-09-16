# jackstraw: Statistical Inference for Unsupervised Learning

This **R** package performs association tests between the observed data and their systematic patterns of variation. Systematic variation can be modeled by latent variables, that can arise from biological processes, experimental conditions, environmental factors, and others. We often estimate these patterns using principal component analysis (PCA), factor analysis (FA), logistic factor analysis (LFA), K-means clustering, partition around medoids (PAM), and related methods. The jackstraw methods learn over-fitting characteristics inherent in unsupervised learning, where the observed data are used to estimate the systematic patterns and to be tested again (see [circular analysis](https://en.wikipedia.org/wiki/Circular_analysis)).

Using a variety of unsupervised learning techniques, the jackstraw provides a resampling strategy and testing scheme to estimate statistical significance of association between the observed data and their systematic patterns of variation. For example, the cell cycle in microarray data may be estimated by principal components (PCs). Then, we can use the jackstraw for PCA to identify genes that are significantly associated with these PCs. On the other hand, cell identities in single cell RNA-seq (scRNA-seq) data are often determined by K-means clustering or other unsupervised clustering algorithms. Then, the jackstraw for clustering can identify single cells that are significant members of a given cluster.

## Use cases

Using `jackstraw_pca`, we can find statistically significant variables with regard to the top `r` principal components (PCs). If we only specify `r`, we conduct association tests with all `r` PCs simultaneously. Alternatively, we could test association with respect to a subset of `r` PCs, using an optional argument `r1`. By specifying `r` (a total number of significant PCs) and `r1` (a numeric vector of target PCs), `jackstraw_pca` helps find statistically significant variables with respect to `r1` PCs, while accounting for the fact that there are `r` significant PCs. The package also supports truncated PCA, using augmented implicitly restarted Lanczos bidiagonalization algorithm (IRLBA; `jackstraw_irlba`) or randomized Singular Value Decomposition (RSVD; `jackstraw_rpca`). 

[Logistic factor analysis (LFA)](https://doi.org/10.1093/bioinformatics/btv641) estimates population structure from genetic data ([single-nucleotide polymorphisms](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism); SNPs). `jackstraw_lfa` provides corresponding association tests between SNPs and population structure, as estimated by LFA. Due to the requirements of a CRAN package, please manually install [`lfa`](https://bioconductor.org/packages/release/bioc/html/lfa.html) from Bioconductor. See the R help on `lfa`. In general, one could directly specify an estimation method for latent variables in `jackstraw_subspace`.

Instead of continuous latent variables that are estimated by PCA, LFA, or others, one may be interested in estimating discrete clusters from a high dimensional data. For K-means clustering, `jackstraw_kmeans` evaluates whether data points are significant members of a given cluster, by testing association between observed data and cluster centers. This can help select data points that are reliable members of clusters and further improve the cluster membership. Note that in order to use the jackstraw for clustering, it's necessary to first apply the clustering algorithm to the data and provide the resulting object (e.g., `kmeans.dat`).

Related algorithms, such as [Partitioning Around Medoids (PAM) or k-medoids](https://doi.org/10.1002/9780470316801.ch2) and [Mini Batch K-means](https://doi.org/10.1145/1772690.1772862) algorithms, are supported by `jackstraw_pam` and `jackstraw_MiniBatchKmeans`, respectively. Generally, `jackstraw_cluster` can be used for other clustering algorithms.

There are few additional functions to support statistical inference for unsupervised learning, such as finding a number of PCs or clusters. Based on p-values, we could estimate posterior inclusion probabilities (PIPs) using `pip`.

# References

*Chung, N.C.* (2020) Statistical significance of cluster membership for unsupervised evaluation of cell identities. Bioinformatics, 36(10): 3107–3114
https://doi.org/10.1093/bioinformatics/btaa087

*Chung, N.C.* and *Storey, J.D.* (2015) Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics, 31(4): 545-554
https://doi.org/10.1093/bioinformatics/btu674

# Short Tutorials

[Association Test with Principal Components with a Gentle Introduction to Latent Variable Models](https://cbml.science/post/association-test-with-principal-components/)

[Statistical Test of Cluster Memberships with a Toy Data Set (`mtcars`)](https://cbml.science/post/test-of-cluster-memberships/)

[Unsupervised Evaluation of Cell Identities in Single Cell Genomics using the 10X Genomics Data](https://cbml.science/post/unsupervised-evaluation-of-cell-identities/)

# Installation

## Bioconductor dependencies

Install Bioconductor dependencies,  [`lfa`](https://bioconductor.org/packages/release/bioc/html/lfa.html), [`gcatest`](https://bioconductor.org/packages/release/bioc/html/gcatest.html), [`qvalue`](https://bioconductor.org/packages/release/bioc/html/qvalue.html), manually first:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c('qvalue', 'lfa', 'gcatest'))
```

The following `jackstraw` functions requires Bioconductor packages:

  - `jackstraw_lfa`, `pseudo_Rsq`, `efron_Rsq` requires [`lfa`](https://bioconductor.org/packages/release/bioc/html/lfa.html).
  - `jackstraw_lfa` and `jackstraw_alstructure` requires [`gcatest`](https://bioconductor.org/packages/release/bioc/html/gcatest.html).
  - `pip` requires the package [`qvalue`](https://bioconductor.org/packages/release/bioc/html/qvalue.html).

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

The stable version **jackstraw v1.3.17** is on CRAN. To install from CRAN:
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

