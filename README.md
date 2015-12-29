# jackstraw

This **R** package performs association tests between variables and estimated latent variables.
Latent variables are unobserved and must be estimated directly from the data, using various techniques such as principal component analysis, logistic factor analysis, and others.
However, treating estimated latent variables as independent variables in a conventional regression framework would results in an anti-conservative bias (i.e., artificially inflated significance).
The jackstraw method account for this fact that latent variables are estimated from the data and to protect association tests from an anti-conservative bias.

*Chung, N.C.* and *Storey, J.D.*	(2015) Statistical significance of variables driving systematic variation in high-dimensional data	Bioinformatics, 31(4): 545-554
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
