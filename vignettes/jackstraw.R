## ----setup,include=FALSE,echo=FALSE-------------
opts_chunk$set(fig.align='center', fig.width=6, fig.height=6, tidy=TRUE, cache=FALSE, warning=FALSE, message=TRUE)
options(keep.source = TRUE, width=50)
desc <- packageDescription("jackstraw")

## ----sim_pca_data-------------------------------
library(jackstraw)
library(corpcor)

set.seed(1)
B = c(runif(100, min=0.1, max=1), rep(0,900))
L = c(rep(1, 10), rep(-1, 10))
L = L / sd(L)
E = matrix(rnorm(1000*20), nrow=1000)
Y = B %*% t(L) + E

dim(Y)
Y[1:5,1:5]

## ----sim_pca_PA, dependson="sim_pca_data"-------
PA = permutationPA(Y, B=10, threshold=0.05)

plot(PA$p, pch=20, main="Permutation Parallel Analysis P-values", ylab="P-values", xlab="Principal Component")

## ----sim_pca_1pc, dependson="sim_pca_data"------
svd.out = fast.svd(Y)

par(mfrow=c(2,1))
plot(svd.out$d^2/sum(svd.out$d^2), pch=20, main="The scree plot", xlab="PC", ylab="Percent Variance Explained")
plot(svd.out$d[1] * svd.out$v[,1], pch=20, main="1st PC", xlab="Observation", ylab="Magnitude")

## ----sim_pca_pvalues, dependson="sim_pca_jackstraw"----
par(mfrow=c(1,2))
hist(js.pca$p.value[1:100], 10, col="black", main="Alternative P-values")
hist(js.pca$p.value[101:1000], 10, col="black", main="Null P-values")

## ----sim_lfa_jackstraw, dependson="sim_lfa_data", cache=TRUE----
js.lfa = jackstraw_lfa(dat$Y, r = 2, FUN = function(x) lfa.corpcor(x, 2)[, , drop = FALSE], s = 200, B = 10, devR = TRUE)


hist(js.lfa$p.value, 10, col="black")

## ----sim_lfa_pvalues, dependson="sim_lfa_jackstraw"----
par(mfrow=c(1,2))
hist(js.lfa$p.value[which(dat$H==1)], 10, col="black", main="Alternative P-values")
hist(js.lfa$p.value[which(dat$H==0)], 10, col="black", main="Null P-values")

