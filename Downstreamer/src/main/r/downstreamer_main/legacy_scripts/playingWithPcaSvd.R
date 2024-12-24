setwd("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\compareDecomposeMethods")

load("testExpression.RData", verbose = T)


str(expSub)




genes <- sample(nrow(expSub), 1000)
samples <- sample(ncol(expSub), 2000)

expSub2 <- expSub[genes, samples]


write.table(expSub2, file = "testData.txt", sep = "\t", quote = FALSE, col.names = NA)

str(expSub2)

#https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
expSubScale = expSub2 - rowMeans(expSub2);
# Standardize each variable
expSubScale = expSubScale / sqrt(rowSums(expSubScale^2));   

write.table(expSubScale, file = "testDataScaled2.txt", sep = "\t", quote = FALSE, col.names = NA)


hist(apply(expSubScale, 1, mean))
hist(apply(expSubScale, 1, sd))


expSvd <- svd(expSubScale, nu = 1000, nv = 1000)

reconstructedSvd <- expSvd$u %*% diag(expSvd$d) %*% t(expSvd$v)
str(reconstructedSvd)
str(expSubScale)

str(expSvd)

sum(expSvd$d^2)

plot(log(expSvd$d))

plot(log(sqrt(expSvd$d)))

sum(sqrt(expSvd$d))


sum((expSvd$d ^ 2) / (1000 - 1))
head((expSvd$d ^ 2) / (1000 - 1))

expCovNormal <- cov(t(expSubScale))
expCor <- cor(t(expSubScale))
expCov = tcrossprod(expSubScale);

range(expCov)
range(expCor)


plot(expCor[,1], expCov[,1])

expEigen <- eigen(expCov)
str(expEigen)

plot(sqrt(expEigen$values[1:1000]), expSvd$d)

plot(expEigen$values[1:1000], expSvd$d^2)

plot(expEigen$vectors[,1]*-1, expSvd$u[,1])
plot(expEigen$vectors[,99]*-1, expSvd$u[,99])
plot(expEigen$vectors[,1000], expSvd$u[,1000])

expPcs <- t(expSubScale) %*% expEigen$vectors


plot(diag(a), diag(b))

str(v)

v <- expPcs[,1:500] %*% diag(1/sqrt(expEigen$values[1:500]))
str(v)
str(expSvd$v)

plot(v[,1], expSvd$v[,1])
plot(v[,100], expSvd$v[,100])

test <- expEigen$vectors[,1:500] %*% diag(sqrt(expEigen$values[1:500])) %*% t(v[1,1:500,drop=F])
str(test)
str(expSubScale)

plot(test[,1], expSubScale[,1])
cor.test(test[,1], expSubScale[,1])
plot(test[1,], expSubScale[1,])
cor.test(test[1,], expSubScale[1,])

dim(expPcs)
dim(expEigen$values)

expPcsSub <- t(expSubScale) %*% expEigen$vectors[,1:100]

str(expPcsSub)
str(expPcs)

plot(expPcsSub[,1], expPcs[,1])
plot(expPcsSub[,100], expPcs[,100])


heatmap3::heatmap3(expEigen$vectors, scale = "none", Rowv = NA, Colv = NA, balanceColor = T)
str(expPcs)

range(expEigen$vectors[,-(1:100)])

str(expEigen$vectors)

str(expSvd)
expPcsSvd <- expSvd$v %*% diag(expSvd$d[1:1000])


str(expPcsSvd)
str(expPcsSvd2)
plot(expPcsSvd, expPcsSvd2)

plot(expPcs[,1]*-1, expPcsSvd[,1])
plot(expPcs[,99]*-1, expPcsSvd[,99])
plot(expPcs[,1000]*-1, expPcsSvd[,1000])

plot(abs(expPcs[2000,]), abs(expPcsSvd[2000,]))






library(readr)
table_tmp <- read_delim("TestWithFloranne3/pc-scores.txt.gzip", delim = "\t", quote = "")
GN_pcs <- as.matrix(table_tmp[,-1])
rownames(GN_pcs) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("TestWithFloranne3/eigenvectors.txt.gzip", delim = "\t", quote = "")
GN_eigen <- as.matrix(table_tmp[,-1])
rownames(GN_eigen) <- table_tmp[,1][[1]]
rm(table_tmp)

GN_variance <- read.delim("TestWithFloranne3/explained_variance.csv", header = FALSE)[,1]
GN_singularValues <- read.delim("TestWithFloranne3/singular_values.csv", header = FALSE)[,1]



plot(GN_pcs[,1], expPcsSvd[,1])
plot(GN_pcs[,99], expPcsSvd[,99])
plot(GN_pcs[,1000]*-1, expPcsSvd[,1000])

plot(GN_eigen[,1], expSvd$u[,1])
plot(GN_eigen[,1000]*-1, expSvd$u[,1000])

plot(GN_variance, expSvd$d^2)
plot(GN_variance, expEigen$values[1:100])


plot(GN_variance, (expSvd$d ^ 2) / (2000 - 1))

plot(expSvd$d, GN_singularValues)
plot(expSvd$d, GN_variance)


######################


expPca <- princomp(t(expSubScale), cor = F)#center = F, scale. = F

str(expPca)

sum(expPca$sdev)
plot(expPca$sdev^2, (expSvd$d^2) / (2000 - 1))

plot(expPca$sdev, sqrt(expEigen$values))


plot(expPca$loadings[,1]*-1, expEigen$vectors[,1])
str(expEigen$vectors)
str(expPca$scores)

plot(expPca$scores[,1], expPcs[,1]*-1)
plot(expPca$scores[,100], expPcs[,100]*-1)


plot(expPca$scores[,1], expPcsSvd[,1])
plot(expPca$scores[,100], expPcsSvd[,100]*-1)



str(expPca$rotation)
str(expPcsSvd)



#below to compare genenetwork results with randomizer
library(readr)
table_tmp <- read_delim("pc-scores.txt.gzip", delim = "\t", quote = "")
GN_pcs <- as.matrix(table_tmp[,-1])
rownames(GN_pcs) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("eigenvectors.txt.gzip", delim = "\t", quote = "")
GN_eigen <- as.matrix(table_tmp[,-1])
rownames(GN_eigen) <- table_tmp[,1][[1]]
rm(table_tmp)

GN_variance <- read.delim("explained_variance.csv", header = FALSE)[,1]

str(GN_variance)

plot(expSvd$v[,1], GN_pcs[,1])
plot(expSvd$v[,10], GN_pcs[,10])

heatmap3::heatmap3(cor(expSvd$v, GN_pcs)^2, balanceColor = T, scale = "none", Rowv = NA, Colv = NA)

heatmap3::heatmap3(cor(expSvd$u, GN_eigen)^2, balanceColor = T, scale = "none", Rowv = NA, Colv = NA)

plot(GN_variance[1:10], expSvd$d[1:10])
sum(GN_variance[1:75])
heatmap3::heatmap3(cor(expEigen$vectors)[1:10,1:10], balanceColor = T, scale = "none", Rowv = NA, Colv = NA)


sum(sqrt(expSvd$d[1:75]))/sum(sqrt(expSvd$d))*100



plot(expSvd$u[,5], GN_eigen[,5])

plot(expSvd$u[,1], expSvd$u[,5])


#below to compare genenetwork results with full SVD
library(readr)
table_tmp <- read_delim("GN2/pc-scores.txt.gzip", delim = "\t", quote = "")
GN_pcs <- as.matrix(table_tmp[,-1])
rownames(GN_pcs) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("GN2/eigenvectors.txt.gzip", delim = "\t", quote = "")
GN_eigen <- as.matrix(table_tmp[,-1])
rownames(GN_eigen) <- table_tmp[,1][[1]]
rm(table_tmp)

GN_variance <- read.delim("GN2/explained_variance.csv", header = FALSE)[,1]

str(GN_variance)

plot(expSvd$v[,1], GN_pcs[,1])
plot(expSvd$v[,10], GN_pcs[,10])

heatmap3::heatmap3(cor(expSvd$v, GN_pcs)^2, balanceColor = T, scale = "none", Rowv = NA, Colv = NA)

heatmap3::heatmap3(cor(expSvd$u, GN_eigen)^2, balanceColor = T, scale = "none", Rowv = NA, Colv = NA)

plot(GN_variance[1:100], expSvd$d[1:100]^2)


plot(GN_variance, expPca$sdev[1:100]^2)

plot(expSvd$u[,5], GN_eigen[,5])

plot(expSvd$u[,1], expSvd$u[,5])

plot(log(expSvd$d))
abline(v=60)



library(corpcor)
expSvdFast <- fast.svd(expSubScale)
plot(expSvdFast$d,  expSvd$d)
plot(expSvdFast$u[,1],  expSvd$u[,1])
plot(expSvdFast$u[1,1:100],  expSvd$u[1,])
plot(expSvdFast$v[,1],  expSvd$v[,1])
plot(expSvdFast$v[1,1:100],  expSvd$v[1,])

str(expSvdFast)
