
library(readr)


colTypes <- cols(
  .default = col_double(),
  `X` = col_character()
)


expFile <- "/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/Filtered_Matrices/TPM_log2_QNorm_QCed_CovCorrected_AllCovariates.txt"

table_tmp <- read_delim(expFile, delim = "\t", quote = "", col_types = colTypes)
exp <- as.matrix(table_tmp[,-1])
rownames(exp) <- table_tmp[,1][[1]]
rm(table_tmp)

str(exp)

#save(exp, file = "/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/Filtered_Matrices/TPM_log2_QNorm_QCed_CovCorrected_AllCovariates.RData")
load(file = "/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/Filtered_Matrices/TPM_log2_QNorm_QCed_CovCorrected_AllCovariates.RData")

#exp contains expression rows genes cols samples

#First center and scale each row
#https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
expScale = exp - rowMeans(exp);
# Standardize each variable
expScale = expScale / sqrt(rowSums(expScale^2));   
expCov = tcrossprod(expScale);#equivalent to correlation due to center scale

expEigen <- eigen(expCov)

eigenVectors <- expEigen$vectors
colnames(eigenVectors) <- paste0("PC_",1:ncol(eigenVectors))
rownames(eigenVectors) <- rownames(expScale)
str(eigenVectors)

eigenValues <- expEigen$values
names(eigenValues) <- paste0("PC_",1:length(eigenValues))
str(eigenValues)

save(eigenVectors, eigenValues, expFile, file = "/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/PCA_Patrick/eigen.RData")

#Here calculate sample principle components. Number needed is arbritary (no more than eigen vectors)
expPcs <- t(expScale) %*% expEigen$vectors[,1:1000]

colnames(expPcs) <- paste0("PC_",1:ncol(expPcs))
str(expPcs)

save(expPcs, expFile, file = "/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/PCA_Patrick/pcs.RData")



library("havok")

medianSingularValue <- sqrt(median(eigenValues))

omega <- optimal_SVHT_coef(nrow(exp) / ncol(exp), sigma_known = F)
threshold <- omega * medianSingularValue

(numberComponentsToInclude <- sum(sqrt(eigenValues) > threshold ))

explainedVariance <- eigenValues * 100 / sum(eigenValues)

cumsum(explainedVariance)[numberComponentsToInclude]

which(cumsum(explainedVariance)>= 80)[1]

rpng()
plot(cumsum(explainedVariance), ylab = "Cumulative % explained variance")
abline(v = numberComponentsToInclude, col = "red3", lwd = 3)
dev.off()


rpng()
plot(cumsum(explainedVariance)[1:10000], ylab = "Cumulative % explained variance")
abline(v = numberComponentsToInclude, col = "red3", lwd = 3)
dev.off()

rpng()
plot(explainedVariance[1:numberComponentsToInclude], ylab = "% explained variance")
abline(v = numberComponentsToInclude, col = "red3", lwd = 3)
dev.off()



rpng()
plot(log(eigenValues), ylab = "Log eigenvalues")
abline(v = numberComponentsToInclude, col = "red3", lwd = 3)
dev.off()


rpng()
plot(log(sqrt(eigenValues)), ylab = "Signular values")
abline(h = log())
abline(v = numberComponentsToInclude, col = "red3", lwd = 3)
dev.off()


rpng()
plot(log(eigenValues)[1:5000], ylab = "Eigenvalues")
abline(v = numberComponentsToInclude, col = "red3", lwd = 3)
dev.off()

#Below is compare to genenetwork pipeline PCA


colTypes <- cols(
  .default = col_double(),
  `12/08/2022` = col_character()
)



table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/PCA_writeS/pc-scores.txt.gzip", delim = "\t", quote = "", col_types = colTypes)
pcsPipeline <- as.matrix(table_tmp[,-1])
rownames(pcsPipeline) <- table_tmp[,1][[1]]
rm(table_tmp)

str(pcsPipeline)

table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/PCA_writeS/eigenvectors.txt.gzip", delim = "\t", quote = "", col_types = colTypes)
eigenVecPipeline <- as.matrix(table_tmp[,-1])
rownames(eigenVecPipeline) <- table_tmp[,1][[1]]
rm(table_tmp)

variancePipeline <- read.csv(file = "/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/PCA_writeS/variances_and_singular_values.csv")
str(variancePipeline)

rpng()
plot(expPcs[,3], pcsPipeline[,3])
dev.off()


rpng()
layout(matrix(1:4,nrow = 2, byrow = TRUE))
for(i in 10^(0:3))
plot(expPcs[,i], pcsPipeline[,i], main = paste0("PC_",i), xlab = "PC in R", ylab = "PC in pipeline")
dev.off()




rpng()
layout(matrix(1:4,nrow = 2, byrow = TRUE))
for(i in 10^(0:3))
plot(expEigen$vectors[,i], eigenVecPipeline[,i], main = paste0("PC_",i), xlab = "Eigen vector in R", ylab = "Eigen vector in pipeline")
dev.off()


rpng()
plot(expEigen$values[1:1000], variancePipeline$Explained_variance, xlab = "Eigen values in R", ylab = "Eigen values in pipeline")
dev.off()

rpng()
plot(log(expEigen$values[1:1000]), log(variancePipeline$Explained_variance), xlab = "Eigen values in R", ylab = "Eigen values in pipeline")
dev.off()


rpng()
plot(log(expEigen$values))
dev.off()
