
library(readr)


colTypes <- cols(
  .default = col_double(),
  `x` = col_character()
)



table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/Filtered_Matrices/TPM_log2_QNorm_QCed_CovCorrected_AllCovariates.txt", delim = "\t", quote = "", col_types = colTypes)
exp <- as.matrix(table_tmp[,-1])
rownames(exp) <- table_tmp[,1][[1]]
rm(table_tmp)

str(exp)

save(exp, file = "/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/Filtered_Matrices/TPM_log2_QNorm_QCed_CovCorrected_AllCovariates.RData")


#https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
expScale = exp - rowMeans(exp);
# Standardize each variable
expScale = expScale / sqrt(rowSums(expScale^2));   
expCov = tcrossprod(expScale);#equevelent to correlation due to center scale
range(expCov)
str(expCov)

expEigen <- eigen(expCov)
expPcs <- t(expScale) %*% expEigen$vectors

