
setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")


library(readr)

colTypes <- cols(
  .default = col_double(),
  `X1` = col_character()
)


table_tmp <- read_delim("Recount3_QC_2ndRun/Filtered_Matrices/TPM_log2.txt", delim = "\t", quote = "", col_types = colTypes)
tmpLog2 <- as.matrix(table_tmp[,-1])
rownames(tmpLog2) <- table_tmp[,1][[1]]
rm(table_tmp)
str(svd)

#save(tmpLog2, file = "Recount3_QC_2ndRun/Filtered_Matrices/TPM_log2.RData")

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

load(file = "TPM_log2.RData")

library(preprocessCore)


normalize.quantiles(tmpLog2,copy=FALSE)

write.table(tmpLog2, file = "TPM_log2_QQnorm.txt", sep = "\t", quote = FALSE, col.names = NA)

#tpmLog2QqCor <- cor(t(tmpLog2))

#save(tpmLog2QqCor, file = "tpmLog2QqCor.RData")

#expSvd <- svd(tmpLog2, nu = 10, nv = 10)

#save(expSvd, file = "qqnormSvd.RData")
