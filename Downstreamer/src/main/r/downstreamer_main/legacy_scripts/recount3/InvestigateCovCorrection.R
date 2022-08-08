#srun --cpus-per-task=1 --mem=200gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)



remoter::client("localhost", port = 55508, password = "laberkak")

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

library(readr)


colTypes <- cols(
  .default = col_double(),
  `X1` = col_character()
)



table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/TPM_log2_QQnorm.txt", delim = "\t", quote = "", col_types = colTypes)
tpmLog2 <- as.matrix(table_tmp[,-1])
rownames(tpmLog2) <- table_tmp[,1][[1]]
rm(table_tmp)

save(tpmLog2, file = "ExploreCovariateCorrection/tpmLog2.RData")

table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/Filtered_Matrices/TPM_log2_QNorm_QCed_CovCorrected_CovariatesInGN.txt", delim = "\t", quote = "", col_types = colTypes)
tpmLog2Cov <- as.matrix(table_tmp[,-1])
rownames(tpmLog2Cov) <- table_tmp[,1][[1]]
rm(table_tmp)

save(tpmLog2Cov, file = "ExploreCovariateCorrection/tpmLog2Cov.RData")


load("ExploreCovariateCorrection/tpmLog2.RData")
load("ExploreCovariateCorrection/tpmLog2Cov.RData")


covariates <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/Covariate_Filtration/Cov_Correction.txt", delim = "\t", quote = "")

dim(covariates)
str(covariates)

samplesPastQc <- colnames(tpmLog2Cov)


tpmLog2 <- tpmLog2[,samplesPastQc]

dim(tpmLog2)
dim(tpmLog2Cov)

covariates <- covariates[match(samplesPastQc,covariates$samples),]
str(covariates)

all(covariates$samples == colnames(tpmLog2))

covariates <- covariates[,-1]
rownames(covariates) <- covariates$samples
str(tpmLog2)
cor.test(covariates$read_count, tpmLog2[1,])
rpng()
plot(covariates$read_count, tpmLog2[1,])
dev.off()


cor.test(covariates$read_count, tpmLog2Cov[1,])
rpng()
plot(covariates$read_count, tpmLog2Cov[1,])
dev.off()


rpng()
plot(tpmLog2[1,], tpmLog2Cov[1,])
dev.off()


rpng()
plot(covariates$deletion, tpmLog2Cov[2,])
dev.off()



rpng()
plot(tpmLog2Cov[3,], tpmLog2Cov[1,])
dev.off()

test <- lm(tpmLog2[1,] ~ covariates$Mapped_reads)
anova(test)
x <- residuals(test)
rpng()
plot(tpmLog2[1,], x)
dev.off()

range(tpmLog2[1,])
range(x)

rpng()
plot(tpmLog2[1000,], tpmLog2Cov[1000,], xlab = "TPM LOG2 Qnorm", ylab = "TPM LOG2 Qnorm Cov correct")
dev.off()
cor.test(tpmLog2[4,], tpmLog2Cov[4,])
