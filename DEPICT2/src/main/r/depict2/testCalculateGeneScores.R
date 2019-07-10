
setwd("C:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

trait = "Alzheimers_v29/Alzheimers"

library(readr)


table_tmp <- read_delim(paste0(trait, "_eigen_Enrichment_betasExHla.txt"), delim = "\t", quote = "")
eigenBeta<- as.matrix(table_tmp[,-1])
rownames(eigenBeta) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim(paste0(trait,"_eigen_Enrichment_normalizedPathwayScores_ExHla.txt"), delim = "\t", quote = "")
eigenNorm <- as.matrix(table_tmp[,-1])
rownames(eigenNorm) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim(paste0(trait,"_eigen_Enrichment_normalizedGwasGeneScores_ExHla.txt"), delim = "\t", quote = "")
zscoresNorm <- as.matrix(table_tmp[,-1])
rownames(zscoresNorm) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim(paste0(trait,"_genePvalues.txt"), delim = "\t", quote = "")
genePvalues <- as.matrix(table_tmp[,-1])
rownames(genePvalues) <- table_tmp[,1][[1]]
rm(table_tmp)

geneScores <- apply(eigenNorm, 1, function(row){row %*% eigenBeta[,1]})

plot(zscoresNorm, geneScores, xlab = "Centered and scaled gene Z-score", ylab = "Reconstructed gene scores")
cor.test(zscoresNorm, geneScores)



write.table(cbind(zscoresNorm, geneScores), file = paste0(trait,"_geneScoresExHla.txt"), sep = "\t", quote = F)
