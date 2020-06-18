setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

library(readr)
table_tmp <- read_delim("pheno_9_66/pheno_9_geneMaxSnpScores.txt", delim = "\t", quote = "")
geneMaxSnpZscores <- as.matrix(table_tmp[,-1])
rownames(geneMaxSnpZscores) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("pheno_9_66/pheno_9_genePvalues.txt", delim = "\t", quote = "")
genePvalues <- as.matrix(table_tmp[,-1])
rownames(genePvalues) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("pheno_9_66/pheno_9_geneMaxSnpZscoresNullGwas.txt", delim = "\t", quote = "")
geneMaxSnpZscoresNullGwas <- as.matrix(table_tmp[,-1])
rownames(geneMaxSnpZscoresNullGwas) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("pheno_9_66/pheno_9_genePvaluesNullGwas.txt", delim = "\t", quote = "")
genePvaluesNullGwas <- as.matrix(table_tmp[,-1])
rownames(genePvaluesNullGwas) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("pheno_9_66/random_pathway_database.txt", delim = "\t", quote = "")
pathways <- as.matrix(table_tmp[,-1])
rownames(pathways) <- table_tmp[,1][[1]]
rm(table_tmp)


table_tmp <- read_delim("pheno_9_66/pheno_9_debugFiles/random_12_q_Enrichment_geneCorrelationsForMetaGenes_ExHla.txt", delim = "\t", quote = "")
ds_cor12qMeta <- as.matrix(table_tmp[,-1])
rownames(cor12qMeta) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("pheno_9_66/pheno_9_debugFiles/random_12_q_Enrichment_geneZscoresForMetaGenes_ExHla.txt", delim = "\t", quote = "")
sd_zScorecor12qMeta <- as.matrix(table_tmp[,-1])
rownames(zScorecor12qMeta) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("pheno_9_66/pheno_9_B_debugFiles/random_Enrichment_geneZscoresForMetaGenes_ExHla.txt", delim = "\t", quote = "")
ds_geneZscoresForCor <- as.matrix(table_tmp[,-1])
rownames(ds_geneZscoresForCor) <- table_tmp[,1][[1]]
rm(table_tmp)


table_tmp <- read_delim("pheno_9_66/pheno_9_B_debugFiles/random_12_q_Enrichment_geneCor.txt", delim = "\t", quote = "")
ds_random_12_q_Enrichment_geneCor <- as.matrix(table_tmp[,-1])
rownames(ds_random_12_q_Enrichment_geneCor) <- table_tmp[,1][[1]]
rm(table_tmp)


table_tmp <- read_delim("pheno_9_66/pheno_9_debugFiles/random_Enrichment_normalizedNullGwasGeneScoresCor_ExHla.txt", delim = "\t", quote = "")
ds_normalizedNullGwasGeneScoresCor <- as.matrix(table_tmp[,-1])
rownames(ds_normalizedNullGwasGeneScoresCor) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("pheno_9_66/pheno_9_debugFiles/random_Enrichment_invCorMatrix_ExHla.txt", delim = "\t", quote = "")
ds_invCorMatrix <- as.matrix(table_tmp[,-1])
rownames(ds_invCorMatrix) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("pheno_9_66/pheno_9_B_debugFiles/random_Enrichment_invCorMatrix_ExHla.txt", delim = "\t", quote = "")
ds_invCorMatrix_B <- as.matrix(table_tmp[,-1])
rownames(ds_invCorMatrix_B) <- table_tmp[,1][[1]]
rm(table_tmp)

geneZscores <- qnorm((genePvalues/2))
geneZscoresNullGwas <- qnorm((genePvaluesNullGwas/2))

sharedGenes <- read.delim("pheno_9_66/pheno_9_debugFiles/random_Enrichment_sharedGenes.txt",stringsAsFactors = F)$Gene

genes <- read.delim("pheno_9_66/ensgR75_protein_coding_subset_12q.txt", stringsAsFactors = F)

sharedGenes <- intersect(genes$Ensembl.Gene.ID, sharedGenes)

geneMaxSnpZscores <- geneMaxSnpZscores[sharedGenes,,drop=F]
genePvalues <- genePvalues[sharedGenes,,drop=F]
geneMaxSnpZscoresNullGwasCorr <- geneMaxSnpZscoresNullGwas[sharedGenes,1:10000,drop=F]
geneMaxSnpZscoresNullGwasFdr <- geneMaxSnpZscoresNullGwas[sharedGenes,10001:10100,drop=F]
genePvaluesNullGwasCorr <- genePvaluesNullGwas[sharedGenes,1:10000,drop=F]
genePvaluesNullGwasFdr <- genePvaluesNullGwas[sharedGenes,10001:10100,drop=F]
pathways <- pathways[sharedGenes,,drop=F]

geneZscoreNullGwasCorr <- qnorm((genePvaluesNullGwasCorr/2)) * -1
geneZscores <- qnorm((genePvalues/2)) * -1

geneMaxSnpZscoresNullGwasCorrScaled <- apply(geneZscoreNullGwasCorr, 2, scale)
rownames(geneMaxSnpZscoresNullGwasCorrScaled) <- rownames(geneZscoreNullGwasCorr)

geneMaxSnpZscoresNullGwasCorrCor <- cor(t(geneMaxSnpZscoresNullGwasCorrScaled))


rThres <- 0.5



metaGenesMap <- rep(0,length(sharedGenes))
names(metaGenesMap) <- sharedGenes
metaGenes <- as.list(as.vector(rep(NA, length(sharedGenes)), mode = "character"))
metaGeneCount = 0;

for(g in 1:(length(sharedGenes))){
#for(g in 1:5){
  
  currentMetaGeneId <- NA
  
  if(metaGenesMap[g] == 0){
    metaGeneCount <- metaGeneCount + 1
    currentMetaGeneId <- metaGeneCount
    metaGenesMap[g] <- currentMetaGeneId
    metaGenes[[currentMetaGeneId]] <- sharedGenes[g]
  } else {
    currentMetaGeneId <- metaGenesMap[g]
  }
  
  if(g < length(sharedGenes)){
    for(g2 in ((g+1):length(sharedGenes))){
      
      
      if(geneMaxSnpZscoresNullGwasCorrCor[g,g2] >= rThres){
        if(metaGenesMap[g2] == 0){
          metaGenes[[currentMetaGeneId]] <- c(metaGenes[[currentMetaGeneId]], sharedGenes[g2])
          metaGenesMap[g2] <- currentMetaGeneId
        } else if(metaGenesMap[g2] != currentMetaGeneId){
          otherMetaGene <- metaGenes[[metaGenesMap[g2]]]
          if(metaGenesMap[g2] == 2){
            print(otherMetaGene)
            print(currentMetaGeneId)
          }
          
          metaGenes[[metaGenesMap[g2]]] <- NA
          metaGenes[[currentMetaGeneId]] <- c(metaGenes[[currentMetaGeneId]], otherMetaGene)
          metaGenesMap[otherMetaGene] <- currentMetaGeneId
          
        }
      }
    }
  }
}


metaGenes <- sapply(metaGenes, sort)
res <- sapply(metaGenes[unique(metaGenesMap)], function(y){return(paste(y,sep="_",collapse="_"))})
cat(res, sep = "\n", file = "test.txt")


geneMaxSnpZscoresNullGwasCorrScaledMetaList <- lapply(metaGenes[unique(metaGenesMap)], function(metaGene){
  
  return(apply(geneMaxSnpZscoresNullGwasCorrScaled[metaGene, ,drop =F], 2, sum) / length(metaGene))
  
})

geneMaxSnpZscoresNullGwasCorrScaledMeta <- do.call(rbind, geneMaxSnpZscoresNullGwasCorrScaledMetaList)

rownames(geneMaxSnpZscoresNullGwasCorrScaledMeta) <- sapply(metaGenes[unique(metaGenesMap)], function(y){return(paste(y,sep="_",collapse="_"))})

geneMaxSnpZscoresNullGwasCorrScaledMetaScaled <- apply(geneMaxSnpZscoresNullGwasCorrScaledMeta, 2, scale)
rownames(geneMaxSnpZscoresNullGwasCorrScaledMetaScaled)  <- rownames(geneMaxSnpZscoresNullGwasCorrScaledMeta) 

str(geneMaxSnpZscoresNullGwasCorrScaledMetaScaled)

str(ds_normalizedNullGwasGeneScoresCor)

all(rownames(ds_normalizedNullGwasGeneScoresCor) %in% rownames(geneMaxSnpZscoresNullGwasCorrScaledMeta))

#plot(geneMaxSnpZscoresNullGwasCorrScaledMetaScaled[rownames(ds_normalizedNullGwasGeneScoresCor),], ds_normalizedNullGwasGeneScoresCor)

range(geneMaxSnpZscoresNullGwasCorrScaledMetaScaled[rownames(ds_normalizedNullGwasGeneScoresCor),] - ds_normalizedNullGwasGeneScoresCor)

geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor <- cor(t(geneMaxSnpZscoresNullGwasCorrScaledMetaScaled))
str(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor)

library(heatmap3)
library(RColorBrewer)
heatmap3(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor, scale= "none", Rowv = NA, Colv = NA, col = colorRampPalette( brewer.pal(9, "YlOrBr"))(100))

tmp<-geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor
diag(tmp) <- 0
range(tmp)

tmp2 <- cor(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor)
diag(tmp) <- 0
range(tmp)

geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor[abs(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor) <= 0.01] <- 0

str(ds_random_12_q_Enrichment_geneCor)
str(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor)

plot(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor[rownames(ds_random_12_q_Enrichment_geneCor), colnames(ds_random_12_q_Enrichment_geneCor)], ds_random_12_q_Enrichment_geneCor)
range(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor[rownames(ds_random_12_q_Enrichment_geneCor), colnames(ds_random_12_q_Enrichment_geneCor)] - ds_random_12_q_Enrichment_geneCor)

heatmap3(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCorInv, scale= "none", Rowv = NA, Colv = NA, col = colorRampPalette( brewer.pal(9, "YlOrBr"))(100))
heatmap3(ds_invCorMatrix[rownames(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCorInv), colnames(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCorInv)], scale= "none", Rowv = NA, Colv = NA, col = colorRampPalette( brewer.pal(9, "YlOrBr"))(100))
heatmap3(ds_invCorMatrix_B[rownames(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCorInv), colnames(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCorInv)], scale= "none", Rowv = NA, Colv = NA, col = colorRampPalette( brewer.pal(9, "YlOrBr"))(100))

geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCorInv <- solve(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor)

range(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCor[rownames(ds_geneZscoresForCor),colnames(ds_geneZscoresForCor)] - ds_random_12_q_Enrichment_geneCor)

plot(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCorInv[rownames(ds_invCorMatrix),colnames(ds_invCorMatrix)], ds_invCorMatrix)

range(geneMaxSnpZscoresNullGwasCorrScaledMetaScaledCorInv[rownames(ds_invCorMatrix_B),colnames(ds_invCorMatrix_B)] - ds_invCorMatrix_B)
