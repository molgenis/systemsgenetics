#srun --cpus-per-task=20 --mem=200gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, sync = T)




remoter::client("localhost", port = 55556)#55501

library("havok")
library(parallel)

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

#load(file = "perTissueNormalization/selectedSamplesRawExpression.RData", verbose = T)
load("tissuePredictions/samplesWithPrediction_16_09_22_noOutliers.RData", verbose = T)


load(file = "Metadata/combinedMeta_2022_09_15.RData", verbose = T)
combinedMeta <- combinedMeta[rownames(samplesWithPredictionNoOutliers),]

dim(combinedMeta)
dim(samplesWithPredictionNoOutliers)



selectedSamples <- rownames(samplesWithPredictionNoOutliers)


sraFiles <- list.files(path="rse-sra/SRA_Files/", pattern="sra*", full.names=TRUE, recursive=FALSE)
gtexFiles <- list.files(path="rse-gtex/rse_gtex", pattern="rse*", full.names=TRUE, recursive=FALSE)
allFiles <- c(sraFiles, gtexFiles, "rse-tcga/rseTCGA.rda", "rse-tcga/rse_ESCA_TCGA.rda")


file = gtexFiles[2]

perChunkExp <- sapply(allFiles, function(file){
  
  loadedObject <- load(file)
  
  sreObjects <- get(loadedObject[1])
  
  #sometimes single RSE is not in list. Put in list of one to make code uniform
  if(!is.list(sreObjects)){
    sreObjects <- list(sreObjects)
  }
  
  #sreObject <- sreObjects[[1]]
  
  perStudyExp <- lapply(sreObjects, function(sreObject){
    studyExp <- sreObject@assays@data@listData$TPM
    return(studyExp[,colnames(studyExp) %in% selectedSamples, drop = F])
  })
  
 # x <- do.call(cbind, perStudyExp)
  #table(notFound %in% colnames(x) )
  
  return(do.call(cbind, perStudyExp))
  
})


recountHealthyExp <- do.call(cbind, perChunkExp)
all(selectedSamples %in% colnames(recountHealthyExp ))
table(selectedSamples %in% colnames(recountHealthyExp ))



#Some samples are duplicated in the chunks, now make sure only one is in the matrix
uniqueSamplesIndex <- match(selectedSamples, colnames(recountHealthyExp))
recountHealthyExp <- recountHealthyExp[,uniqueSamplesIndex]
str(recountHealthyExp)



includedGenes <- apply(recountHealthyExp, 1, function(x){(sum(x==0)/ncol(recountHealthyExp)) <= 0.5})
table(includedGenes)


recountHealthyExp <- recountHealthyExp[includedGenes,]

recountHealthyExp <- log2(recountHealthyExp + 1)

#BiocManager::install("preprocessCore", configure.args="--disable-threading")
library(preprocessCore, lib.loc = .libPaths()[2])
normalize.quantiles(recountHealthyExp,copy=FALSE)

str(recountHealthyExp)


all(colnames(recountHealthyExp) == rownames(combinedMeta))

geneMean <- rowMeans(recountHealthyExp)





#TCGA samples don't have sra.sample_spots instead use recount_qc.bc_frag.count
missingSampleSpots <- is.na(combinedMeta[,"sra.sample_spots"])
combinedMeta[missingSampleSpots,"sra.sample_spots"] <- combinedMeta[missingSampleSpots,"recount_qc.bc_frag.count"]

#TCGA and Gtex don't report layout but is all paired
missingLayout <- is.na(combinedMeta[,"sra.library_layout"])
combinedMeta[missingLayout,"sra.library_layout"] <- "paired"

covariatesToCorrectFor <- read.delim("CovariateNames.txt", header = F)$V1
combinedMetaSelection <- combinedMeta[,covariatesToCorrectFor]
combinedMetaSelection$sra.library_layout <- as.factor(combinedMetaSelection$sra.library_layout)




cl <- makeCluster(5)
recountHealthyExpNorm <- parApply(cl, recountHealthyExp, 1 ,function(geneExp, combinedMetaSelection){
  return(residuals(lm(geneExp ~ . ,data = combinedMetaSelection)))
}, combinedMetaSelection = combinedMetaSelection)
stopCluster(cl)


recountHealthyExpNorm <- t(recountHealthyExpNorm)

recountHealthyExpNorm <- recountHealthyExpNorm + geneMean

#save(recountHealthyExpNorm, file = paste0("CombinedHealthyTissue/combinedHealthyTissue_TPM_Log2_QQ_CovCor_Exp.RData"))


#https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
expScale = recountHealthyExpNorm - rowMeans(recountHealthyExpNorm);
# Standardize each variable
expScale = expScale / sqrt(rowSums(expScale^2));   
expCov = tcrossprod(expScale);#equevelent to correlation due to center scale
expEigen <- eigen(expCov)
eigenVectors <- expEigen$vectors
colnames(eigenVectors) <- paste0("Comp_",1:ncol(eigenVectors))
rownames(eigenVectors) <- rownames(expScale)

eigenValues <- expEigen$values
names(eigenValues) <- paste0("Comp_",1:length(eigenValues))

#Here calculate sample principle components. Number needed is arbritary (no more than eigen vectors)
expPcs <- t(expScale) %*% expEigen$vectors[,1:50]
colnames(expPcs) <- paste0("Comp_",1:ncol(expPcs))

#nrSamples <- ncol(expScale)

#expSvd <- svd(expScale, nu = nrSamples, nv = min(nrSamples, 50))


#eigenValues <- expSvd$d^2
#eigenVectors <- expSvd$u
#colnames(eigenVectors) <- paste0("Comp_",1:ncol(eigenVectors))
#rownames(eigenVectors) <- rownames(expScale)

#expPcs <- expSvd$v %*% diag(expSvd$d[1:ncol(expSvd$v)])
#colnames(expPcs) <- paste0("Comp_",1:ncol(expPcs))
#rownames(expPcs) <- colnames(expScale)

explainedVariance <- eigenValues * 100 / sum(eigenValues)


combinedHealtyTissuePca <- list(eigenVectors = eigenVectors, eigenValues = eigenValues, expPcs = expPcs, explainedVariance = explainedVariance)

#save(combinedHealtyTissuePca, file = paste0("CombinedHealthyTissue/combinedHealthyTissue_PCA.RData"))

medianSingularValue <- median(sqrt(eigenValues))

omega <- optimal_SVHT_coef(nrow(expScale) / ncol(expScale), sigma_known = F)
threshold <- omega * medianSingularValue
numberComponentsToInclude <- sum(sqrt(eigenValues) > threshold )

h <- cumsum(explainedVariance)[numberComponentsToInclude ]


png(paste0("CombinedHealthyTissue/explainedVar.png"),width = 1000, height = 1000)
plot(cumsum(explainedVariance), pch = 16, cex = 0.5, xlab = "Component", ylab = "Cumulative explained %")
abline(h = h, lwd = 2, col = "darkred")
text(0,h+1,numberComponentsToInclude, adj = 0)
dev.off()

#write.table(eigenVectors[,1:numberComponentsToInclude], file = gzfile(paste0("CombinedHealthyTissue/combinedHealthyTissue_eigenVec.txt.gz")), sep = "\t", quote = F, col.names = NA)






str(samplesWithPredictionNoOutliers)
tissueClasses <- unique(samplesWithPredictionNoOutliers$predictedTissue     )

tissue <- tissueClasses[1]
medianPerTissueList <- lapply(tissueClasses, function(tissue){
  
  tissueSamples <- rownames(samplesWithPredictionNoOutliers)[samplesWithPredictionNoOutliers$predictedTissue == tissue]
  
  tissueExp <- recountHealthyExpNorm[,tissueSamples]

  tissueMedian <- apply(tissueExp, 1, median)
  
  return(tissueMedian)
  
})
str(medianPerTissueList)

medianPerTissue <- do.call(cbind, medianPerTissueList)
colnames(medianPerTissue) <- tissueClasses

colHeatmap <- rev(c(colorRampPalette(c("#f03b20", "#feb24c", "#ffeda0"))(99), rep("white",40), colorRampPalette(c("#e0ecf4", "#9ebcda", "#8856a7"))(99)))
colBreaks <- seq(-1,1,length.out= length(colHeatmap)+1)

shortTissue <- ifelse(nchar(tissueClasses) > 20, paste0(substr(tissueClasses,0,17),"..."), tissueClasses)



library(pheatmap)
library(viridisLite, lib.loc = .libPaths()[2])



medianPerTissueCor <- cor(medianPerTissue)
diag(medianPerTissueCor) <- NA
clust <- hclust(as.dist(1 - medianPerTissueCor))
rpng(width = 1000, height = 1000)
pheatmap(medianPerTissueCor, scale = "none", labels_row = shortTissue, labels_col = shortTissue, cluster_rows = clust, cluster_cols = clust, border_color = NA)
dev.off()

medianPerTissueScaled <- t(scale(t(medianPerTissue)))
medianPerTissueScaledCor <- cor(medianPerTissueScaled)
diag(medianPerTissueScaledCor) <- 0
clust2 <- hclust(as.dist(1 - medianPerTissueScaledCor))
rpng(width = 1000, height = 1000)
pheatmap(medianPerTissueScaledCor, scale = "none", col = colHeatmap, breaks = colBreaks, labels_row = shortTissue, labels_col = shortTissue, cluster_rows = clust2, cluster_cols = clust2, border_color = NA)
dev.off()


medianOfmedianPerTissue <- apply(medianPerTissue, 1, median)

medianOfmedianPerTissueUnique <- medianPerTissue - medianOfmedianPerTissue
medianOfmedianPerTissueUniqueCor <- cor(medianOfmedianPerTissueUnique)
diag(medianOfmedianPerTissueUniqueCor) <- 0
clust3 <- hclust(as.dist(1 - medianOfmedianPerTissueUniqueCor))
rpng(width = 1000, height = 1000)
pheatmap(medianOfmedianPerTissueUniqueCor, scale = "none", col = colHeatmap, breaks = colBreaks, labels_row = shortTissue, labels_col = shortTissue, cluster_rows = clust3, cluster_cols = clust3, border_color = NA)
dev.off()
