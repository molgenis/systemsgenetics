#srun --cpus-per-task=20 --mem=200gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, sync = T)




remoter::client("localhost", port = 55504)

library(DESeq2)
library(parallel)
library(viridisLite, lib.loc = .libPaths()[2])
library(preprocessCore)




setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Recount3\\")

load(file = "perTissueNormalization/selectedSamplesRawExpression.RData", verbose = T)
load("tissuePredictions/samplesWithPrediction_16_09_22.RData", verbose = T)

samplesWithPrediction$sra.library_layout[samplesWithPrediction$study=="TCGA"] <- "paired"
samplesWithPrediction$sra.library_layout[samplesWithPrediction$study=="GTEx"] <- "paired"



sort(table(samplesWithPrediction$predictedTissue))
tissueClasses <- unique(samplesWithPrediction$predictedTissue)

mclapply(tissueClasses,  mc.cores = 10, function(tissue){
  
  tissueSamples <- rownames(samplesWithPrediction)[samplesWithPrediction$predictedTissue == tissue]
  tissueExp <- selectedSamplesExp[,tissueSamples]
  numberOfSamples <- length(tissueSamples)
  
  includedGenes <- apply(tissueExp, 1, function(x){(sum(x==0)/numberOfSamples) <= 0.5})
  
  tissueExp <- tissueExp[includedGenes,]
  
  tissueExp <- log2(tissueExp + 1)
  
  normalize.quantiles(tissueExp,copy=FALSE)

  save(tissueExp, file = paste0("perTissueNormalization/perTissueQq/",make.names(tissue),".RData"))
  
})

  
perTissuePca <- mclapply(tissueClasses,  mc.cores = 10, function(tissue){
  
  load(file = paste0("perTissueNormalization/perTissueQq/",make.names(tissue),".RData"))
  
  #https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
  expScale = tissueExp - rowMeans(tissueExp);
  # Standardize each variable
  expScale = expScale / sqrt(rowSums(expScale^2));   
  #expCov = tcrossprod(expScale);#equevelent to correlation due to center scale
  #expEigen <- eigen(expCov)
  #eigenVectors <- expEigen$vectors
  #colnames(eigenVectors) <- paste0("PC_",1:ncol(eigenVectors))
  #rownames(eigenVectors) <- rownames(expScale)
  
  #eigenValues <- expEigen$values
  #names(eigenValues) <- paste0("PC_",1:length(eigenValues))
  
  #Here calculate sample principle components. Number needed is arbritary (no more than eigen vectors)
  #expPcs <- t(expScale) %*% expEigen$vectors[,1:10]
  #colnames(expPcs) <- paste0("PC_",1:ncol(expPcs))
  
  expSvd <- svd(expScale, nu = 50, nv = 50)
  
  eigenValues <- expSvd$d^2
  eigenVectors <- expSvd$u
  colnames(eigenVectors) <- paste0("PC_",1:ncol(eigenVectors))
  rownames(eigenVectors) <- rownames(expScale)
  
  expPcs <- expSvd$v[,1:50] %*% diag(expSvd$d[1:50])
  colnames(expPcs) <- paste0("PC_",1:ncol(expPcs))
  rownames(expPcs) <- colnames(expScale)
  
  pcaRes <- list(eigenVectors, eigenValues, expPcs)
  
  save(pcaRes, file = paste0("perTissueNormalization/perTissueQqPca/",make.names(tissue),".RData"))
  
  return(pcaRes)
  
})


sink <- mclapply(tissueClasses,  mc.cores = 10, function(tissue, samplesWithPrediction){
  
  load(file = paste0("perTissueNormalization/perTissueQqPca/",make.names(tissue),".RData"))

  expPcs <- pcaRes$expPcs
  
  tissueSamplesInfo <- samplesWithPrediction[rownames(expPcs),]
  studies <- length(unique(tissueSamplesInfo$study))
  
  
  palette(adjustcolor(viridis(studies, option = "H"), alpha.f = 0.5))
  
  pchMap <- rep(c(15,16,17), length.out = studies)
  

  plot(expPcs[,1],expPcs[,2], col = as.factor(tissueSamplesInfo$study), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 1)
  
  
  return(NULL)
  
}, samplesWithPrediction = samplesWithPrediction)





















save(perTissueExp, perTissuePca, file = "perTissueNormalization/tmpTestSession.RData")
#load(file = "perTissueNormalization/tmpTestRlog.RData")

save(expPcs, samplesWithPrediction, file = "perTissueNormalization/tmpTest2.RData")
load("perTissueNormalization/tmpTest2.RData")



tissueSamplesInfo <- samplesWithPrediction[rownames(expPcs),]
str(tissueSamplesInfo)

#Put TCGA and GTEx to paired end

studies <- length(unique(tissueSamplesInfo$study))



palette(adjustcolor(viridis(studies, option = "H"), alpha.f = 0.5))

pchMap <- rep(c(15,16,17), length.out = studies)

rpng(width = 1000, height = 1000)
plot(expPcs[,1],expPcs[,2], col = as.factor(tissueSamplesInfo$study), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 1)
pairs(expPcs[,1:5], col = as.factor(tissueSamplesInfo$study), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 1, upper.panel = NULL)
dev.off()
View(tissueSamplesInfo)

palette(adjustcolor(c("dodgerblue1", "maroon2"), alpha.f = 0.5))
plot(expPcs[,1],expPcs[,2], col = as.factor(tissueSamplesInfo$sra.library_layout), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 1)




breakPoints <- seq(0.5,1,by = 0.05)
breakCols <- (adjustcolor(viridis(length(breakPoints), option = "inferno"), alpha.f = 0.5))


plot(expPcs[,1],expPcs[,2], col = breakCols[as.numeric(cut(tissueSamplesInfo$predictedTissueScore, breaks = breakPoints ))], pch = 16, cex = 1)
legend("bottomright",title="PredictionScore",legend=seq(0.5,1,by = 0.05),col = breakCols,pch=16)


plot(expPcs[,1],expPcs[,4], col = breakCols[as.numeric(cut(tissueSamplesInfo$predictedTissueScore, breaks = breakPoints ))], pch = 16, cex = 1)

legend("topleft",title="PredictionScore",legend=seq(0.5,1,by = 0.05),col = breakCols,pch=16)

pairs(expPcs[,1:10], col = breakCols[as.numeric(cut(tissueSamplesInfo$predictedTissueScore, breaks = breakPoints ))], cex = 1, upper.panel = NULL, pch = 16)



sum(expPcs[,2]>10)
x <- cbind(expPcs, tissueSamplesInfo)
View(x)
