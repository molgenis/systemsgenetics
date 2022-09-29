#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, sync = T)


remoter::client("localhost", port = 55501)

library(DESeq2)


setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")


load(file = "perTissueNormalization/selectedSamplesRawExpression.RData", verbose = T)



tissueClasses <- unique(samplesWithPrediction$predictedTissue)

tissueClasses <- tissueClasses[c(1,2,6,14,55)]


perTissueExp <- lapply(tissueClasses, function(tissue){
  tissueSamples <- rownames(samplesWithPrediction)[samplesWithPrediction$predictedTissue == tissue]
  tissueExp <- selectedSamplesExp[,tissueSamples]
  numberOfSamples <- length(tissueSamples)
  
  includedGenes <- apply(tissueExp, 1, function(x){(sum(x==0)/numberOfSamples) <= 0.5})
  
  
  tissueExp <- tissueExp[includedGenes,]
  
  mode(tissueExp) <- "integer"
  
  rlogExp <- rlog(tissueExp)

  return(rlogExp)
  
})
names(perTissueExp) <- tissueClasses
#save(perTissueExp, file = "perTissueNormalization/selectedSamplesRawExpressionPerTissue.RData")


perTissuePca <- lapply(perTissueExp, function(exp){
  
  #https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
  expScale = exp - rowMeans(exp);
  # Standardize each variable
  expScale = expScale / sqrt(rowSums(expScale^2));   
  expCov = tcrossprod(expScale);#equevelent to correlation due to center scale
  
  expEigen <- eigen(expCov)
  
  eigenVectors <- expEigen$vectors
  colnames(eigenVectors) <- paste0("PC_",1:ncol(eigenVectors))
  rownames(eigenVectors) <- rownames(expScale)
  
  eigenValues <- expEigen$values
  names(eigenValues) <- paste0("PC_",1:length(eigenValues))
  
  #Here calculate sample principle components. Number needed is arbritary (no more than eigen vectors)
  expPcs <- t(expScale) %*% expEigen$vectors[,1:10]
  colnames(expPcs) <- paste0("PC_",1:ncol(expPcs))
  
  return(list(eigenVectors, eigenValues, expPcs))
  
})

#save(rlogExp, file = "perTissueNormalization/tmpTestRlog.RData")
#load(file = "perTissueNormalization/tmpTestRlog.RData")

exp <- rlogExp





tissueSamplesInfo <- samplesWithPrediction[rownames(expPcs),]
str(tissueSamplesInfo)

studies <- length(unique(tissueSamplesInfo$study))


library(viridisLite, lib.loc = .libPaths()[2])
palette(adjustcolor(viridis(studies, option = "H"), alpha.f = 0.5))

pchMap <- rep(c(15,17,19), length.out = studies)

rpng(width = 1000, height = 1000)
plot(expPcs[,1],expPcs[,2], col = as.factor(tissueSamplesInfo$study), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 2)
dev.off()
