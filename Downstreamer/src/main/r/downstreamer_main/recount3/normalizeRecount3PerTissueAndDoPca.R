#srun --cpus-per-task=20 --mem=200gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, sync = T)


remoter::client("localhost", port = 55504)

library(DESeq2)
library(parallel)
library(viridisLite, lib.loc = .libPaths()[2])

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")


load(file = "perTissueNormalization/selectedSamplesRawExpression.RData", verbose = T)
load("tissuePredictions/samplesWithPrediction_16_09_22.RData", verbose = T)


tissueClasses <- unique(samplesWithPrediction$predictedTissue)

#tissueClasses <- tissueClasses[1:29]
#tissueClasses <- tissueClasses[30:57]

#tissueClasses <- tissueClasses[c(1,2,6,14,55)]

#limit expression to max int
selectedSamplesExp[selectedSamplesExp > .Machine$integer.max] <- .Machine$integer.max


mclapply(tissueClasses,  mc.cores = 20, function(tissue){
  
  tissueSamples <- rownames(samplesWithPrediction)[samplesWithPrediction$predictedTissue == tissue]
  tissueExp <- selectedSamplesExp[,tissueSamples]
  numberOfSamples <- length(tissueSamples)
  
  includedGenes <- apply(tissueExp, 1, function(x){(sum(x==0)/numberOfSamples) <= 0.5})
  
  tissueExp <- tissueExp[includedGenes,]
  
  mode(tissueExp) <- "integer"
  
  save(tissueExp, file = paste0("perTissueNormalization/raw/",make.names(tissue),".RData"))
  
})




#cl <- makeCluster(20)

#clusterExport(cl, c("samplesWithPrediction", "selectedSamplesExp"))
#tissue <- "Prostate"
#
perTissueExp <- mclapply(tissueClasses,  mc.cores = 20, function(tissue){
  
  load(file = paste0("perTissueNormalization/raw/",make.names(tissue),".RData"))
  rlogExp <- rlog(tissueExp)
  save(rlogExp, file = paste0("perTissueNormalization/rlogExp/",make.names(tissue),".RData"))
  return(NULL)
  
})

#names(perTissueExp) <- tissueClasses
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

save(perTissueExp, perTissuePca, file = "perTissueNormalization/tmpTestSession.RData")
#load(file = "perTissueNormalization/tmpTestRlog.RData")

save(expPcs, samplesWithPrediction, file = "perTissueNormalization/tmpTest2.RData")
load("perTissueNormalization/tmpTest2.RData")


samplesWithPrediction$sra.library_layout[samplesWithPrediction$study=="TCGA"] <- "paired"
samplesWithPrediction$sra.library_layout[samplesWithPrediction$study=="GTEx"] <- "paired"


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


plot(expPcs[,1],expPcs[,5], col = breakCols[as.numeric(cut(tissueSamplesInfo$predictedTissueScore, breaks = breakPoints ))], pch = 16, cex = 1)


pairs(expPcs[,1:10], col = breakCols[as.numeric(cut(tissueSamplesInfo$predictedTissueScore, breaks = breakPoints ))], cex = 1, upper.panel = NULL, pch = 16)



  sum(expPcs[,2]>10)
x <- cbind(expPcs, tissueSamplesInfo)
View(x)
