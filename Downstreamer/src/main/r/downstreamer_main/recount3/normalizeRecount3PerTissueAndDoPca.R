#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, sync = T)


remoter::client("localhost", port = 55501)

load(file = "perTissueNormalization/selectedSamplesRawExpression.RData")


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
table(is.wholenumber(selectedSamplesExp))


tissueClasses <- unique(samplesWithPrediction$predictedTissue)

tissue <- tissueClasses[1]

perTissueExp <- lapply(tissueClasses, function(tissue){
  tissueSamples <- rownames(samplesWithPrediction)[samplesWithPrediction$predictedTissue == tissue]
  tissueExp <- selectedSamplesExp[,tissueSamples]
  numberOfSamples <- length(tissueSamples)
  
  includedGenes <- apply(tissueExp, 1, function(x){(sum(x==0)/numberOfSamples) >= 0.5})
  table(includedGenes)
  
})
names(perTissueExp) <- tissueClasses
#save(perTissueExp, file = "perTissueNormalization/selectedSamplesRawExpressionPerTissue.RData")




library(DESeq2)

perTissueExpRlog <- lapply(tissueClasses, function(tissue){
  
  
  rawCounts <- perTissueExp[[tissue]]
  rlogExp <- rlog(rawCounts)
  
  x <- apply(log2(rawCounts), 1, mean)
  hist(x)
  dev.off()
  
  
  
})


perTissueExpRlog <- lapply(tissueClasses, function(tissue){
  
  
  
  
})


#https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
expScale = exp - rowMeans(exp);
# Standardize each variable
expScale = expScale / sqrt(rowSums(expScale^2));   
expCov = tcrossprod(expScale);#equevelent to correlation due to center scale
range(expCov)
str(expCov)

expEigen <- eigen(expCov)

eigenVectors <- expEigen$vectors
colnames(eigenVectors) <- paste0("PC_",1:ncol(eigenVectors))
rownames(eigenVectors) <- rownames(expScale)
str(eigenVectors)

eigenValues <- expEigen$values
names(eigenValues) <- paste0("PC_",1:length(eigenValues))
str(eigenValues)















