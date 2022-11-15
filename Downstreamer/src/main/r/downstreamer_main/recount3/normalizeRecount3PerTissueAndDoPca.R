#srun --cpus-per-task=20 --mem=200gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, sync = T)




remoter::client("localhost", port = 55506)

library(DESeq2)
library(parallel)
library(viridisLite, lib.loc = .libPaths()[2])

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Recount3\\")

load(file = "perTissueNormalization/selectedSamplesRawExpression.RData", verbose = T)
load("tissuePredictions/samplesWithPrediction_16_09_22_noOutliers.RData", verbose = T)


sort(table(samplesWithPredictionNoOutliers$predictedTissue))
tissueClasses <- unique(samplesWithPredictionNoOutliers$predictedTissue)


#tissueClasses <- tissueClasses[1:29]
#tissueClasses <- tissueClasses[30:57]

#tissueClasses <- tissueClasses[c(1,2,6,14,55)]

#limit expression to max int
#selectedSamplesExp[selectedSamplesExp > .Machine$integer.max] <- .Machine$integer.max
tissue = "Kidney"

lapply(tissueClasses, function(tissue){
  
  tissueSamples <- rownames(samplesWithPredictionNoOutliers)[samplesWithPredictionNoOutliers$predictedTissue == tissue]
  tissueExp <- selectedSamplesExp[,tissueSamples]
  numberOfSamples <- length(tissueSamples)
  
  includedGenes <- apply(tissueExp, 1, function(x){(sum(x==0)/numberOfSamples) <= 0.5})
  
  tissueExp <- tissueExp[includedGenes,]
  
  print(paste(tissue,numberOfSamples,sum(includedGenes)))
  
  save(tissueExp, file = paste0("perTissueNormalization/raw/",make.names(tissue),".RData"))
  
})


perTissueExp <- lapply(tissueClasses, function(tissue){
  
  load(file = paste0("perTissueNormalization/raw/",make.names(tissue),".RData"))
  
  #limit expression to max int
  tissueExp[tissueExp > .Machine$integer.max] <- .Machine$integer.max

  vstExp <- vst(tissueExp, nsub = 2000)
  
  save(vstExp, file = paste0("perTissueNormalization/vstExp/",make.names(tissue),".RData"))
  return(NULL)
  
})


#names(perTissueExp) <- tissueClasses
#save(perTissueExp, file = "perTissueNormalization/selectedSamplesRawExpressionPerTissue.RData")
tissue = "Kidney"
load(file = "Metadata/combinedMeta_2022_09_15.RData", verbose = T)

covariatesToCorrectFor <- read.delim("CovariateNames.txt", header = F)$V1

str(sum(is.na(combinedMeta[,"sra.sample_spots"])))

#TCGA samples don't have sra.sample_spots instead use recount_qc.bc_frag.count
missingSampleSpots <- is.na(combinedMeta[,"sra.sample_spots"])
combinedMeta[missingSampleSpots,"sra.sample_spots"] <- combinedMeta[missingSampleSpots,"recount_qc.bc_frag.count"]

#TCGA and Gtex don't report layout but is all paired
missingLayout <- is.na(combinedMeta[,"sra.library_layout"])
combinedMeta[missingLayout,"sra.library_layout"] <- "paired"

combinedMetaSelection <- combinedMeta[,covariatesToCorrectFor]
combinedMetaSelection$sra.library_layout <- as.factor(combinedMetaSelection$sra.library_layout)
table(combinedMeta[,"sra.library_layout"],useNA = "a")

cl <- makeCluster(20)

tissue <- "Salivary Gland-Minor Salivary Gland"

sink <- lapply(tissueClasses, function(tissue){
  cat(tissue, "\n")
  
  
  load(file = paste0("perTissueNormalization/vstExp/",make.names(tissue),".RData"))

  combinedMetaTissue <- combinedMetaSelection[colnames(vstExp),]
  if(min(table(combinedMetaTissue$sra.library_layout)) < 10){
    combinedMetaTissue <- combinedMetaTissue[,-match("sra.library_layout",colnames(combinedMetaTissue))]
  }
  
  vstExpCovCor <- parApply(cl, vstExp, 1 ,function(geneExp, combinedMetaTissue){
    return(residuals(lm(geneExp ~ . ,data = combinedMetaTissue)))
  }, combinedMetaTissue = combinedMetaTissue)
  vstExpCovCor <- t(vstExpCovCor)
  
  
  save(vstExpCovCor, file = paste0("perTissueNormalization/vstExpCovCor/",make.names(tissue),".RData"))

})

stopCluster(cl)


sink <- lapply(tissueClasses, function(tissue){
  cat(tissue, "\n")

  
  load(file = paste0("perTissueNormalization/vstExpCovCor/",make.names(tissue),".RData"))
  
  #https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
  expScale = vstExpCovCor - rowMeans(vstExpCovCor);
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
  
  nrSamples <- ncol(expScale)
  
  expSvd <- svd(expScale, nu = nrSamples, nv = min(nrSamples, 50))
  
  
  eigenValues <- expSvd$d^2
  eigenVectors <- expSvd$u
  colnames(eigenVectors) <- paste0("Comp_",1:ncol(eigenVectors))
  rownames(eigenVectors) <- rownames(expScale)
  
  expPcs <- expSvd$v %*% diag(expSvd$d[1:ncol(expSvd$v)])
  colnames(expPcs) <- paste0("Comp_",1:ncol(expPcs))
  rownames(expPcs) <- colnames(expScale)
  
  explainedVariance <- eigenValues * 100 / sum(eigenValues)
  
  
  tissueVstPca <- list(eigenVectors = eigenVectors, eigenValues = eigenValues, expPcs = expPcs, explainedVariance = explainedVariance)
  
  save(tissueVstPca, file = paste0("perTissueNormalization/vstCovCorPca/",make.names(tissue),".RData"))
  
})


sink <- lapply(tissueClasses, function(tissue){
  cat(tissue, "\n")
  
  load(file = paste0("perTissueNormalization/vstExpCovCor/",make.names(tissue),".RData"))
  
  write.table(vstExpCovCor, file = gzfile(paste0("perTissueNormalization/dataForLude/",make.names(tissue),"_vstCovCor.txt.gz")))
  

})


sink <- lapply(tissueClasses, function(tissue){
  cat(tissue, "\n")
  
  load(file = paste0("perTissueNormalization/vstCovCorPca/",make.names(tissue),".RData"))
  
  write.table(tissueVstPca$eigenVectors, file = gzfile(paste0("perTissueNormalization/dataForLude/",make.names(tissue),"_eigenVec.txt.gz")))
  
  
})


str(vstExp)
rownames(vstExp) <- (gsub("\\..+", "", rownames(vstExp)))
write.table(vstExp, file = gzfile(paste0("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/perTissueNormalization/vstExp/",make.names(tissue),".txt.gz")), sep = "\t", col.names = NA, quote = F)

eigenVec <- tissueVstPca$eigenVectors
rownames(vstExp) <- (gsub("\\..+", "", rownames(eigenVec)))
write.table(eigenVec, file = gzfile(paste0("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/perTissueNormalization/vstPca/",make.names(tissue),"_eigenVec.txt.gz")), sep = "\t", col.names = NA, quote = F)


metaNumeric <- as.matrix(combinedMeta[,sapply(combinedMeta, is.numeric)])



marCorZPerTissue <- lapply(tissueClasses, function(tissue){

  load(file = paste0("perTissueNormalization/vstPca/",make.names(tissue),".RData"), verbose = T)
  
  pcs <- tissueVstPca$expPcs
  metaTest <- metaNumeric[rownames(pcs),]

  metaVPcsZ <- apply(pcs, 2, function(x){
    apply(metaTest, 2, function(y, x){
      z <- !is.na(x) & !is.na(y)
      if(sum(z) < 10){
        return(0)
      }
      p <- cor.test(x[z],y[z])$p.value
      return(-qnorm(p/2))
    }, x = x)
  })
  
  maxCorZ <- apply(abs(metaVPcsZ), 1, max, na.rm = T)
  maxCorZ[is.infinite(maxCorZ)] <- 0

  return(maxCorZ)

})

str(marCorZPerTissue)

marCorZPerTissue2 <- do.call(cbind, marCorZPerTissue)
str(marCorZPerTissue2)

covariateZscores <- apply(marCorZPerTissue2,1,mean)

marCorZPerTissue2["recount_qc.star.%_of_reads_mapped_to_multiple_loci",]
library(vioplot)
vioplot(marCorZPerTissue2["recount_qc.star.%_of_reads_mapped_to_multiple_loci",])

library(beeswarm)
beeswarm(marCorZPerTissue2["recount_qc.star.%_of_reads_mapped_to_multiple_loci",],
         method = 'swarm',
         pch = 16, ylim = c(0,30), ylab = "z-score", main = "Highest correlation between components and\n % non uniquely mapping reads.")

dev.off()

cor.test(pcs[,1], metaTest[,"recount_qc.star.%_of_reads_mapped_to_multiple_loci"])
plot(pcs[,1], metaTest[,"recount_qc.star.%_of_reads_mapped_to_multiple_loci"], pch = 16, col=adjustcolor("grey", alpha.f = 0.5))
dev.off()

pcs <- expPcs


tissueSamplesInfo <- samplesWithPredictionNoOutliers[rownames(pcs),]
studies <- length(unique(tissueSamplesInfo$study))

library(viridisLite, lib.loc = .libPaths()[2])

rpng(width = 1500, height = 1500)
palette(adjustcolor(viridis(studies, option = "H"), alpha.f = 0.5))
pchMap <- rep(c(15,16,17), length.out = studies)
plot(pcs[,1],metaTest[,"recount_qc.star.%_of_reads_mapped_to_multiple_loci"], col = as.factor(tissueSamplesInfo$study), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 2, main = paste0("Studies (", studies,")"), xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = "Percentage read map multiple loci", bty = "n")
plot(pcs[,1],pcs[,2], col = as.factor(tissueSamplesInfo$study), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 2, main = paste0("Studies (", studies,")"), xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = "Percentage read map multiple loci", bty = "n")
dev.off()






hist(metaVPcsZ, breaks = 100)
dev.off()
library(heatmap3)
rpng(width = 1000, height = 1000)
heatmap3(metaVPcsZ, balanceColor = T, scale = "none", Rowv = NA, Colv =NA)
dev.off()




#OLD



#save(perTissueExp, perTissuePca, file = "perTissueNormalization/tmpTestSession.RData")
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

legend("topleft",title="PredictionScore",legend=seq(0.5,1,by = 0.05),col = breakCols,pch=16)

pairs(expPcs[,1:10], col = breakCols[as.numeric(cut(tissueSamplesInfo$predictedTissueScore, breaks = breakPoints ))], cex = 1, upper.panel = NULL, pch = 16)



  sum(expPcs[,2]>10)
x <- cbind(expPcs, tissueSamplesInfo)
View(x)
