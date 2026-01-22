#srun --cpus-per-task=20 --mem=200gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, sync = T)




remoter::client("localhost", port = 55001)#55556 

library(DESeq2)
library(parallel)
library(viridisLite, lib.loc = .libPaths()[2])
library("havok")

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Recount3\\")

#load(file = "perTissueNormalization/selectedSamplesRawExpression.RData", verbose = T)
load("tissuePredictions/samplesWithPrediction_16_09_22_noOutliers.RData", verbose = T)
str(samplesWithPredictionNoOutliers)

#write.table(samplesWithPredictionNoOutliers, file = "tissuePredictions/samplesWithPrediction_16_09_22_noOutliers.txt", quote = F, col.names = NA, sep = "\t")


sort(table(samplesWithPredictionNoOutliers$predictedTissue))
tissueClasses <- unique(samplesWithPredictionNoOutliers$predictedTissue)
sort(tissueClasses)

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

tissue <- "Microglia"

sink <- lapply(tissueClasses, function(tissue){
  cat(tissue, "\n")
  
  load(file = paste0("perTissueNormalization/vstExp/",make.names(tissue),".RData"))

  combinedMetaTissue <- combinedMetaSelection[colnames(vstExp),]
  if(min(table(combinedMetaTissue$sra.library_layout)) < 10){
    combinedMetaTissue <- combinedMetaTissue[,-match("sra.library_layout",colnames(combinedMetaTissue))]
  }
  
  geneMean <- rowMeans(vstExp)
  
  vstExpCovCor <- parApply(cl, vstExp, 1 ,function(geneExp, combinedMetaTissue){
    return(residuals(lm(geneExp ~ . ,data = combinedMetaTissue)))
  }, combinedMetaTissue = combinedMetaTissue)
  vstExpCovCor <- t(vstExpCovCor)
  
  vstExpCovCor <- vstExpCovCor + geneMean
  
  save(vstExpCovCor, file = paste0("perTissueNormalization/vstExpCovCor/",make.names(tissue),".RData"))
  write.table(vstExpCovCor, file = gzfile(paste0("perTissueNormalization/dataForLude/",make.names(tissue),"_vstCovCor.txt.gz")), quote = F, sep = "\t", col.names = NA)
  

})

stopCluster(cl)



cor.test(vstExpCovCor[,1], vstExp[,1])
cor.test(vstExpCovCorOld[,1], vstExp[,1])

cor.test(vstExpCovCor[,1], vstExpCovCorOld[,1])

plot(vstExpCovCor[,1], vstExpCovCorOld[,1])
dev.off()

#vstExpCovCorOld <- vstExpCovCor
sort(tissueClasses)
tissue <- "Prostate"

numberOfComps <- lapply(tissueClasses, function(tissue){

  
  
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
  
  
  numberGenes <- nrow(expScale)
  
  
  numberComponentsToInclude <- which.max(cumsum(explainedVariance) >= 80)
  #cumsum(explainedVariance)[numberComponentsToInclude]
    
  #   
  # cronbachAlpha <- lapply(as.list(1:min(2000, nrSamples)), function(comp,numberGenes, expScale,eigenVectors){
  #   
  #   geneVariance <- sapply(1:nrow(expScale), function(r){
  #     var(expScale[r,] * eigenVectors[r,comp])
  #   })
  #   return((numberGenes / (numberGenes - 1))*(1 - (sum(geneVariance) / var(t(expScale) %*% eigenVectors[,comp]))))
  # }, numberGenes = numberGenes, expScale = expScale,eigenVectors = eigenVectors)
  # #,  mc.cores = 20
  # cronbachAlpha <- unlist(cronbachAlpha)
  # 
  # cronbachAlpha <- cronbachAlpha[cronbachAlpha>0]
  # 
  # (numberComponentsToInclude <- min(which(cronbachAlpha < 0.7))-1)
  
  
  #medianSingularValue <- median(expSvd$d)
  
  #omega <- optimal_SVHT_coef(ncol(expScale) / nrow(expScale), sigma_known = F)
  #threshold <- omega * medianSingularValue
  #numberComponentsToInclude <- sum(expSvd$d > threshold )
  
  cat(paste0(tissue," ",numberComponentsToInclude) , "\n")
  h <- cumsum(explainedVariance)[numberComponentsToInclude ]

  #rpng(width = 1000, height = 1000)
  png(paste0("perTissueNormalization/vstCovCorPca/plots/",make.names(tissue),"_explainedVar.png"),width = 1000, height = 1000)
  plot(cumsum(explainedVariance), pch = 16, cex = 0.5, xlab = "Component", ylab = "Cumulative explained %", main = tissue)
  abline(h = h, lwd = 2, col = "darkred")
  text(0,h+1,numberComponentsToInclude, adj = 0)
  dev.off()

  #rpng(width = 1000, height = 1000)
  png(paste0("perTissueNormalization/vstCovCorPca/plots/",make.names(tissue),"_eigenvalues.png"),width = 1000, height = 1000)
  plot(eigenValues, log = "y", ylab = "Eigenvalues")
  abline(v = numberComponentsToInclude, lwd = 2, col = "darkred")
  dev.off()
  
  
  write.table(eigenVectors[,1:numberComponentsToInclude], file = gzfile(paste0("perTissueNormalization/vstCovCorPca/",make.names(tissue),"_eigenVec.txt.gz")), sep = "\t", quote = F, col.names = NA)
  
  
  tissueVstPca <- list(eigenVectors = eigenVectors, eigenValues = eigenValues, expPcs = expPcs, explainedVariance = explainedVariance)
  
  save(tissueVstPca, file = paste0("perTissueNormalization/vstCovCorPca/",make.names(tissue),".RData"))
  return(numberComponentsToInclude)
})


tissue <- "Prostate"

numberOfComps <- lapply(tissueClasses, function(tissue){
  
  load( file = paste0("perTissueNormalization/vstCovCorPca/",make.names(tissue),".RData"))
  
  numberSamples <- ncol(tissueVstPca$eigenVectors)

  
  explainedVariance <- tissueVstPca$eigenValues * 100 / sum( tissueVstPca$eigenValues)
  
  
  numberComponentsToIncludeVariance <- which.max(cumsum(explainedVariance) >= 80)
  
  sampleEigen <- explainedVariance * numberSamples / 100
  
  
  numberComponentsToIncludeSampleEigen <- sum( sampleEigen >= 1)
  
  sampleEigen[numberComponentsToIncludeSampleEigen]
  
  numberComponentsToInclude <- min(numberComponentsToIncludeVariance, numberComponentsToIncludeSampleEigen)
  

  
  write.table(tissueVstPca$eigenVectors[,1:numberComponentsToInclude], file = gzfile(paste0("perTissueNormalization/vstCovCorPca/",make.names(tissue),"_eigenVec2.txt.gz")), sep = "\t", quote = F, col.names = NA)
  
  
  return(numberComponentsToInclude)

})

numberOfComps <- do.call("c", numberOfComps)
names(numberOfComps) <- tissueClasses





#as.data.frame(numberOfComps)

nrSamplesCombined <- 1#nrow(samplesWithPredictionNoOutliers)
nrComponentsCominbedNetwork <- 848

sampleCounts <- table(samplesWithPredictionNoOutliers$predictedTissue)

rpng(width = 1000, height = 1000)
plot(as.numeric(sampleCounts[names(numberOfComps)]), numberOfComps, 
     xlim = c(min(sampleCounts),nrSamplesCombined), 
     ylim = c(min(numberOfComps,nrComponentsCominbedNetwork),max(numberOfComps, nrComponentsCominbedNetwork)), 
     log = "xy", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.7) ,
     xlab = "Number of samples", ylab = "Number of components")
points(nrSamplesCombined, nrComponentsCominbedNetwork, pch = 16, col=adjustcolor("magenta1", alpha.f = 0.5))
dev.off()

write.table(as.data.frame(numberOfComps), file = "perTissueNormalization/vstCovCorPca/compsPerTissue.txt", col.names = NA, quote = F, sep = "\t")
write.table(as.data.frame(numberOfSamples), file = "perTissueNormalization/vstCovCorPca/samplesPerTissue.txt", col.names = NA, quote = F, sep = "\t")

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
plot(pcs[,1],pcs[,2], col = as.factor(tifssueSamplesInfo$study), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 2, main = paste0("Studies (", studies,")"), xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = "Percentage read map multiple loci", bty = "n")
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




tissue <- "Kidney"
sink <- lapply(tissueClasses, function(tissue){
  
  
  
  load(file = paste0("perTissueNormalization/vstExpCovCor/",make.names(tissue),".RData"))
  
  corMatrix <- cor(t(vstExpCovCor))
  write.table(corMatrix, file = gzfile(paste0("perTissueNormalization/dataForLude/",make.names(tissue),"_vstCovCor_coexp.txt.gz")), quote = F, sep = "\t", col.names = NA)
  
})