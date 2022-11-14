#srun --cpus-per-task=20 --mem=100gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, sync = T)




remoter::client("localhost", port = 55508)#55506

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

table(samplesWithPrediction$predictedTissue)


sort(table(samplesWithPrediction$predictedTissue))
tissueClasses <- unique(samplesWithPrediction$predictedTissue)

#not used currently, we now use the expression data used for the primary QC and sample predictions.
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

load(file = "/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/Filtered_Matrices/TPM_log2_QNorm_QCed_CovCorrected_AllCovariates.RData", verbose = T)

  
mclapply(tissueClasses,  mc.cores = 10, function(tissue, exp){
  
  #load(file = paste0("perTissueNormalization/perTissueQq/",make.names(tissue),".RData"))
  
  tissueSamples <- rownames(samplesWithPrediction)[samplesWithPrediction$predictedTissue == tissue]
  tissueExp <- exp[,tissueSamples]
  
  
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
  
  
  explainedVariance <- eigenValues * 100 / nrow(expScale)
  
  
  pcaRes <- list(eigenVectors = eigenVectors, eigenValues = eigenValues, expPcs = expPcs, explainedVariance = explainedVariance)
  
  save(pcaRes, file = paste0("perTissueNormalization/perTissueQqPca/",make.names(tissue),".RData"))
  
  return(NULL)
  
}, exp = exp)

tissue = "Kidney"
tissue = "Brain-Nucleus accumbens (basal ganglia)"

ERP009290
ERP009290
samplesWithPrediction[samplesWithPrediction$study=="ERP009290",]

nonOutlierSampleList <- lapply(tissueClasses, function(tissue, samplesWithPrediction){
  
  load(file = paste0("perTissueNormalization/perTissueQqPca/",make.names(tissue),".RData"))

  expPcs <- pcaRes$expPcs[,1:10]
  explainedVariance <- pcaRes$explainedVariance
  tissueSamplesInfo <- samplesWithPrediction[rownames(expPcs),]
  studies <- length(unique(tissueSamplesInfo$study))
  
  #are in the same order
  write.table(cbind(tissueSamplesInfo, expPcs), col.names = T, row.names = F, sep = "\t", quote = F, file = paste0("perTissueNormalization/qcPlots/",make.names(tissue),".txt"))
  
  shortTissue <- ifelse(nchar(tissue) > 20, paste0(substr(tissue,0,17),"..."), tissue)
  
  breakPoints <- seq(0.5,1,by = 0.05)
  breakCols <- (adjustcolor(viridis(length(breakPoints), option = "inferno"), alpha.f = 0.5))
  
  expPcsMeans <- apply(expPcs, 2, mean)
  expPcsSds <- apply(expPcs, 2, sd)
  
  #Larger threshold is needed for wholeblood
  sdThreshold <- ifelse(tissue == "Whole Blood", 4,3)
  
  threshold <- expPcsMeans + sdThreshold * expPcsSds
  
  outlierPerComp <- sapply(1:10, function(i){
   
    abs(expPcs[,i]) > threshold[i]
  })
  tissueSamplesInfo$outlier <- apply(outlierPerComp, 1, any)
  sum(tissueSamplesInfo$outlier)
  
  ### Do some manual corrections
  if(tissue == "Airway basal cells"){
    tissueSamplesInfo$outlier[expPcs[,1] < 2] <- TRUE #Checked annotation, these are wrongly predicted
  }
  tissueSamplesInfo$outlier[tissueSamplesInfo$study == "ERP009290"] <- TRUE #Mixed tissue samples
      
  if(tissue == "Brain-Hindbrain-Fetal"){
    tissueSamplesInfo$outlier[expPcs[,1] < -10] <- TRUE #Checked annotation, these are wrongly predicted
    tissueSamplesInfo$outlier[expPcs[,3] > 3] <- TRUE #Checked annotation, these are wrongly predicted
  }
  if(tissue == "Brain-Nucleus accumbens (basal ganglia)"){
    tissueSamplesInfo$outlier[expPcs[,1] < -5] <- TRUE #Checked annotation, these are wrongly predicted
    tissueSamplesInfo$outlier[expPcs[,3] < -10] <- TRUE #Checked annotation, these are wrongly predicted
  }
  if(tissue == "Kidney"){
    tissueSamplesInfo$outlier[expPcs[,3] < -6] <- TRUE #Checked annotation, these are wrongly predicted
  }
  if(tissue == "Macrophages-iPSC"){
    tissueSamplesInfo$outlier[tissueSamplesInfo$study == "ERP020977"] <- FALSE #not real outlier, strange distribution due to stimulations.
  }
  if(tissue == "Monocytes"){
    tissueSamplesInfo$outlier[expPcs[,2] > 3] <- TRUE #Checked annotation, these are wrongly predicted
  }
  if(tissue == "Nasal Lavage"){
    tissueSamplesInfo$outlier[expPcs[,1] < -5] <- TRUE #Checked annotation, these are wrongly predicted
  }
  if(tissue == "Vagina"){
    tissueSamplesInfo$outlier[expPcs[,1] > 10] <- TRUE #Checked annotation, these are wrongly predicted
    tissueSamplesInfo$outlier[expPcs[,2] < -8] <- TRUE #Checked annotation, these are wrongly predicted
  }
  
  

  colnames(expPcs) <- paste0("Comp ",1:10, " (", round(explainedVariance[1:10],2) ,"%)")
  write.table(cbind(tissueSamplesInfo, expPcs), col.names = NA, row.names = T, sep = "\t", quote = F, file = paste0("perTissueNormalization/qcPlots/",make.names(tissue),".txt"))
  
  
  png(file = paste0("perTissueNormalization/qcPlots/",make.names(tissue),"1.png"), width = 1200, height = 900)
  #rpng(width = 1000, height = 1000)
  layout(matrix(c(1,1,1,1,2,3,4,8,5,6,7,8),ncol = 4, byrow = T), heights = c(0.1,1,1), widths = c(1,1,1,0.1))
  par(mar = c(0,0,0,0), xpd = NA, cex = 1.2)
  plot.new()
  plot.window(xlim = 0:1, ylim = 0:1)
  text(0.5,0.5, paste0(tissue, " (", nrow(tissueSamplesInfo) ,")"), cex = 2, font = 2)
  
  par(mar = c(4,4,3,0.5), xpd = NA)
  
  palette(adjustcolor(viridis(studies, option = "H"), alpha.f = 0.5))
  pchMap <- rep(c(15,16,17), length.out = studies)
  plot(expPcs[,1],expPcs[,2], col = as.factor(tissueSamplesInfo$study), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 1, main = paste0("Studies (", studies,")"), xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = paste0("Comp 2 (", round(explainedVariance[2],2) ,"%)"), bty = "n")
  
  
  palette(adjustcolor(c("lemonchiffon3", "darkorange1", "springgreen2"), alpha.f = 0.5))
  annotated <- factor(rep("Unkown", nrow(tissueSamplesInfo)), levels = c("Unkown", "Other", "Current"))
  annotated[!is.na(tissueSamplesInfo$annotatedTissue) & tissueSamplesInfo$annotatedTissue != tissue] <- "Other"
  annotated[!is.na(tissueSamplesInfo$annotatedTissue) & tissueSamplesInfo$annotatedTissue == tissue] <- "Current"
  plot(expPcs[,1],expPcs[,2], col = annotated, pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 1, main = paste0("Annotated as " , shortTissue), xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = paste0("Comp 2 (", round(explainedVariance[2],2) ,"%)"), bty = "n")
  
  
  palette(adjustcolor(c("dodgerblue1", "maroon2"), alpha.f = 0.5))
  plot(expPcs[,1],expPcs[,2], col = factor(tissueSamplesInfo$sra.library_layout, levels = c("paired", "single")), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 1, main = "Sequencing layout", xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = paste0("Comp 2 (", round(explainedVariance[2],2) ,"%)"), bty = "n")
  
  plot(expPcs[,1],expPcs[,2], col = breakCols[as.numeric(cut(tissueSamplesInfo$predictedTissueScore, breaks = breakPoints ))], pch = 16, cex = 1, main = "Prediction posterior probability", xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = paste0("Comp 2 (", round(explainedVariance[2],2) ,"%)"), bty = "n")
  
  plot(cumsum(explainedVariance)[1:10], bty = "n", pch = 16, xlab = "Components", ylab = "Cumulative explained variance (%)", main = "Explained variance", ylim = c(0,100), xlim = c(0,10))
  
  
  
  
  
  par(mar = c(0,2,3,1), xpd = NA)
  plot.new()
  plot.window(xlim = 0:1, ylim = 0:1)
  legend("topleft",title="Annotation",legend=c("Unkown", "Other", shortTissue), col = c("lemonchiffon3", "darkorange1", "springgreen2") , pch = 16, bty = "n")
  legend("top",title="Layout",legend=c("Single", "Paired"), col = c("maroon2", "dodgerblue1") , pch = 16, bty = "n")
  legend("topright",title="Probability",legend=seq(0.5,1,by = 0.05),col = breakCols,pch=16, bty = "n")
  
  
  dev.off()
  
  
  
  #png(file = paste0("perTissueNormalization/qcPlots/",make.names(tissue),"2.png"), width = 2000, height = 2000)
  #pairs(expPcs[,1:10], col = breakCols[as.numeric(cut(tissueSamplesInfo$predictedTissueScore, breaks = breakPoints ))], cex = 2, upper.panel = NULL, pch = 16)
  #dev.off()
  
  #png(file = paste0("perTissueNormalization/qcPlots/",make.names(tissue),"3.png"), width = 2000, height = 2000)
  #palette(adjustcolor(viridis(studies, option = "H"), alpha.f = 0.5))
  #pchMap <- rep(c(15,16,17), length.out = studies)
  #pairs(expPcs[,1:10], col = as.factor(tissueSamplesInfo$study), pch = pchMap[as.factor(tissueSamplesInfo$study)], cex = 2, upper.panel = NULL)
  #dev.off()
  
  #png(file = paste0("perTissueNormalization/qcPlots/",make.names(tissue),"3.png"), width = 2000, height = 2000)
  #palette(adjustcolor(c("grey", "firebrick3"), alpha.f = 0.5))
  #pairs(expPcs[,1:10], col = tissueSamplesInfo$outlier + 1, pch = 16, cex = 2, upper.panel = NULL)
  #dev.off()
  
  png(file = paste0("perTissueNormalization/qcPlots/",make.names(tissue),"4.png"), width = 1500, height = 700)
  #rpng(width = 1000, height = 1000)
  palette(adjustcolor(c("grey", "firebrick3"), alpha.f = 0.5))
  layout(matrix(c(1,1,1,1,1,2:11),ncol = 5, byrow = T), heights = c(0.1,1,1))
  par(mar = c(0,0,0,0), xpd = NA, cex = 1.2)
  plot.new()
  plot.window(xlim = 0:1, ylim = 0:1)
  text(0.5,0.5, paste0(tissue, " (", nrow(tissueSamplesInfo) ,")"), cex = 2, font = 2)
  
  par(mar = c(4,4,3,0.5), xpd = NA)
  
  for(i in 2:10){
    plot(expPcs[,1],expPcs[,i], col = tissueSamplesInfo$outlier + 1, pch = 16, cex = 1, main = paste0("Comp ", i), xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = paste0("Comp ", i," (", round(explainedVariance[i],2) ,"%)"), bty = "n")
    abline(v=c(-threshold[1],threshold[1]), lwd = 2, col = "firebrick3", xpd = FALSE)
    abline(h=c(-threshold[i],threshold[i]), lwd = 2, col = "firebrick3", xpd = FALSE)
  }
  
  par(mar = c(0,2,3,1), xpd = NA)
  plot.new()
  plot.window(xlim = 0:1, ylim = 0:1)
  legend("top",title="Outliers",legend=c("Included", "Excluded"), col = c("grey", "firebrick3") , pch = 16, bty = "n")
  
  
  
  dev.off()

  
  png(file = paste0("perTissueNormalization/qcPlots/",make.names(tissue),"5.png"), width = 1500, height = 700)
  #rpng()
  layout(matrix(c(1,1,1,1,1,2:11),ncol = 5, byrow = T), heights = c(0.1,1,1))
  par(mar = c(0,0,0,0), xpd = NA, cex = 1.2)
  plot.new()
  plot.window(xlim = 0:1, ylim = 0:1)
  text(0.5,0.5, paste0(tissue, " (", nrow(tissueSamplesInfo) ,")"), cex = 2, font = 2)
  
  par(mar = c(4,4,3,0.5), xpd = NA)
  
  for(i in 2:10){
    plot(expPcs[,1],expPcs[,i], col =  breakCols[as.numeric(cut(tissueSamplesInfo$predictedTissueScore, breaks = breakPoints ))], pch = 16, cex = 1, main = paste0("Comp ", i), xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = paste0("Comp ", i," (", round(explainedVariance[i],2) ,"%)"), bty = "n")
    abline(v=c(-threshold[1],threshold[1]), lwd = 2, col = "firebrick3", xpd = FALSE)
    abline(h=c(-threshold[i],threshold[i]), lwd = 2, col = "firebrick3", xpd = FALSE)
  }
  
  par(mar = c(0,2,3,1), xpd = NA)
  plot.new()
  plot.window(xlim = 0:1, ylim = 0:1)
  legend("top",title="Probability",legend=seq(0.5,1,by = 0.05),col = breakCols,pch=16, bty = "n")
  
  
  dev.off()
  
  
  
  return(tissueSamplesInfo[!tissueSamplesInfo$outlier,1:(ncol(tissueSamplesInfo)-1)])
  
}, samplesWithPrediction = samplesWithPrediction)

samplesWithPredictionNoOutliers <- do.call(rbind, nonOutlierSampleList)
#save(samplesWithPredictionNoOutliers, file = "tissuePredictions/samplesWithPrediction_16_09_22_noOutliers.RData", verbose = T)

load("tissuePredictions/samplesWithPrediction_16_09_22_noOutliers.RData", verbose = T)

tissue = "fibroblasts_cell-lines_smooth-muscle-cell_mesenchymal-stem-cells"
tissue = "HUVEC"

sink <- lapply(tissueClasses, function(tissue, exp){
  
  #load(file = paste0("perTissueNormalization/perTissueQq/",make.names(tissue),".RData"))
  
  tissueSamples <- rownames(samplesWithPredictionNoOutliers)[samplesWithPredictionNoOutliers$predictedTissue == tissue]
  tissueExp <- exp[,tissueSamples]
  
  save(tissueExp, file = paste0("perTissueNormalization/globalQqnorm/",make.names(tissue),".RData")) 
  
}, exp = exp)

#Create txt matrices
sink <- lapply(tissueClasses, function(tissue){
  
  load(file = paste0("perTissueNormalization/globalQqnorm/",make.names(tissue),".RData"), verbose = T) 
  
  write.table(tissueExp, file = gzfile(paste0("perTissueNormalization/qqExp/",make.names(tissue),".txt.gz")), sep = "\t", quote = F, col.names = NA)
  
})


cat(make.names(tissueClasses), sep = "\n")

sink <- lapply(tissueClasses, function(tissue){
  
  #load(file = paste0("perTissueNormalization/perTissueQq/",make.names(tissue),".RData"))
  
  #tissueSamples <- rownames(samplesWithPredictionNoOutliers)[samplesWithPredictionNoOutliers$predictedTissue == tissue]
  #tissueExp <- exp[,tissueSamples]
  
  load(file = paste0("perTissueNormalization/globalQqnorm/",make.names(tissue),".RData"), verbose = T) 
  
  write.table(tissueExp, file = gzfile("huvec.txt.gz"), sep = "\t", quote = F, col.names = NA)
  
  
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
  
  
  explainedVariance <- eigenValues * 100 / nrow(expScale)
  
  
  pcaRes <- list(eigenVectors = eigenVectors, eigenValues = eigenValues, expPcs = expPcs, explainedVariance = explainedVariance)
      
  save(pcaRes, file = paste0("perTissueNormalization/perTissueQqPcaNoOutliers/",make.names(tissue),".RData"))
  
  return(expPcs)
  
})

write.table(pcaRes$eigenVectors, file = "huvecEigenVectors.txt", sep = "\t", quote = F, col.names = NA)
write.table(pcaRes$eigenValues, file = "huvecEigenValues.txt", sep = "\t", quote = F, col.names = NA)
write.table(pcaRes$expPcs, file = "huvecPcs.txt", sep = "\t", quote = F, col.names = NA)

pcsPerTissue <- lapply(tissueClasses, function(tissue){
  load(file = paste0("perTissueNormalization/perTissueQqPcaNoOutliers/",make.names(tissue),".RData"))
  eigenvectors <- pcaRes$eigenVectors
  colnames(eigenvectors) <- paste0(tissue,"_",colnames(eigenvectors))
  
  return(eigenvectors)
})
#str(pcsPerTissue)
pcsPerTissue2 <- do.call(cbind, pcsPerTissue)

str(pcsPerTissue2)

rownames(pcsPerTissue2) <- (gsub("\\..+", "", rownames(pcsPerTissue2)))
write.table(pcsPerTissue2, file = "perTissueNormalization/perTissueQqPcaNoOutliers/combinedComponents.txt", sep = "\t", quote = FALSE, col.names = NA)

pcsPerTissue2t <- t(pcsPerTissue2)
pcsPerTissue2Scale = pcsPerTissue2t - rowMeans(pcsPerTissue2t)
# Standardize each variable
pcsPerTissue2Scale = pcsPerTissue2Scale / sqrt(rowSums(pcsPerTissue2Scale^2))

pcCorMatrix <- pcsPerTissue2Scale %*% t(pcsPerTissue2Scale)

range(pcCorMatrix)
range(diag(pcCorMatrix))

sum(pcCorMatrix[lower.tri(pcCorMatrix)] >= 0.8)

identicalPerPc <- apply(pcCorMatrix, 2, function(x){sum(x>=0.7)})
tail(sort(identicalPerPc))

hist(pcCorMatrix[,"Brain-Cortex_PC_3"])
dev.off()

pcCorMatrix[,"Whole Blood Fetal_PC_1"][pcCorMatrix[,"Whole Blood Fetal_PC_1"] >= 0.7]


compEigen <- eigen(pcCorMatrix)
str(compEigen)
sum(compEigen$values)

(numberOfCompsEigenvalue1 <- sum(as.numeric(compEigen$values) >= 1))

str(compEigen)
pcsOfComps <- t(pcsPerTissue2t) %*% compEigen$vectors[,1:numberOfCompsEigenvalue1]
colnames(pcsOfComps) <- paste0("PC_",1:ncol(pcsOfComps))
rownames(pcsOfComps) <- (gsub("\\..+", "", rownames(pcsOfComps)))
write.table(pcsOfComps, col.names = NA, sep = "\t", quote = F, file = gzfile("perTissueNormalization/perTissueQqPcaNoOutliers/pcaCombinedComponents.txt.gz"))
str(pcsOfComps)

rpng()
plot(cumsum(as.numeric(compEigen$values) * 100 / sum(as.numeric(compEigen$values))))
dev.off()

rpng()
plot(as.numeric(compEigen$values))
dev.off()
















head(as.numeric(compEigen$values))

sum(eigenValues >= 1)

compSvd <- svd(t(pcsPerTissue2Scale))
str(compSvd)

(numberOfCompsEigenvalue1 <- sum(compSvd$d^2>=1))

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

load("problem.RData")
load("solution.RData")


combinedCompsPcs <- compSvd$u[,1:numberOfCompsEigenvalue1] %*% diag(compSvd$d[1:numberOfCompsEigenvalue1])
combinedCompsPcs2 <- compSvd$u[,1:numberOfCompsEigenvalue1] %*% diag(compSvd$d[1:numberOfCompsEigenvalue1])
combinedCompsPcs3 <- compSvd$u[,1:numberOfCompsEigenvalue1] %*% diag(compSvd$d[1:numberOfCompsEigenvalue1])
combinedCompsPcs4 <- compSvd$u[,1:numberOfCompsEigenvalue1] %*% diag(compSvd$d[1:numberOfCompsEigenvalue1])
combinedCompsPcs5 <- compSvd$u[,1:numberOfCompsEigenvalue1] %*% diag(compSvd$d[1:numberOfCompsEigenvalue1])

cor.test(combinedCompsPcs[,1], combinedCompsPcsCluser[,1])

cor.test(combinedCompsPcs[,1], combinedCompsPcs2[,1])
cor.test(combinedCompsPcs[,1], combinedCompsPcs3[,1])
cor.test(combinedCompsPcs[,1], combinedCompsPcs4[,1])  
cor.test(combinedCompsPcs[,1], combinedCompsPcs5[,1])  

range(abs(combinedCompsPcs[,1]) - abs(combinedCompsPcs2[,1]))

plot(combinedCompsPcs[,1], combinedCompsPcs2[,1])
plot(combinedCompsPcs[,1], combinedCompsPcs3[,1])
plot(combinedCompsPcs[,1], combinedCompsPcs4[,1])
plot(combinedCompsPcs[,1], combinedCompsPcs5[,1])
dev.off()

save(compSvd, numberOfCompsEigenvalue1, file = "problem.RData")

str(combinedCompsPcs)
range(combinedCompsPcs)

plot(as.numeric(teest[,1]),combinedCompsPcs[,1])
dev.off()

  
head(compSvd$d^2)
head(compSvd$d^2)

eigenValues <- compSvd$d^2

rpng()
plot(cumsum(eigenValues * 100 / sum(eigenValues)))
dev.off()




rpng()
plot(eigenValues)
dev.off()

