#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)


remoter::client("localhost", port = 55556)#55556  55507

library(havok)



setwd("/groups/umcg-fg/tmp02/projects/genenetwork/recount3/")


load(file = "Metadata/combinedMeta_2022_09_15.RData", verbose = T)

combinedMeta <- combinedMeta[!combinedMeta$exclude,]

samplesPassInitialQc <- read.delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/Filtered_Matrices/TPM_log2_QNorm_QCed_CovCorrected_AllCovariates_samples.txt.gz", header = F)$V1

combinedMeta <- combinedMeta[rownames(combinedMeta) %in% samplesPassInitialQc,]


tcgaAll <- combinedMeta[combinedMeta$study == "TCGA",]

tcgaCancer <- tcgaAll[tcgaAll$Cancer,]

table(tcgaCancer$tcga.cgc_sample_sample_type)
sum(table(tcgaCancer$tcga.cgc_sample_sample_type))

tcgaCancer <- tcgaCancer[tcgaCancer$tcga.cgc_sample_sample_type == "Primary Tumor",]
table(tcgaCancer$tcga.gdc_cases.project.primary_site,useNA= "always")

rpng(width = 1000, height = 1000)
par(mar = c(9,4,3,1))
barplot(table(tcgaCancer$tcga.gdc_cases.project.primary_site), las = 2, main = "TCGA primary tumors")
dev.off()

dim(tcgaCancer)

rownames(tcgaCancer)




allFiles <- c("rse-tcga/rseTCGA.rda", "rse-tcga/rse_ESCA_TCGA.rda")


perChunkExp <- sapply(allFiles, function(file){
  
  loadedObject <- load(file)
  
  sreObjects <- get(loadedObject[1])
  
  #sometimes single RSE is not in list. Put in list of one to make code uniform
  if(!is.list(sreObjects)){
    sreObjects <- list(sreObjects)
  }
  
  #sreObject <- sreObjects[[1]]
  
  perStudyExp <- lapply(sreObjects, function(sreObject){
    studyExp <- sreObject@assays@data@listData$raw_counts
    return(studyExp[,colnames(studyExp) %in% rownames(tcgaCancer), drop = F])
  })
  
  return(do.call(cbind, perStudyExp))
  
})

tcgaExpRaw <- do.call(cbind, perChunkExp)
str(tcgaExpRaw)


tcgaCancer <- tcgaCancer[rownames(tcgaCancer) %in% colnames(tcgaExpRaw ),]

#use same order of samples
tcgaExpRaw <- tcgaExpRaw[,rownames(tcgaCancer) ]




tcgaExpRawNonZeroFraction <- apply(tcgaExpRaw, 1, function(x){
  return(sum(x > 0) / length(x))
})

tcgaExpRawNonZeroFraction[grepl("ENSG00000167751", names(tcgaExpRawNonZeroFraction))]
tcgaExpRawNonZeroFraction[grepl("ENSG00000142515", names(tcgaExpRawNonZeroFraction))]

#Expression in 50% of samples
tcgaExpRaw <- tcgaExpRaw[tcgaExpRawNonZeroFraction >= 0.5,]


library(DESeq2)
#limit expression to max int
tcgaExpRaw[tcgaExpRaw > .Machine$integer.max] <- .Machine$integer.max
tcgaExpRawDeSeq <- DESeqDataSetFromMatrix(tcgaExpRaw, tcgaCancer[,"tcga.gdc_cases.project.primary_site",drop = F], ~ tcga.gdc_cases.project.primary_site)

tcgaExpVst <- assay(vst(tcgaExpRawDeSeq, blind = F))


covariatesToCorrectFor <- read.delim("CovariateNames.txt", header = F)$V1
#TCGA is all paired so don't correct for that
covariatesToCorrectFor <- covariatesToCorrectFor[covariatesToCorrectFor != "sra.library_layout"]


#TCGA samples don't have sra.sample_spots instead use recount_qc.bc_frag.count
missingSampleSpots <- is.na(tcgaCancer[,"sra.sample_spots"])
tcgaCancer[missingSampleSpots,"sra.sample_spots"] <- tcgaCancer[missingSampleSpots,"recount_qc.bc_frag.count"]
tcgaCancerForCorrection <- tcgaCancer[,covariatesToCorrectFor]



tcgaExpVstGeneMean <- apply(tcgaExpVst, 1, mean)
rpng(width = 1000, height = 1000)
hist(tcgaExpVstGeneMean, main = "35523 genes expressed in 50% of the 8739 TCGA tumor samples", xlab = "Mean expression after VST (no correction)", breaks = 50)
dev.off()

cl <- makeCluster(20)  
tcgaExpVstCovCor <- parApply(cl, tcgaExpVst, 1 ,function(geneExp, tcgaCancerForCorrection){
  return(residuals(lm(geneExp ~ . ,data = tcgaCancerForCorrection)))
}, tcgaCancerForCorrection = tcgaCancerForCorrection)
stopCluster(cl)

tcgaExpVstCovCor <- t(tcgaExpVstCovCor)
tcgaExpVstCovCor2 <- tcgaExpVstCovCor + tcgaExpVstGeneMean

str(tcgaExpVstCovCor2)

#write.table(tcgaExpVstCovCor2, file = gzfile("tcga/tcgaCancerVstCovCor.txt.gz"), sep = "\t", quote = F, col.names = NA)
#write.table(tcgaCancer, file = gzfile("tcga/metaData.txt.gz"), sep = "\t", quote = F, col.names = NA)

#save(tcgaCancer, file = "tcga/metaData.RData")
#save(tcgaExpVstCovCor2, file = "tcga/tcgaCancerVstCovCor.RData")

load("tcga/metaData.RData")
load("tcga/tcgaCancerVstCovCor.RData", verbose =T)

expScale = tcgaExpVstCovCor2 - rowMeans(tcgaExpVstCovCor2);
#expScale = tcgaExpVst - rowMeans(tcgaExpVst);

# Standardize each variable
expScale = expScale / sqrt(rowSums(expScale^2));   

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

rpng(width = 1000, height = 1000)
plot(cumsum(explainedVariance), pch = 16, cex = 0.5, xlab = "component", ylab = "Cumulative explained %", main = "PCA on all TCGA samples")
abline(h = 80, lwd = 2, col = "darkred")
text(0,81,sum(cumsum(explainedVariance)<=80), adj = 0)
abline(h = 90, lwd = 2, col = "darkred")
text(0,91,sum(cumsum(explainedVariance)<=90), adj = 0)
dev.off()

(numberComponentsToInclude <- which.max(cumsum(explainedVariance) >= 85))




# library("havok")
# 
# medianSingularValue <- median(expSvd$d)
# 
# omega <- optimal_SVHT_coef(ncol(tcgaExpVstCovCor2) / nrow(tcgaExpVstCovCor2), sigma_known = F)
# threshold <- omega * medianSingularValue
# (numberComponentsToInclude <- sum(expSvd$d > threshold ))



write.table(eigenVectors[,1:numberComponentsToInclude], file = gzfile("tcga/tcgaCancerVstCovCorEigen.txt.gz"), sep = "\t", quote = F, col.names = NA)


library(viridisLite, lib.loc = .libPaths()[2])

rpng(width = 1000, height = 1000)
palette(adjustcolor(viridis(length(unique(tcgaCancer[,"tcga.gdc_cases.project.primary_site"])), option = "H"), alpha.f = 0.5))
plot(expPcs[,1],expPcs[,2], col = as.factor(tcgaCancer[,"tcga.gdc_cases.project.primary_site"]), pch = 16, cex = 1, xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = paste0("Comp 2 (", round(explainedVariance[2],2) ,"%)"), bty = "n")
legend("topright", fill = 1:length(levels(as.factor(tcgaCancer[,"tcga.gdc_cases.project.primary_site"]))), legend = levels(as.factor(tcgaCancer[,"tcga.gdc_cases.project.primary_site"])))
dev.off()

rpng(width = 1000, height = 1000)
palette(adjustcolor(viridis(length(unique(tcgaCancer[,"tcga.gdc_cases.project.primary_site"])), option = "H"), alpha.f = 0.5))
plot(expPcs[,1],expPcs[,3], col = as.factor(tcgaCancer[,"tcga.gdc_cases.project.primary_site"]), pch = 16, cex = 1, xlab = paste0("Comp 1 (", round(explainedVariance[1],2) ,"%)"), ylab = paste0("Comp 3 (", round(explainedVariance[3],2) ,"%)"), bty = "n")
legend("topright", fill = 1:length(levels(as.factor(tcgaCancer[,"tcga.gdc_cases.project.primary_site"]))), legend = levels(as.factor(tcgaCancer[,"tcga.gdc_cases.project.primary_site"])))
dev.off()




#Per tissue networks
tissues <- c("Prostate", "Breast", "Skin", "Colorectal", "Stomach", "Ovary")

tissue <- tissues[1]


sapply(tissues, function(tissue){
  
  length(rownames(tcgaCancer)[tcgaCancer$tcga.gdc_cases.project.primary_site == tissue])
})




numberOfComps <- lapply(tissues, function(tissue){
  
  tissueSamples <- rownames(tcgaCancer)[tcgaCancer$tcga.gdc_cases.project.primary_site == tissue]
  
  tissueExp <- tcgaExpVstCovCor2[,tissueSamples]
  
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
  
  (numberComponentsToInclude <- which.max(cumsum(explainedVariance) >= 80))
  
  
  
  # 
  # medianSingularValue <- median(expSvd$d)
  # 
  # omega <- optimal_SVHT_coef(ncol(expScale) / nrow(expScale), sigma_known = F)
  # threshold <- omega * medianSingularValue
  # numberComponentsToInclude <- sum(expSvd$d > threshold )
  # 
  # cat(paste0(tissue," ",numberComponentsToInclude) , "\n")
  h <- cumsum(explainedVariance)[numberComponentsToInclude ]
  
  #rpng(width = 1000, height = 1000)
  png(paste0("tcga/tissuePca/plots/",make.names(tissue),"_explainedVar.png"),width = 1000, height = 1000)
  plot(cumsum(explainedVariance), pch = 16, cex = 0.5, xlab = "Component", ylab = "Cumulative explained %", main = paste0("TCGA ", tissue))
  abline(h = h, lwd = 2, col = "darkred")
  text(0,h+1,numberComponentsToInclude, adj = 0)
  dev.off()
  
  
  write.table(eigenVectors[,1:numberComponentsToInclude], file = gzfile(paste0("tcga/tissuePca/",make.names(tissue),"_eigenVec.txt.gz")), sep = "\t", quote = F, col.names = NA)
  
  
  tissuePca <- list(eigenVectors = eigenVectors, eigenValues = eigenValues, expPcs = expPcs, explainedVariance = explainedVariance)
  
  save(tissuePca, file = paste0("tcga/tissuePca/",make.names(tissue),".RData"))
  return(numberComponentsToInclude)
})

tissue <- tissues[1]
numberOfComps <- lapply(tissues, function(tissue){
  
  load(file = paste0("tcga/tissuePca/",make.names(tissue),".RData"))
  
  
  str(tissuePca)
  
  numberSamples <- ncol(tissuePca$eigenVectors)
  
  explainedVariance <- tissuePca$eigenValues * 100 / sum( tissuePca$eigenValues)
  
  
  numberComponentsToIncludeVariance <- which.max(cumsum(explainedVariance) >= 80)
  
  sampleEigen <- explainedVariance * numberSamples / 100
  
  
  numberComponentsToIncludeSampleEigen <- sum( sampleEigen >= 1)
  
  
  numberComponentsToInclude <- min(numberComponentsToIncludeVariance, numberComponentsToIncludeSampleEigen)
  
   write.table(tissuePca$eigenVectors[,1:numberComponentsToInclude], file = gzfile(paste0("tcga/tissuePca/",make.names(tissue),"_eigenVec2.txt.gz")), sep = "\t", quote = F, col.names = NA)
  
  
  return(numberComponentsToInclude)
  

  
})
