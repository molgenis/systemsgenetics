library(parallel)

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

load(file = "combinedMeta_2022_08_14.RData", verbose = T)
load(file = "Recount3_QC_2ndRun/PCA_Patrick/pcs.RData", verbose = T)
load(file = "Recount3_QC_2ndRun/PCA_Patrick/eigen.RData")



pcsAndMeta <- merge(expPcs[,1:570], combinedMeta, by = 0, all.x = T)
rownames(pcsAndMeta) <- pcsAndMeta$Row.names
#save(pcsAndMeta, file = "DataForLasso.RData")

#write.table(merge(combinedMeta, expPcs[,1:100], by = 0, all.y = T), file = "Metadata/pcsAndAnnotations.txt", sep = "\t", quote = FALSE, col.names = NA)

dim(pcsAndMeta)
str(combinedMeta)

defaultCol <- adjustcolor("grey", alpha.f = 0.6)
tissueCol <- read.delim("Recount3_QC_2ndRun/SRA_Studies_Annotations_Patrick/Annotations_color2.txt", row.names = 1)



pcsAndMeta$col <- defaultCol
tissueAndCol <- pcsAndMeta[,"Tissue"] %in% tissueCol$PlotClass
pcsAndMeta$col[tissueAndCol] <-  adjustcolor(tissueCol$col[match(pcsAndMeta[tissueAndCol,"Tissue"], tissueCol$PlotClass)], alpha.f = 0.6)
tissue2AndCol <- pcsAndMeta[,"Tissue2"] %in% tissueCol$PlotClass
pcsAndMeta$col[tissue2AndCol] <-  adjustcolor(tissueCol$col[match(pcsAndMeta[tissue2AndCol,"Tissue2"], tissueCol$PlotClass)], alpha.f = 0.6)
plotOrderTissues <- order((pcsAndMeta$col != defaultCol) + 1)


rpng(width = 800, height = 800)
#pdf(file = "test.pdf")
plot(pcsAndMeta[plotOrderTissues,"PC_1"], pcsAndMeta[plotOrderTissues,"PC_2"], col = pcsAndMeta$col[plotOrderTissues], cex = 0.3, pch = 16)
dev.off()



defaultCol <- adjustcolor("grey", alpha.f = 0.3)
pcsAndMeta$colCelline <- defaultCol
pcsAndMeta$colCelline[!is.na(pcsAndMeta[,"Cellline"]) & pcsAndMeta[,"Cellline"]] <-  adjustcolor("magenta", alpha.f = 0.6)
pcsAndMeta$colCelline[!is.na(pcsAndMeta[,"Cellline"]) & !pcsAndMeta[,"Cellline"]] <-  adjustcolor("royalblue1", alpha.f = 0.6)
pcsAndMeta$colCelline[!is.na(pcsAndMeta[,"Cancer"]) & pcsAndMeta[,"Cancer"]] <-  adjustcolor("forestgreen", alpha.f = 0.6)

plotOrderCellline <- rep(0,nrow(pcsAndMeta))
plotOrderCellline[!is.na(pcsAndMeta[,"Cellline"]) & pcsAndMeta[,"Cellline"]] <- 2
plotOrderCellline[!is.na(pcsAndMeta[,"Cellline"]) & !pcsAndMeta[,"Cellline"]] <- 1
plotOrderCellline[!is.na(pcsAndMeta[,"Cancer"]) & pcsAndMeta[,"Cancer"]] <-  3
plotOrderCellline <- order(plotOrderCellline)

table(plotOrderCellline)

table(pcsAndMeta$Tissue2, pcsAndMeta$Cellline, useNA = "a")

  rpng(width = 800, height = 800)
#pdf(file = "test.pdf")
plot(pcsAndMeta[plotOrderCellline,"PC_1"], pcsAndMeta[plotOrderCellline,"PC_2"], col = pcsAndMeta$colCelline[plotOrderCellline], cex = 0.3, pch = 16)
#plot(pcsAndMeta[,"PC_1"], pcsAndMeta[,"PC_2"], col = pcsAndMeta$colCelline[], cex = 0.3, pch = 16)
dev.off()


table(pcsAndMeta$Cellline[!is.na(pcsAndMeta$CelllineName)&pcsAndMeta$CelllineName=="iPSC"], useNA = "always")


threshold <- c("PC_1" = 0.2, "PC_2" = -0.025, "PC_24" = -0.04, "PC_27" = 0.09, "PC_32" = -0.05, "PC_62" = 0.03, "PC_56" = -0.07, "PC_63" = 0.04)
rescueThreshold <- c("PC_10" = -0.15, "PC_18" = -0.1, "PC_20" = -0.125, "PC_21" = -0.18)

pcsAndMeta$excludeBasedOnPredictionCellline <- FALSE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_2 > threshold["PC_2"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_1 > threshold["PC_1"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_24 < threshold["PC_24"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_27 > threshold["PC_27"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_32 < threshold["PC_32"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_62 > threshold["PC_62"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_56 < threshold["PC_56"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_63 > threshold["PC_63"]] <- TRUE

pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_10 < rescueThreshold["PC_10"]] <- FALSE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_18 < rescueThreshold["PC_18"]] <- FALSE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_20 < rescueThreshold["PC_20"]] <- FALSE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_21 < rescueThreshold["PC_21"]] <- FALSE

pcsAndMeta$colExclude <- adjustcolor("goldenrod2", alpha.f = 0.3)
pcsAndMeta$colExclude[!pcsAndMeta$excludeBasedOnPredictionCellline2] <- adjustcolor("aquamarine", alpha.f = 0.3)
plotOrderExclude <- order(pcsAndMeta$excludeBasedOnPredictionCellline2)
set.seed(42)
plotOrderExclude <- sample(nrow(pcsAndMeta), nrow(pcsAndMeta))
pcsAndMeta$cexExclude <- 0.6
pcsAndMeta$cexExclude[!pcsAndMeta$excludeBasedOnPredictionCellline2] <- 0.6

table(pcsAndMeta$colExclude)

table(pcsAndMeta$Tissue[!pcsAndMeta$excludeBasedOnPredictionCellline], pcsAndMeta$Cancer[!pcsAndMeta$excludeBasedOnPredictionCellline])
table(pcsAndMeta$Tissue[pcsAndMeta$excludeBasedOnPredictionCellline], pcsAndMeta$Cancer[pcsAndMeta$excludeBasedOnPredictionCellline])

table(pcsAndMeta$Cancer[!pcsAndMeta$excludeBasedOnPredictionCellline])
table(pcsAndMeta$Cellline[!pcsAndMeta$excludeBasedOnPredictionCellline])

#rpng(width = 800, height = 800)
#pdf(file = "test.pdf")
#plot(pcsAndMeta[plotOrderExclude,"PC_1"], pcsAndMeta[plotOrderExclude,"PC_2"], col = pcsAndMeta$colExclude[plotOrderExclude], cex = 0.3, pch = 16)
#dev.off()

#rpng(width = 1200, height = 400)

pc<-453
for(pc in c(1,3:50)){

  pcName <- paste0("PC_",pc)

  png(paste0("Recount3_QC_2ndRun/IdentifyCancerCelllines/",pcName,".png"), width = 2100, height = 700)  
  
  
  layout(matrix(1:3,nrow =1))
  plot(pcsAndMeta[plotOrderTissues,pcName], pcsAndMeta[plotOrderTissues,"PC_2"], col = pcsAndMeta$col[plotOrderTissues], cex = 0.7, pch = 16, xlab = pcName, ylab = "PC_2")
 # abline(h = threshold["PC_2"], col = "goldenrod2", lwd = 3)
  #abline(v = threshold[pcName], col = "goldenrod2", lwd = 3)
  #abline(v = rescueThreshold[pcName], col = "aquamarine", lwd = 3)
  plot(pcsAndMeta[plotOrderCellline,pcName], pcsAndMeta[plotOrderCellline,"PC_2"], col = pcsAndMeta$colCelline[plotOrderCellline], cex = 0.7, pch = 16, xlab = pcName, ylab = "PC_2")
  #abline(h = threshold["PC_2"], col = "goldenrod2", lwd = 3)
  #abline(v = threshold[pcName], col = "goldenrod2", lwd = 3)
  #abline(v = rescueThreshold[pcName], col = "aquamarine", lwd = 3)
  plot(pcsAndMeta[plotOrderExclude,pcName], pcsAndMeta[plotOrderExclude,"PC_2"], col = pcsAndMeta$colExclude[plotOrderExclude], cex = pcsAndMeta$cexExclude, pch = 16, xlab = pcName, ylab = "PC_2")
  #abline(h = threshold["PC_2"], col = "goldenrod2", lwd = 3)
  #abline(v = threshold[pcName], col = "goldenrod2", lwd = 3)
  #abline(v = rescueThreshold[pcName], col = "aquamarine", lwd = 3)
  dev.off()

}


sum(is.na(pcsAndMeta$Tissue))
sum(is.na(pcsAndMeta$Tissue2))
sum(is.na(pcsAndMeta$Tissue2))

anyTissue <- pcsAndMeta$Tissue != "" | pcsAndMeta$Tissue2 != ""

table(pcsAndMeta$Cellline[anyTissue], useNA = "a")
table(pcsAndMeta$Cellline[!anyTissue], useNA = "a")

table(pcsAndMeta$Cancer[!is.na(pcsAndMeta$Cellline)], useNA = "a")
table(pcsAndMeta$Cancer[is.na(pcsAndMeta$Cellline)], useNA = "a")


sum(pcsAndMeta$Cellline, na.rm = T)

pcsAndMetaTissueVsCellline <- pcsAndMeta[!is.na(pcsAndMeta$Cellline) & !is.na(pcsAndMeta$Cancer),]


table(pcsAndMetaTissueVsCellline$Cellline)
table(pcsAndMetaTissueVsCellline$Cancer)

pcsAndMetaTissueVsCellline$class <- NA
pcsAndMetaTissueVsCellline$class[pcsAndMetaTissueVsCellline$Cellline] <- "Cell line"
pcsAndMetaTissueVsCellline$class[!pcsAndMetaTissueVsCellline$Cellline & !pcsAndMetaTissueVsCellline$Cancer] <- "Cancer"
pcsAndMetaTissueVsCellline$class[!pcsAndMetaTissueVsCellline$Cellline & pcsAndMetaTissueVsCellline$Cancer] <- "Tissue"
pcsAndMetaTissueVsCellline$class <- as.factor(pcsAndMetaTissueVsCellline$class)

table(pcsAndMetaTissueVsCellline$class, useNA = "a")


#train <- sample(nrow(pcsAndMetaTissueVsCellline), size = round(nrow(pcsAndMetaTissueVsCellline) * 0.1))

length(unique(pcsAndMetaTissueVsCellline$study[pcsAndMetaTissueVsCellline$Cellline]))
length(unique(pcsAndMetaTissueVsCellline$study[!pcsAndMetaTissueVsCellline$Cellline]))

study <- unique(pcsAndMetaTissueVsCellline$study[pcsAndMetaTissueVsCellline$Cellline])[1]

set.seed(42)
selectedCellineTraining <- sapply(unique(pcsAndMetaTissueVsCellline$study[pcsAndMetaTissueVsCellline$Cellline]), function(study){
  thisStudyCellinesSamples <- rownames(pcsAndMetaTissueVsCellline)[pcsAndMetaTissueVsCellline$Cellline & pcsAndMetaTissueVsCellline$study == study]
  selected <- thisStudyCellinesSamples[sample(length(thisStudyCellinesSamples), min(1,length(thisStudyCellinesSamples)))]
  return(selected)
})

selectedTissueTraining <- sapply(unique(pcsAndMetaTissueVsCellline$study[!pcsAndMetaTissueVsCellline$Cellline]), function(study){
  thisStudyNonCellinesSamples <- rownames(pcsAndMetaTissueVsCellline)[!pcsAndMetaTissueVsCellline$Cellline & pcsAndMetaTissueVsCellline$study == study]
  selected <- thisStudyNonCellinesSamples[sample(length(thisStudyNonCellinesSamples), min(1,length(thisStudyNonCellinesSamples)))]
  return(selected)
})

selectedCellineTraining <- do.call("c",selectedCellineTraining)
selectedTissueTraining <- do.call("c",selectedTissueTraining)

str(selectedCellineTraining)
str(selectedTissueTraining)

train <- c(selectedCellineTraining, selectedTissueTraining)
test <- rownames(pcsAndMetaTissueVsCellline)[!rownames(pcsAndMetaTissueVsCellline) %in% train]


library(glmnet)


cfit <- cv.glmnet(x = as.matrix(pcsAndMetaTissueVsCellline[train,paste0("PC_",1:570)]), y = pcsAndMetaTissueVsCellline[train, "Cellline"], family = "binomial", type.measure = "auc", keep = TRUE, nfolds = 10, parallel = F)
rpng()
plot(cfit)
dev.off()

max(abs(coef(cfit)))

cfit


rocs <- roc.glmnet(cfit$fit.preval, newy = pcsAndMetaTissueVsCellline$Cellline)
best <- cfit$index["1se",]
rpng()
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs[[best]], lwd = 2,col = "red")
dev.off()

str(rocs)

assess.glmnet(cfit, s = "lambda.1se", newx = as.matrix(pcsAndMetaTissueVsCellline[test,paste0("PC_",1:570)]),  newy = pcsAndMetaTissueVsCellline[test, "Cellline"])
  


pcsAndMeta$excludeBasedOnPredictionCellline2 <- as.logical(predict(cfit, type = "class", s = "lambda.1se", newx = as.matrix(pcsAndMeta[,paste0("PC_",1:570)]),  newy = pcsAndMeta$Cellline)[,1])

pcsAndMeta$excludeBasedOnPredictionCellline2Score <- predict(cfit, s = "lambda.1se", newx = as.matrix(pcsAndMeta[,paste0("PC_",1:570)]),  newy = pcsAndMeta$Cellline)[,1]

rpng()
boxplot(pcsAndMeta$excludeBasedOnPredictionCellline2 ~ pcsAndMeta$excludeBasedOnPredictionCellline2Score)
dev.off()

table(pcsAndMeta$excludeBasedOnPredictionCellline2, pcsAndMeta$Cellline, useNA = "a")

table(pcsAndMeta$excludeBasedOnPredictionCellline2)




getwd()
geneInfo <- read.delim("genes_Ensembl94.txt.gz")
rownames(geneInfo) <- geneInfo$Gene.stable.ID



eigenVectors2 <- eigenVectors
rownames(eigenVectors2) <- gsub("\\..+","",rownames(eigenVectors))

eigenVectors3 <- merge(geneInfo, eigenVectors2, by = 0)
dim(eigenVectors3)

eigenVectors4 <- eigenVectors3[!eigenVectors3$Chromosome.scaffold.name %in% c("X", "Y", "MT"),]
str(eigenVectors4)

eigenVectors4 <- eigenVectors4[order(eigenVectors4$Chromosome.scaffold.name, eigenVectors4$Gene.start..bp),]
str(eigenVectors4)

str(eigenVectors4)


eigenvectorAutoCor2 <- sapply(paste0("PC_",1:570), function(x){
  y <- acf(eigenVectors4[,x],  type = "correlation", plot = F, lag.max = 1000)$acf
  return(sum(y[-1]^2))
})


which.max(eigenvectorAutoCor2)
rpng()
plot(eigenvectorAutoCor2, xlab = "Component", ylab = "sum(r^2 of first 1000 lags of auto-correlation)")
dev.off()

rpng()
plot(eigenVectors4[,"PC_143"], pch = 16, cex = 0.5, col=adjustcolor("grey", alpha.f = 0.5))
dev.off()



str(eigenValues)

v <- expPcs[,1:570] %*% diag(1/sqrt(eigenValues[1:570]))
str(v)



#Sort gene info
geneInfo2 <- geneInfo[order(geneInfo$Chromosome.scaffold.name, geneInfo$Gene.start..bp.),]
#Filter gene info
genesInOrder <- geneInfo2[geneInfo2$Gene.stable.ID %in% rownames(eigenVectors2),"Gene.stable.ID"]
str(genesInOrder)




sample <- "SRR094185"

#sampleAutoCor <- matrix(nrow = nrow(expPcs), ncol = 2, dimnames = list(rownames(expPcs), c("OriginalAutoCor", "CnvAutoCor")))

eigenVectorXeigenValues <- eigenVectors2[,1:570] %*% diag(sqrt(eigenValues[1:570]))
eigenVectorXautoCor <- eigenVectors2[,1:570] %*% diag(eigenvectorAutoCor2) 

#for(sample in rownames(expPcs)[1:10]){

cl <- makeCluster(20)

clusterExport(cl, c("eigenVectorXeigenValues", "eigenVectorXautoCor", "v", "genesInOrder"))


sampleAutoCor <- parSapply(cl, rownames(expPcs), function(sample){

  originalReconstruct <- eigenVectorXeigenValues %*% t(v[sample,1:570,drop=F])
  cnvReconstruct <- eigenVectorXautoCor %*% t(v[sample,1:570,drop=F])
  
  originalAutoCor <- sum(acf(originalReconstruct[genesInOrder,1],  type = "correlation", plot = F, lag.max = 1000)$acf[-1]^2)
  cnvAutoCor <- sum(acf(cnvReconstruct[genesInOrder,1],  type = "correlation", plot = F, lag.max = 1000)$acf[-1]^2)

  return(c("OriginalAutoCor" = originalAutoCor, "CnvAutoCor" = cnvAutoCor))
  
})


sampleAutoCor <- t(sampleAutoCor)
str(sampleAutoCor)





#save(sampleAutoCor, file = "sampleAutoCor.RData")


rpng()
plot(sampleAutoCor[,"OriginalAutoCor"], sampleAutoCor[,"CnvAutoCor"])
dev.off()

table(pcsAndMeta$Tissue2, useNA = "a")
table(pcsAndMeta$Tissue, useNA = "a")

table(pcsAndMeta[(pcsAndMeta$Tissue != "" | pcsAndMeta$Tissue2 != ""), "Cancer"],useNA = "a")

pcsAndMeta$class <- "Unkown"
pcsAndMeta$class[!is.na(pcsAndMeta$Cellline) & pcsAndMeta$Cellline] <- "Cell line"
pcsAndMeta$class[pcsAndMeta$CelllineName == "iPSC"] <- "iPSC"
pcsAndMeta$class[(pcsAndMeta$Tissue != "" | pcsAndMeta$Tissue2 != "") & !is.na(pcsAndMeta$Cancer) & pcsAndMeta$Cancer] <- "Cancer"
pcsAndMeta$class[(pcsAndMeta$Tissue != "" | pcsAndMeta$Tissue2 != "") & !is.na(pcsAndMeta$Cancer) & !pcsAndMeta$Cancer] <- "Tissue"
pcsAndMeta$class <- as.factor(pcsAndMeta$class)

table(pcsAndMeta$class, useNA = "a")

write.table(pcsAndMeta, file = "tmp.txt", sep = "\t", quote = FALSE, col.names = NA)


pcsAndMeta$OriginalAutoCor <- NA
pcsAndMeta$OriginalAutoCor <- sampleAutoCor[rownames(pcsAndMeta),"OriginalAutoCor"]

pcsAndMeta$CnvAutoCor <- NA
pcsAndMeta$CnvAutoCor <- sampleAutoCor[rownames(pcsAndMeta),"CnvAutoCor"]
str(pcsAndMeta$CnvAutoCor)

library(vioplot)
pdf("autoCorrelations.pdf", width = 14)
layout(matrix(1:2, ncol = 2))
vioplot(log(pcsAndMeta[,"OriginalAutoCor"]) ~ pcsAndMeta$class, main = "Normal expression auto correlation")
vioplot( log(pcsAndMeta[,"CnvAutoCor"]) ~ pcsAndMeta$class, main = "CNV components auto correlation")
dev.off()





#Use eigenvector2 to have gene names without version number for later sorting
test <- eigenVectors2[,1:570] %*% diag(sqrt(eigenValues[1:570])) %*% t(v[sample,1:570,drop=F])
str(test)
rpng()
plot(expScale[,sample], test[,1])
dev.off()
cor.test(expScale[,sample], test[,1])



test2 <- eigenVectors2[,1:570] %*% diag(eigenvectorAutoCor2) %*% t(v[sample,1:570,drop=F])
str(test2)
rpng()
plot(expScale[,sample], test[,1])
dev.off()
cor.test(expScale[,sample], test[,1])



rpng()
x <- acf(test[genesInOrder,1],  type = "correlation", plot = T, lag.max = 1000)
dev.off()
sum(x$acf[-1]^2)

rpng()
y <- acf(test2[genesInOrder,1],  type = "correlation", plot = T, lag.max = 1000)
dev.off()
sum(y$acf[-1]^2)


rpng()
layout(matrix(1:2,nrow = 2))
plot(test[genesInOrder,1], pch = 16, cex = 0.5, col=adjustcolor("grey", alpha.f = 0.5))
plot(test2[genesInOrder,1], pch = 16, cex = 0.5, col=adjustcolor("grey", alpha.f = 0.5))
dev.off()


rpng()

dev.off()
    