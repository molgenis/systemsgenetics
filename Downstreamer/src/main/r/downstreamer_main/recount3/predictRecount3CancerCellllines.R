library(parallel)

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

load(file = "Metadata/combinedMeta_2022_09_15.RData", verbose = T)
load(file = "Recount3_QC_2ndRun/PCA_Patrick/pcs.RData", verbose = T)
load(file = "Recount3_QC_2ndRun/PCA_Patrick/eigen.RData", verbose = T)

#save.image("predictionSession.RData")
load("predictionSession.RData")

explainedVariance <- eigenValues * 100 / sum(eigenValues)
(compsToUse <- which(cumsum(explainedVariance)>= 80)[1])


sum(rownames(expPcs[,1:compsToUse]) %in% rownames(combinedMeta))
sum(!rownames(expPcs[,1:compsToUse]) %in% rownames(combinedMeta))


pcsAndMeta <- merge(expPcs[,1:compsToUse], combinedMeta, by = 0)
rownames(pcsAndMeta) <- pcsAndMeta$Row.names
colnames(pcsAndMeta)

#save(pcsAndMeta, compsToUse, file = "DataForPredictions.RData")
load(file = "DataForPredictions.RData")

colnames(pcsAndMeta)
pcsAndMeta <- pcsAndMeta[!pcsAndMeta$exclude,]
dim(pcsAndMeta)

table(pcsAndMeta$excludeBasedOnPredictionCancer, pcsAndMeta$excludeBasedOnPredictionCellline2)

#write.table(merge(combinedMeta, expPcs[,1:100], by = 0, all.y = T), file = "Metadata/pcsAndAnnotations.txt", sep = "\t", quote = FALSE, col.names = NA)

dim(pcsAndMeta)


defaultCol <- adjustcolor("grey", alpha.f = 0.6)
tissueCol <- read.delim("Recount3_QC_2ndRun/SRA_Studies_Annotations_Patrick/Annotations_color2.txt")

#write.table(tissueCol, file = "Recount3_QC_2ndRun/SRA_Studies_Annotations_Patrick/Annotations_color2.txt", quote = F, row.names = F, sep = "\t")

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

pc<453
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







###################### Celline predictions


pcsAndMetaTissueVsCellline <- pcsAndMeta[(!is.na(pcsAndMeta$Cellline) & pcsAndMeta$Cellline) | ((pcsAndMeta$Tissue != "" | pcsAndMeta$Tissue2 != "") & !is.na(pcsAndMeta$Cancer) & !pcsAndMeta$Cancer ),]



table(pcsAndMetaTissueVsCellline$Cellline, pcsAndMetaTissueVsCellline$Cancer,useNA = "a")



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


cfit <- cv.glmnet(x = as.matrix(pcsAndMetaTissueVsCellline[train,c(paste0("PC_",1:compsToUse))]), y = pcsAndMetaTissueVsCellline[train, "Cellline"], family = "binomial", type.measure = "auc", keep = TRUE, nfolds = 10, parallel = F)
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

assess.glmnet(cfit, s = "lambda.1se", newx = as.matrix(pcsAndMetaTissueVsCellline[test,c(paste0("PC_",1:compsToUse))]),  newy = pcsAndMetaTissueVsCellline[test, "Cellline"])
  


pcsAndMeta$excludeBasedOnPredictionCellline2 <- as.logical(predict(cfit, type = "class", s = "lambda.1se", newx = as.matrix(pcsAndMeta[,c(paste0("PC_",1:compsToUse))]),  newy = pcsAndMeta$Cellline)[,1])

pcsAndMeta$excludeBasedOnPredictionCellline2Score <- predict(cfit, s = "lambda.1se", type = "response", newx = as.matrix(pcsAndMeta[,c(paste0("PC_",1:compsToUse))]),  newy = pcsAndMeta$Cellline)[,1]
pcsAndMeta$excludeBasedOnPredictionCellline2Scoreb <- predict(cfit, s = "lambda.1se", newx = as.matrix(pcsAndMeta[,c(paste0("PC_",1:compsToUse))]),  newy = pcsAndMeta$Cellline)[,1]









############Cancer prediction


library(pROC)


pcsAndMetaTissueVsCancer <- pcsAndMeta[(pcsAndMeta$Tissue != "" | pcsAndMeta$Tissue2 != "") & !is.na(pcsAndMeta$Cancer),]

roc(pcsAndMetaTissueVsCancer$Cancer, pcsAndMetaTissueVsCancer$OriginalAutoCor )
roc(pcsAndMetaTissueVsCancer$Cancer, pcsAndMetaTissueVsCancer$CnvAutoCor )
table(pcsAndMetaTissueVsCancer$Cancer, useNA = "a")

pcsAndMetaTissueVsCancer$Cancer2 <- as.factor(pcsAndMetaTissueVsCancer$Cancer)
table(pcsAndMetaTissueVsCancer$Cancer2, useNA = "a")

cancerModel <- glm(Cancer2 ~ CnvAutoCor, family=binomial(link='logit'), data = pcsAndMetaTissueVsCancer)



pcsAndMeta$excludeBasedOnPredictionCancer <- predict(cancerModel, newdata = data.frame(pcsAndMeta[,"CnvAutoCor",drop=F]), type = "response") > 0.5
pcsAndMeta$excludeBasedOnPredictionCancerScore <- predict(cancerModel, newdata = data.frame(pcsAndMeta[,"CnvAutoCor",drop=F]), type = "response")









##############Plotting
pcsAndMeta$colCelline <- defaultCol
pcsAndMeta$colCelline[!is.na(pcsAndMeta[,"Cellline"]) & pcsAndMeta[,"Cellline"]] <-  adjustcolor("magenta", alpha.f = 0.6)
pcsAndMeta$colCelline[!is.na(pcsAndMeta[,"Cellline"]) & !pcsAndMeta[,"Cellline"]] <-  adjustcolor("royalblue1", alpha.f = 0.6)
pcsAndMeta$colCelline[!is.na(pcsAndMeta[,"Cancer"]) & pcsAndMeta[,"Cancer"]] <-  adjustcolor("forestgreen", alpha.f = 0.6)



defaultCol <- adjustcolor("grey", alpha.f = 0.3)
pcsAndMeta$classCol <- defaultCol
pcsAndMeta$classCol[!is.na(pcsAndMeta$Cellline) & pcsAndMeta$Cellline] <- adjustcolor("magenta", alpha.f = 0.6)
pcsAndMeta$classCol[pcsAndMeta$CelllineName == "iPSC"] <- adjustcolor("grey", alpha.f = 0.6)
pcsAndMeta$classCol[(pcsAndMeta$Tissue != "" | pcsAndMeta$Tissue2 != "") & !is.na(pcsAndMeta$Cancer) & pcsAndMeta$Cancer] <- adjustcolor("forestgreen", alpha.f = 0.6)
pcsAndMeta$classCol[(pcsAndMeta$Tissue != "" | pcsAndMeta$Tissue2 != "") & !is.na(pcsAndMeta$Cancer) & !pcsAndMeta$Cancer] <- adjustcolor("royalblue1", alpha.f = 0.6)


pcsAndMeta$classCex <- 0.5
pcsAndMeta$classCex[!is.na(pcsAndMeta$Cellline) & pcsAndMeta$Cellline] <- 0.6
pcsAndMeta$classCex[pcsAndMeta$CelllineName == "iPSC"] <- 0.8
pcsAndMeta$classCex[(pcsAndMeta$Tissue != "" | pcsAndMeta$Tissue2 != "") & !is.na(pcsAndMeta$Cancer) & pcsAndMeta$Cancer] <- 0.6
pcsAndMeta$classCex[(pcsAndMeta$Tissue != "" | pcsAndMeta$Tissue2 != "") & !is.na(pcsAndMeta$Cancer) & !pcsAndMeta$Cancer] <- 0.6

set.seed(42)
plotOrderClass <- sample(nrow(pcsAndMeta), nrow(pcsAndMeta))
plotOrderClass[pcsAndMeta$classCol == defaultCol] <- 1

rpng(width = 800, height = 800)
 plot(pcsAndMeta[plotOrderClass,"excludeBasedOnPredictionCellline2Score"], pcsAndMeta[plotOrderClass,"excludeBasedOnPredictionCancerScore"], col = pcsAndMeta$classCol[plotOrderClass], cex = pcsAndMeta$classCex[plotOrderClass], pch = 16, 
     xlab = "Cell line logistic regression score using components", 
     ylab= "Cancer logistic regression scores using auto correlation")
abline(v = 0.5, lwd = 2, col = adjustcolor("magenta", alpha.f = 0.4))
abline(h = 0.5, lwd = 2, col = adjustcolor("forestgreen", alpha.f = 0.4))
dev.off()

rpng(width = 800, height = 800)
pdf("cancerCellLine.pdf", useDingbats = F)
par(xpd=NA)
plot(pcsAndMeta[plotOrderClass,"excludeBasedOnPredictionCellline2Scoreb"], log(pcsAndMeta[plotOrderClass,"CnvAutoCor"]), col = pcsAndMeta$classCol[plotOrderClass], cex = pcsAndMeta$classCex[plotOrderClass], pch = 16, 
     xlab = "Sampe cell line prediction", 
     ylab= "Sample copy number variant load",
     frame = F)
abline(v = min(pcsAndMeta[pcsAndMeta$excludeBasedOnPredictionCellline2,"excludeBasedOnPredictionCellline2Scoreb"]), lwd = 2, col = adjustcolor("magenta", alpha.f = 0.4))
abline(h = log(min(pcsAndMeta[pcsAndMeta$excludeBasedOnPredictionCancer,"CnvAutoCor"]))
       , lwd = 2, col = adjustcolor("forestgreen", alpha.f = 0.4))
dev.off()











rpng(width = 800, height = 800)
plot(pcsAndMeta[plotOrderClass,"excludeBasedOnPredictionCellline2Scoreb"], log(pcsAndMeta[plotOrderClass,"CnvAutoCor"]), col = pcsAndMeta$classCol[plotOrderClass], cex = pcsAndMeta$classCex[plotOrderClass], pch = 16, 
     xlab = "Sampe cell line prediction", 
     ylab= "Sample copy number variant load",
     frame = F)
abline(v = min(pcsAndMeta[pcsAndMeta$excludeBasedOnPredictionCellline2,"excludeBasedOnPredictionCellline2Scoreb"]), lwd = 2, col = adjustcolor("magenta", alpha.f = 0.4))
abline(h = log(min(pcsAndMeta[pcsAndMeta$excludeBasedOnPredictionCancer,"CnvAutoCor"]))
       , lwd = 2, col = adjustcolor("forestgreen", alpha.f = 0.4))
dev.off()













rpng(width = 800, height = 800)
plot(pcsAndMeta[plotOrderTissues,"excludeBasedOnPredictionCellline2Score"], pcsAndMeta[plotOrderTissues,"excludeBasedOnPredictionCancerScore"], col = pcsAndMeta$col[plotOrderTissues], cex =0.6  , pch = 16, 
     xlab = "Cell line logistic regression score using components", 
     ylab= "Cancer logistic regression scores using auto correlation")
abline(v = 0.5, lwd = 2, col = adjustcolor("magenta", alpha.f = 0.4))
abline(h = 0.5, lwd = 2, col = adjustcolor("forestgreen", alpha.f = 0.4))
dev.off()



table(pcsAndMeta$Tissue[pcsAndMeta$excludeBasedOnPredictionCellline2 & !(!is.na(pcsAndMeta$Cancer) & pcsAndMeta$Cancer)])
sort(table(pcsAndMeta$study[pcsAndMeta$excludeBasedOnPredictionCellline2 & !(!is.na(pcsAndMeta$Cancer) & pcsAndMeta$Cancer)]))

table(pcsAndMeta$excludeBasedOnPredictionCellline2, pcsAndMeta$excludeBasedOnPredictionCancer)


table(pcsAndMeta[pcsAndMeta[,"excludeBasedOnPredictionCellline2Score"] > 1 & pcsAndMeta[,"study"] == "GTEx","gtex.smtsd"])


pcsAndMeta$selectedSamples <- !pcsAndMeta$excludeBasedOnPredictionCellline2 & !pcsAndMeta$excludeBasedOnPredictionCancer & !(!is.na(pcsAndMeta$Cancer) & pcsAndMeta$Cancer) & !(!is.na(pcsAndMeta$Cellline) & pcsAndMeta$Cellline)
table(pcsAndMeta$selectedSamples)
table(pcsAndMeta$selectedSamples, useNA = "a")



table(pcsAndMeta$excludeBasedOnPredictionCellline2, pcsAndMeta$Cellline, useNA = "a")

table(pcsAndMeta$excludeBasedOnPredictionCellline2)
table(pcsAndMeta$excludeBasedOnPredictionCancer)




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


eigenvectorAutoCor2 <- sapply(paste0("PC_",1:compsToUse), function(x){
  y <- acf(eigenVectors4[,x],  type = "correlation", plot = F, lag.max = 1000)$acf
  return(sum(y[-1]^2))
})


which.max(eigenvectorAutoCor2)
rpng()
plot(eigenvectorAutoCor2, xlab = "Component", ylab = "sum(r^2 of first 1000 lags of auto-correlation)")
dev.off()

rpng(width = 800, height = 400)
plot(eigenVectors4[,"PC_143"], pch = 16, cex = 0.5, col=adjustcolor("grey", alpha.f = 0.5), ylab = "Eigen vector score comp 143")
dev.off()



str(eigenValues)

v <- expPcs[,1:compsToUse] %*% diag(1/sqrt(eigenValues[1:compsToUse]))
str(v)



#Sort gene info
geneInfo2 <- geneInfo[order(geneInfo$Chromosome.scaffold.name, geneInfo$Gene.start..bp.),]
#Filter gene info
genesInOrder <- geneInfo2[geneInfo2$Gene.stable.ID %in% rownames(eigenVectors2),"Gene.stable.ID"]
str(genesInOrder)




sample <- "SRR5553774"

#sampleAutoCor <- matrix(nrow = nrow(expPcs), ncol = 2, dimnames = list(rownames(expPcs), c("OriginalAutoCor", "CnvAutoCor")))

#Use eigenvector2 to have gene names without version number for later sorting
eigenVectorXeigenValues <- eigenVectors2[,1:compsToUse] %*% diag(sqrt(eigenValues[1:compsToUse]))
eigenVectorXautoCor <- eigenVectors2[,1:compsToUse] %*% diag(eigenvectorAutoCor2) 

#for(sample in rownames(expPcs)[1:10]){

cl <- makeCluster(20)

clusterExport(cl, c("eigenVectorXeigenValues", "eigenVectorXautoCor", "v", "genesInOrder", "compsToUse"))


sampleAutoCor <- parSapply(cl, rownames(expPcs), function(sample){

  originalReconstruct <- eigenVectorXeigenValues %*% t(v[sample,1:compsToUse,drop=F])
  cnvReconstruct <- eigenVectorXautoCor %*% t(v[sample,1:compsToUse,drop=F])
  
  originalAutoCor <- sum(acf(originalReconstruct[genesInOrder,1],  type = "correlation", plot = F, lag.max = 1000)$acf[-1]^2)
  cnvAutoCor <- sum(acf(cnvReconstruct[genesInOrder,1],  type = "correlation", plot = F, lag.max = 1000)$acf[-1]^2)

  return(c("OriginalAutoCor" = originalAutoCor, "CnvAutoCor" = cnvAutoCor))
  
})

rpng()
par(mar=c(3,5,0.1,0.1))
layout(matrix(1:2, nrow = 2))
plot(originalReconstruct, pch = 16, cex = 0.2, col=adjustcolor("grey", alpha.f = 0.5))
plot(cnvReconstruct, pch = 16, cex = 0.2, col=adjustcolor("grey", alpha.f = 0.5))
dev.off()
sampleAutoCor[sample,]

sampleAutoCor <- t(sampleAutoCor)
str(sampleAutoCor)

tail(sort(sampleAutoCor[sampleAutoCor[,2] >10 & sampleAutoCor[,1] <2,2]))

 

#save(sampleAutoCor, file = "sampleAutoCor.RData")
load(file = "sampleAutoCor.RData", verbose = T)


rpng()
plot(sampleAutoCor[,"OriginalAutoCor"], sampleAutoCor[,"CnvAutoCor"], pch = 16, cex = 0.4, col=adjustcolor("grey", alpha.f = 0.5))
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

#write.table(pcsAndMeta, file = "tmp.txt", sep = "\t", quote = FALSE, col.names = NA)


pcsAndMeta$OriginalAutoCor <- NA
pcsAndMeta$OriginalAutoCor <- sampleAutoCor[rownames(pcsAndMeta),"OriginalAutoCor"]

pcsAndMeta$CnvAutoCor <- NA
pcsAndMeta$CnvAutoCor <- sampleAutoCor[rownames(pcsAndMeta),"CnvAutoCor"]
range(pcsAndMeta$CnvAutoCor)

library(vioplot)
pdf("autoCorrelations.pdf", width = 14)
rpng()
layout(matrix(1:2, ncol = 2))
vioplot(log(pcsAndMeta[,"OriginalAutoCor"]) ~ pcsAndMeta$class, main = "Normal expression auto correlation")
vioplot( log(pcsAndMeta[,"CnvAutoCor"]) ~ pcsAndMeta$class, main = "CNV components auto correlation")
dev.off()



