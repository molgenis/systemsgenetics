

load(file = "combinedMeta_2022_08_14.RData", verbose = T)
load(file = "Recount3_QC_2ndRun/PCA_Patrick/pcs.RData", verbose = T)

pcsAndMeta <- merge(expPcs[,1:570], combinedMeta, by = 0, all.x = T)
save(pcsAndMeta, file = "DataForLasso.RData")

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

pc<-1
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


train <- sample(nrow(pcsAndMetaTissueVsCellline), size = round(nrow(pcsAndMetaTissueVsCellline) * 0.8))



library(glmnet)
#cvfit <- cv.glmnet(x = as.matrix(pcsAndMetaTissueVsCellline[,paste0("PC_",1:570)]), y = pcsAndMetaTissueVsCellline$class, family = "multinomial")
summary(cvfit)

rpng()
plot(cvfit, xvar = "lambda", label = TRUE, type.coef = "coef")
dev.off()

rpng()
plot(cvfit)
dev.off()
str(cvfit)


cnf <- confusion.glmnet(cvfit,newx = as.matrix(pcsAndMetaTissueVsCellline[,paste0("PC_",1:10)]), newy = pcsAndMetaTissueVsCellline$class)
cnf








require(doMC)
registerDoMC(cores = 20)

cfit <- cv.glmnet(x = as.matrix(pcsAndMetaTissueVsCellline[train,paste0("PC_",1:570)]), y = pcsAndMetaTissueVsCellline$Cellline[train], family = "binomial", type.measure = "auc", keep = TRUE, nfolds = 10, parallel = F)
rpng()
plot(cfit)
dev.off()

coef(cfit, s = exp(-3))
sum(!coef(cfit, s = exp(-3))==0)
log(0.0005)
rocs <- roc.glmnet(cfit$fit.preval, newy = pcsAndMetaTissueVsCellline$Cellline)
log(0.01)
cfit$cvm

best <- cvfit$index["min",]
rpng()
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs[[best]], lwd = 2,col = "red")
dev.off()

str(rocs)

log(0.1)
exp(-6)
assess.glmnet(cfit, s = exp(-3), newx = as.matrix(pcsAndMetaTissueVsCellline[-train,paste0("PC_",1:570)]),  newy = pcsAndMetaTissueVsCellline$Cellline[-train])
  


pcsAndMeta$excludeBasedOnPredictionCellline2 <- as.logical(predict(cfit, type = "class", s = "lambda.1se", newx = as.matrix(pcsAndMeta[,paste0("PC_",1:570)]),  newy = pcsAndMeta$Cellline)[,1])


table(pcsAndMeta$excludeBasedOnPredictionCellline2, pcsAndMeta$Cellline, useNA = "a")
