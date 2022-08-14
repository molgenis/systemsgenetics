

load(file = "combinedMeta_2022_08_14.RData", verbose = T)
load(file = "Recount3_QC_2ndRun/PCA_Patrick/pcs.RData", verbose = T)

pcsAndMeta <- merge(expPcs[,1:100], combinedMeta, by = 0, all.x = T)


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
plotOrderCellline <- order((pcsAndMeta$colCelline != defaultCol) + 1)

table(pcsAndMeta$Tissue2, pcsAndMeta$Cellline, useNA = "a")

  rpng(width = 800, height = 800)
#pdf(file = "test.pdf")
plot(pcsAndMeta[plotOrderCellline,"PC_1"], pcsAndMeta[plotOrderCellline,"PC_2"], col = pcsAndMeta$colCelline[plotOrderCellline], cex = 0.3, pch = 16)
dev.off()


table(pcsAndMeta$Cellline[!is.na(pcsAndMeta$CelllineName)&pcsAndMeta$CelllineName=="iPSC"], useNA = "always")


threshold <- c("PC_1" = 0.2, "PC_2" = -0.025, "PC_24" = -0.04, "PC_27" = 0.09, "PC_32" = -0.05, "PC_62" = 0.03)
rescueThreshold <- c("PC_10" = -0.15, "PC_18" = -0.1, "PC_20" = -0.125, "PC_21" = -0.18)

pcsAndMeta$excludeBasedOnPredictionCellline <- FALSE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_2 > threshold["PC_2"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_1 > threshold["PC_1"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_24 < threshold["PC_24"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_27 > threshold["PC_27"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_32 < threshold["PC_32"]] <- TRUE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_62 > threshold["PC_62"]] <- TRUE

pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_10 < rescueThreshold["PC_10"]] <- FALSE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_18 < rescueThreshold["PC_18"]] <- FALSE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_20 < rescueThreshold["PC_20"]] <- FALSE
pcsAndMeta$excludeBasedOnPredictionCellline[pcsAndMeta$PC_21 < rescueThreshold["PC_21"]] <- FALSE

pcsAndMeta$colExclude <- adjustcolor("goldenrod2", alpha.f = 0.3)
pcsAndMeta$colExclude[!pcsAndMeta$excludeBasedOnPredictionCellline] <- adjustcolor("aquamarine", alpha.f = 0.6)
plotOrderExclude <- order(pcsAndMeta$excludeBasedOnPredictionCellline * -1)
pcsAndMeta$cexExclude <- 0.3
pcsAndMeta$cexExclude[!pcsAndMeta$excludeBasedOnPredictionCellline] <- 0.7

table(pcsAndMeta$excludeBasedOnPredictionCellline)


table(pcsAndMeta$Cancer[!pcsAndMeta$excludeBasedOnPredictionCellline])
table(pcsAndMeta$Cellline[!pcsAndMeta$excludeBasedOnPredictionCellline])

#rpng(width = 800, height = 800)
#pdf(file = "test.pdf")
#plot(pcsAndMeta[plotOrderExclude,"PC_1"], pcsAndMeta[plotOrderExclude,"PC_2"], col = pcsAndMeta$colExclude[plotOrderExclude], cex = 0.3, pch = 16)
#dev.off()

#rpng(width = 1200, height = 400)

pc<-62
for(pc in c(1,3:100)){

  pcName <- paste0("PC_",pc)

  png(paste0("Recount3_QC_2ndRun/IdentifyCancerCelllines/",pcName,".png"), width = 2100, height = 700)  
  
  
  layout(matrix(1:3,nrow =1))
  plot(pcsAndMeta[plotOrderTissues,pcName], pcsAndMeta[plotOrderTissues,"PC_2"], col = pcsAndMeta$col[plotOrderTissues], cex = 0.7, pch = 16, xlab = pcName, ylab = "PC_2")
  abline(h = threshold["PC_2"], col = "goldenrod2", lwd = 3)
  abline(v = threshold[pcName], col = "goldenrod2", lwd = 3)
  abline(v = rescueThreshold[pcName], col = "aquamarine", lwd = 3)
  plot(pcsAndMeta[plotOrderCellline,pcName], pcsAndMeta[plotOrderCellline,"PC_2"], col = pcsAndMeta$colCelline[plotOrderCellline], cex = 0.7, pch = 16, xlab = pcName, ylab = "PC_2")
  abline(h = threshold["PC_2"], col = "goldenrod2", lwd = 3)
  abline(v = threshold[pcName], col = "goldenrod2", lwd = 3)
  abline(v = rescueThreshold[pcName], col = "aquamarine", lwd = 3)
  plot(pcsAndMeta[plotOrderExclude,pcName], pcsAndMeta[plotOrderExclude,"PC_2"], col = pcsAndMeta$colExclude[plotOrderExclude], cex = pcsAndMeta$cexExclude, pch = 16, xlab = pcName, ylab = "PC_2")
  abline(h = threshold["PC_2"], col = "goldenrod2", lwd = 3)
  abline(v = threshold[pcName], col = "goldenrod2", lwd = 3)
  abline(v = rescueThreshold[pcName], col = "aquamarine", lwd = 3)
  dev.off()

}
