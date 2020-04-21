setwd("C:\\UMCG\\Genetica\\Projects\\Depict2Pgs")




library(readr)

hpoMatrix <- readRDS("../GeneNetwork/Data31995Genes05-12-2017/PCA_01_02_2018/PathwayMatrix/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.rds")

hpoTerm = "HP:0000002"#HP:0001511
hpoMatrix <- hpoMatrix[,hpoTerm, drop = F]

table_tmp <- read_delim("hpo_selection.txt", delim = "\t", quote = "")
hpoPredictions <- as.matrix(table_tmp[,-1])
rownames(hpoPredictions) <- table_tmp[,1][[1]]
rm(table_tmsp)

table_tmp <- read_delim("./height_2018_30124842_hg19_48/height_2018_30124842_hg19_genePvalues.txt", delim = "\t", quote = "")
heightGeneP <- as.matrix(table_tmp[,-1])
rownames(heightGeneP) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./height_2018_30124842_hg19_48/height_2018_30124842_hg19_intermediates_Coregulation_Enrichment_zscoreExHla.txt", delim = "\t", quote = "")
heightCoregulation <- as.matrix(table_tmp[,-1])
rownames(heightCoregulation) <- table_tmp[,1][[1]]
rm(table_tmp)


pli <- read_delim("./gnomad.v2.1.1.lof_metrics.by_gene.txt.gz", delim = "\t", quote = "")

ensgGeneTranscript <- read_delim("./ensgV83_GeneTranscript.txt", delim = "\t", quote = "")

all(pli$transcript %in% ensgGeneTranscript$`Ensembl Transcript ID`)
pli$geneEnsg <- ensgGeneTranscript$`Ensembl Gene ID`[match(pli$transcript, ensgGeneTranscript$`Ensembl Transcript ID`)]

all(rownames(hpoMatrix) ==rownames(hpoPredictions))


genes <- intersect(rownames(heightGeneP)[!is.na(heightGeneP)], rownames(heightCoregulation))
str(genes)

all(genes %in% rownames(heightCoregulation))

hpoMatrix2 <- hpoMatrix[match(genes, rownames(hpoMatrix)),,drop = F]
hpoPredictions2 <- hpoPredictions[match(genes, rownames(hpoPredictions)),]
heightGeneP2 <- heightGeneP[match(genes, rownames(heightGeneP)),]
heightCoregulation2 <- heightCoregulation[match(genes, rownames(heightCoregulation)),]
#heightLude2 <- heightLude[match(genes, rownames(heightLude)),1]

all(rownames(hpoMatrix2) ==names(heightCoregulation2))
all(row.names(hpoPredictions2) ==names(heightCoregulation2))
#all(row.names(hpoPredictions2) ==names(heightLude2))



layout(matrix(1:3, nrow = 1))
boxplot(hpoPredictions2[,hpoTerm] ~ as.factor(hpoMatrix2[,hpoTerm]), main = "GADO Tall stature")
boxplot(-log10(heightGeneP2) ~ as.factor(hpoMatrix2[,hpoTerm]), main = "-log10 Height gene p-value")
boxplot(heightCoregulation2 ~ as.factor(hpoMatrix2[,hpoTerm]), main = "Height gene co-regulation Z-score")
layout(1)



library(pROC)
gado <- roc(as.factor(hpoMatrix2[,hpoTerm]), hpoPredictions2[,hpoTerm], ci = T)
geneP <- roc(as.factor(hpoMatrix2[,hpoTerm]), -log10(heightGeneP2))
(geneCoReg <- roc(as.factor(hpoMatrix2[,hpoTerm]), heightCoregulation2))
#LudeRoc <- roc(as.factor(hpoMatrix2[,hpoTerm]), heightLude2)
#geneGadoAndCoReg <- roc(as.factor(hpoMatrix2[,hpoTerm]), (heightCoregulation2 + hpoPredictions2[,hpoTerm]))

roc.test(gado, geneP)
roc.test(geneP, geneCoReg)
roc.test(gado, geneCoReg)
  
  png("height_2018_30124842_hg19_48/HP_0001519.png", width = 1000, height = 1000)
  par(pty="s", bty = "n", cex = 2.2, cex.axis = 1, las = 1, cex.lab = 1.2)
  plot.roc(geneP, col = "goldenrod2", main = "Prediction of 'Disproportionate tall stature' genes\n using a 'height' GWAS", mgp=c(2.6, 0.7, 0), lwd = 3)
  lines.roc(gado, col = "springgreen2", lwd = 3)
  lines.roc(geneCoReg, col = "dodgerblue3", lwd = 3)
  legend("bottomright", legend=c(paste0("DEPICT2 core-gene prioritization (AUC: ", round(geneCoReg$auc,2),")"), paste0("GWAS gene p-values (AUC: ", round(geneP$auc,2),")"), paste0("GADO (AUC: ", round(gado$auc,2),")")), col=c("dodgerblue3", "goldenrod2", "springgreen2"), lwd=3, bty="n")
  dev.off()
  
  
  
  pdf("height_2018_30124842_hg19_48/HP_0001519.pdf", width = 10, height = 10)
  par(pty="s", bty = "n", cex = 2.2, cex.axis = 1, las = 1, cex.lab = 1.2, xpd = NA)
  plot.roc(geneP, col = "goldenrod2", main = "Prediction of 'Disproportionate tall stature' genes\n using a 'height' GWAS", mgp=c(2.6, 0.7, 0), lwd = 3)
  lines.roc(gado, col = "springgreen2", lwd = 3)
  lines.roc(geneCoReg, col = "dodgerblue3", lwd = 3)
  legend("bottomright", legend=c(paste0("DEPICT2 core-gene prioritization (AUC: ", round(geneCoReg$auc,2),")"), paste0("GWAS gene p-values (AUC: ", round(geneP$auc,2),")"), paste0("GADO (AUC: ", round(gado$auc,2),")")), col=c("dodgerblue3", "goldenrod2", "springgreen2"), lwd=3, bty="n")
  dev.off()
  
  
  
  
  png("height_2018_30124842_hg19_48/HP_0001519_A.png", width = 1000, height = 1000)
  par(pty="s", bty = "n", cex = 2.2, cex.axis = 1, las = 1, cex.lab = 1.2)
  plot.roc(geneP, col = "goldenrod2", main = "Prediction of 'Disproportionate tall stature' genes\n using a 'height' GWAS", mgp=c(2.6, 0.7, 0), lwd = 3)
  legend("bottomright", legend=c(paste0("GWAS gene p-values (AUC: ", round(geneP$auc,2),")")), col=c("goldenrod2", "dodgerblue3", "springgreen2"), lwd=3, bty="n")
  dev.off()
  
  png("height_2018_30124842_hg19_48/HP_0001519_B.png", width = 1000, height = 1000)
  par(pty="s", bty = "n", cex = 2.2, cex.axis = 1, las = 1, cex.lab = 1.2)
  plot.roc(geneP, col = "goldenrod2", main = "Prediction of 'Disproportionate tall stature' genes\n using a 'height' GWAS", mgp=c(2.6, 0.7, 0), lwd = 3)
  lines.roc(geneCoReg, col = "dodgerblue3", lwd = 3)
  legend("bottomright", legend=c(paste0("DEPICT2 core-gene prioritization (AUC: ", round(geneCoReg$auc,2),")"), paste0("GWAS gene p-values (AUC: ", round(geneP$auc,2),")")), col=c("dodgerblue3", "goldenrod2", "springgreen2"), lwd=3, bty="n")
  dev.off()
  



wilcox.test(hpoPredictions2 ~ as.factor(hpoMatrix2[,hpoTerm]))
wilcox.test(-log10(heightGeneP2) ~ as.factor(hpoMatrix2[,hpoTerm]))
x<- wilcox.test(heightCoregulation2 ~ as.factor(hpoMatrix2[,hpoTerm]))

plot(hpoPredictions2[hpoMatrix2[,hpoTerm] == 0 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 0], xlab = "GADO score tall stature", ylab = "DEPICT2 prioritization height GWAS", main = "known tall stature genes only")
plot(hpoPredictions2[hpoMatrix2[,hpoTerm] == 1 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 1], xlab = "GADO score tall stature", ylab = "DEPICT2 prioritization height GWAS", main = "known tall stature genes only")
cor.test(hpoPredictions2[hpoMatrix2[,hpoTerm] == 1 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 1])



#wilcox.test((heightCoregulation + hpoPredictions2[,"HP:0000098"]) ~ as.factor(hpoMatrix2[,"HP:0000098"]))

str(x)

str(hpoMatrix)

heightCoregulation2[hpoMatrix2[,hpoTerm]==1]
heightCoregulation2[hpoMatrix2[,hpoTerm]==0]

range(heightCoregulation2)

#hpoMatrix2["ENSG00000185049",hpoTerm] <- 0
#hpoMatrix2["ENSG00000168924",hpoTerm] <- 0
#hpoMatrix2["ENSG00000109685",hpoTerm] <- 0


zcoreBreaks <- seq(-7,7,0.2)

diseaseHist <- hist(heightCoregulation2[hpoMatrix2[,hpoTerm]==1], plot = F, breaks = zcoreBreaks)
otherHist <- hist(heightCoregulation2[hpoMatrix2[,hpoTerm]==0], plot = F, breaks = zcoreBreaks)

par(mar = c(10,5,1,1), las = 1)
plot.new()
plot.window(ylim = range(diseaseHist$density, otherHist$density), xlim = range(zcoreBreaks))
axis(side =2)
axis(side =1)
title(xlab = "Core-gene score for height GWAS")
title(ylab = "Density", mgp = c(3.5,1,0))
plot.window(ylim = range(diseaseHist$density, otherHist$density), xlim = c(0,length(otherHist$density)))
barplot((diseaseHist$density), space = 0, horiz = F, col = adjustcolor("darkorange", alpha.f = 0.5), border = NA, add = T, axes = F)
barplot((otherHist$density), space = 0, horiz = F, col = adjustcolor("dodgerblue2", alpha.f = 0.5), border = NA, add = T, axes = F)

tail(sort(heightCoregulation2), n = 10)





genesWithPli <- intersect(names(heightCoregulation2), pli$geneEnsg)

pli2 <- pli[match(genesWithPli, pli$geneEnsg),]
heightCoregulation2WithPli <- heightCoregulation2[genesWithPli]

plot(heightCoregulation2WithPli, pli2$mis_z)
cor.test(heightCoregulation2WithPli, pli2$mis_z)



boxplot(heightCoregulation2WithPli[pli2$mis_z <= 0.05], heightCoregulation2WithPli[pli2$mis_z >= 0.95], names = c("pLI <= 0.05", "pLI => 0.95"), ylab = "Height GWAs core-gene predicion")

range(heightCoregulation2WithPli)

x <- cut(heightCoregulation2WithPli, seq(-7,7,2))

boxplot(pli2$pLI ~ x, ylab = "pLI", xlab = "Height GWAs core-gene predicion bins")

layout(matrix(1:3,nrow= 1))

heightCoregulation2["ENSG00000163501"]

boxplot(pli2$syn_z ~ x, ylab = "Synonymous Z", xlab = "Height GWAs core-gene predicion bins")

boxplot(pli2$mis_z ~ x, ylab = "Missense Z", xlab = "Height GWAs core-gene predicion bins")

boxplot(pli2$lof_z ~ x, ylab = "Loss of funtion Z", xlab = "Height GWAs core-gene predicion bins")





heightCoregulation2Groups <- cut(heightCoregulation2, seq(min(heightCoregulation2),max(heightCoregulation2),length.out = 11))
heightGeneP2Groups <- cut(-log10(heightGeneP2), seq(0,12,length.out = 11))
  

heightCoregulation2Groups <- cut(heightCoregulation2, quantile(heightCoregulation2, seq(0,1,0.025)))
heightGeneP2Groups <- cut(-log10(heightGeneP2), quantile(-log10(heightGeneP2), seq(0,1,0.025))[1:40])


coreGeneEnrichment <- aggregate(hpoMatrix2[,hpoTerm], by = list(heightCoregulation2Groups), FUN = function(x){(sum(x) * 100)/length(x)})
gwasEnrichment <- aggregate(hpoMatrix2[,hpoTerm], by = list(heightGeneP2Groups), FUN = function(x){(sum(x) * 100)/length(x)})

aggregate(hpoMatrix2[,hpoTerm], by = list(heightCoregulation2Groups), FUN = length)
aggregate(hpoMatrix2[,hpoTerm], by = list(heightGeneP2Groups), FUN = length)

layout(matrix(1:2,nrow= 1))
barplot(gwasEnrichment$x, names.arg = gwasEnrichment$Group.1, ylab = "Percentage of genes assocaited to Abnormality of body height", xlab = "-log10(GWAS gene p-value", ylim = c(0,15))
barplot(coreGeneEnrichment$x, names.arg = coreGeneEnrichment$Group.1, ylab = "Percentage of genes assocaited to Abnormality of body height", xlab = "Core gene prioritzation", ylim = c(0,15))



