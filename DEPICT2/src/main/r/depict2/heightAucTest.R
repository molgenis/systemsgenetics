setwd("C:\\UMCG\\Genetica\\Projects\\Depict2Pgs")




library(readr)

hpoMatrix <- readRDS("../GeneNetwork/Data31995Genes05-12-2017/PCA_01_02_2018/PathwayMatrix/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.rds")

hpoTerm = "HP:0001519"
hpoMatrix <- hpoMatrix[,hpoTerm, drop = F]


table_tmp <- read_delim("hpo_selection.txt", delim = "\t", quote = "")
hpoPredictions <- as.matrix(table_tmp[,-1])
rownames(hpoPredictions) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./height_2018_30124842_hg19_48/height_2018_30124842_hg19_genePvalues.txt", delim = "\t", quote = "")
heightGeneP <- as.matrix(table_tmp[,-1])
rownames(heightGeneP) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./height_2018_30124842_hg19_48/height_2018_30124842_hg19_intermediates_Coregulation_Enrichment_zscoreExHla.txt", delim = "\t", quote = "")
heightCoregulation <- as.matrix(table_tmp[,-1])
rownames(heightCoregulation) <- table_tmp[,1][[1]]
rm(table_tmp)




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
geneCoReg <- roc(as.factor(hpoMatrix2[,hpoTerm]), heightCoregulation2)
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
  



wilcox.test(hpoPredictions2[,hpoTerm] ~ as.factor(hpoMatrix2[,hpoTerm]))
wilcox.test(-log10(heightGeneP2) ~ as.factor(hpoMatrix2[,hpoTerm]))
x<- wilcox.test(heightCoregulation2 ~ as.factor(hpoMatrix2[,hpoTerm]))

plot(hpoPredictions2[hpoMatrix2[,hpoTerm] == 0 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 0], xlab = "GADO score tall stature", ylab = "DEPICT2 prioritization height GWAS", main = "known tall stature genes only")
plot(hpoPredictions2[hpoMatrix2[,hpoTerm] == 1 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 1], xlab = "GADO score tall stature", ylab = "DEPICT2 prioritization height GWAS", main = "known tall stature genes only")
cor.test(hpoPredictions2[hpoMatrix2[,hpoTerm] == 1 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 1])



#wilcox.test((heightCoregulation + hpoPredictions2[,"HP:0000098"]) ~ as.factor(hpoMatrix2[,"HP:0000098"]))

str(x)

str(hpoMatrix)
