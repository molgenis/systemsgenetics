setwd("C:\\UMCG\\Genetica\\Projects\\Depict2Pgs")




library(readr)

#saveRDS(hpoMatrix, "../GeneNetwork/Data31995Genes05-12-2017/PCA_01_02_2018/PathwayMatrix/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.rds")


hpoMatrix <- readRDS("../GeneNetwork/Data31995Genes05-12-2017/PCA_01_02_2018/PathwayMatrix/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.rds")
str(hpoMatrix)

hpoMatrix <- hpoMatrix[,c("HP:0011153", "HP:0000622") , drop = F]
str(hpoMatrix)

table_tmp <- read_delim("dd.txt", delim = "\t", quote = "")
hpoPredictions <- as.matrix(table_tmp[,-1])
rownames(hpoPredictions) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./educational_attainment_2018_30038396_hg19_v41/educational_attainment_genePvalues.txt", delim = "\t", quote = "")
heightGeneP <- as.matrix(table_tmp[,-1])
rownames(heightGeneP) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./educational_attainment_2018_30038396_hg19_v41/educational_attainment_Coregulation_Enrichment_zscoreExHla.txt", delim = "\t", quote = "")
heightCoregulation <- as.matrix(table_tmp[,-1])
rownames(heightCoregulation) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./educational_attainment_2018_30038396_hg19_v41/educational_attainment_Coregulation_brain_Enrichment_zscoreExHla.txt", delim = "\t", quote = "")
coregulationBrain <- as.matrix(table_tmp[,-1])
rownames(coregulationBrain) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./educational_attainment_2018_30038396_hg19_v41/SRR1237983.txt", delim = "\t", quote = "")
SRR1237983 <- as.matrix(table_tmp[,-1])
rownames(SRR1237983) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./educational_attainment_2018_30038396_hg19_v41/coexp.txt", delim = "\t", quote = "")
coexp <- as.matrix(table_tmp[,-1])
rownames(coexp) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./educational_attainment_2018_30038396_hg19_v41/coexpBrain.txt", delim = "\t", quote = "")
coexpBrain <- as.matrix(table_tmp[,-1])
rownames(coexpBrain) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./educational_attainment_2018_30038396_hg19_v41/EducationAttainmentFinemappingResults.txt", delim = "\t", quote = "")
finemap <- as.matrix(table_tmp[,-1])
rownames(finemap) <- table_tmp[,1][[1]]
rm(table_tmp)


all(rownames(hpoMatrix) ==rownames(hpoPredictions))


genes <- intersect(rownames(heightGeneP)[!is.na(heightGeneP)], rownames(heightCoregulation))
str(genes)

all(genes %in% rownames(heightCoregulation))

hpoMatrix2 <- hpoMatrix[match(genes, rownames(hpoMatrix)),,drop=F]
hpoPredictions2 <- hpoPredictions[match(genes, rownames(hpoPredictions)),]
heightGeneP2 <- heightGeneP[match(genes, rownames(heightGeneP)),]
heightCoregulation2 <- heightCoregulation[match(genes, rownames(heightCoregulation)),]
coregulationBrain2 <- coregulationBrain[match(genes, rownames(coregulationBrain)),]
SRR1237983_2 <- SRR1237983[match(genes, rownames(SRR1237983)),]
coexp2 <- coexp[match(genes, rownames(coexp)),]
coexpBrain2 <- coexpBrain[match(genes, rownames(coexpBrain)),]
finemap2 <- finemap[match(genes, rownames(finemap)),]

all(rownames(hpoMatrix2) ==names(heightCoregulation2))
all(row.names(hpoPredictions2) ==names(heightCoregulation2))
all(row.names(hpoPredictions2) ==names(coregulationBrain2))
all(row.names(hpoPredictions2) ==names(SRR1237983_2))

hpoTerm = "HP:0011153"

cat(rownames(hpoMatrix2)[hpoMatrix2[,hpoTerm] == 1])

layout(matrix(1:3, nrow = 1))
boxplot(hpoPredictions2[,hpoTerm] ~ as.factor(hpoMatrix2[,hpoTerm]), main = "GADO Tall stature")
boxplot(-log10(heightGeneP2) ~ as.factor(hpoMatrix2[,hpoTerm]), main = "-log10 Height gene p-value")
boxplot(heightCoregulation2 ~ as.factor(hpoMatrix2[,hpoTerm]), main = "Gene co-regulation Z-score")
layout(1)

library(pROC)
gado <- roc(as.factor(hpoMatrix2[,hpoTerm]), hpoPredictions2[,hpoTerm])
geneP <- roc(as.factor(hpoMatrix2[,hpoTerm]), -log10(heightGeneP2))
geneCoReg <- roc(as.factor(hpoMatrix2[,hpoTerm]), heightCoregulation2)
geneCoRegBrain <- roc(as.factor(hpoMatrix2[,hpoTerm]), coregulationBrain2)
SRR1237983_2_roc <- roc(as.factor(hpoMatrix2[,hpoTerm]), SRR1237983_2)
coexp_roc <- roc(as.factor(hpoMatrix2[,hpoTerm]), coexp2)
coexpBrain_roc <- roc(as.factor(hpoMatrix2[,hpoTerm]), coexpBrain2)
finemap2_roc <- roc(as.factor(hpoMatrix2[,hpoTerm]), finemap2)
#geneGadoAndCoReg <- roc(as.factor(hpoMatrix2[,hpoTerm]), (heightCoregulation2 + hpoPredictions2[,hpoTerm]))

roc.test(gado, geneProc.test(geneP, geneCoReg)
roc.test(geneCoRegBrain, geneCoReg)
roc.test(gado, geneCoRegBrain)
roc.test(geneCoReg, SRR1237983_2_roc)
roc.test(geneCoRegBrain, coexpBrain_roc)



png("educational_attainment_2018_30038396_hg19_v41/HP_0011153.png", width = 1000, height = 1000)
par(pty="s", bty = "n", cex = 2.2, cex.axis = 1, las = 1, cex.lab = 1.2)
plot.roc(geneP, col = "goldenrod2", main = "Prediction of 'focal motor seizures' genes\n using a 'educational attainment' GWAS", mgp=c(2.6, 0.7, 0), lwd =3)
lines.roc(gado, col = "springgreen2", lwd = 3)
lines.roc(geneCoReg, col = "dodgerblue3", lwd = 3)
lines.roc(SRR1237983_2_roc, col = "magenta4")
legend("bottomright", legend=c(paste0("DEPICT2 core-gene prioritization (AUC: ", round(geneCoReg$auc,2),")"), paste0("GWAS gene p-values (AUC: ", round(geneP$auc,2),")"), paste0("GADO (AUC: ", round(gado$auc,2),")")), col=c("dodgerblue3", "goldenrod2", "springgreen2"), lwd=3, bty="n")
dev.off()


png("educational_attainment_2018_30038396_hg19_v41/HP_0011153_npc.png", width = 1000, height = 1000)
par(pty="s", bty = "n", cex = 2.2, cex.axis = 1, las = 1, cex.lab = 1.2)
plot.roc(geneP, col = "goldenrod2", main = "Prediction of 'focal motor seizures' genes\n using a 'educational attainment' GWAS", mgp=c(2.6, 0.7, 0), lwd =3)
lines.roc(gado, col = "springgreen2", lwd = 3)
lines.roc(geneCoReg, col = "dodgerblue3", lwd = 3)
lines.roc(SRR1237983_2_roc, col = "magenta4")
legend("bottomright", legend=c(paste0("DEPICT2 core-gene prioritization (AUC: ", round(geneCoReg$auc,2),")"), paste0("GWAS gene p-values (AUC: ", round(geneP$auc,2),")"), paste0("GADO (AUC: ", round(gado$auc,2),")"),paste0("N (AUC: ", round(SRR1237983_2_roc$auc,2),")")), col=c("dodgerblue3", "goldenrod2", "springgreen2", "magnenta4"), lwd=3, bty="n")
dev.off()



par(pty="s")
plot.roc(gado, col = "springgreen2", main = "Prediction of medelian genes associated to 'Intellectual disability' using educational attainment GWAS", xlim = c(0,1))
plot.roc(geneP, col = "goldenrod2", xlim = c(1,0), ylim = c(0,1))
lines.roc(geneCoReg, col = "dodgerblue3")
#lines.roc(geneCoRegBrain, col = "magenta4")
#lines.roc(SRR1237983_2_roc, col = "black")
#lines.roc(coexp_roc, col = "black")
#lines.roc(coexpBrain_roc, col = "grey")
#lines.roc(finemap2_roc, col = "black")
legend("bottomright", legend=c("GWAS gene p-values", "GWAS coregulation"), col=c( "goldenrod2", "dodgerblue3", "magenta4", "black", "grey"), lwd=2)
  #legend("bottomright", legend=c("GADO", "GWAS gene p-values", "GWAS coregulation", "GWAS coregulation brain"), col=c("springgreen2", "goldenrod2", "dodgerblue3", "magenta4", "black", "grey"), lwd=2)

plot(coexp2, finemap2, xlab = "Gene network coexpression enrichment", ylab = "Finemappnig with gene network data")
plot(-log10(heightGeneP2), finemap2, xlab = "-log10(gene p-value)", ylab = "Finemappnig with gene network data")
plot(-log10(heightGeneP2), heightCoregulation2, xlab = "-log10(gene p-value)", ylab = "Co-regulation")
plot(-log10(heightGeneP2), coregulationBrain2, xlab = "-log10(gene p-value)", ylab = "Co-regulation brain")

cor.test(-log10(heightGeneP2), finemap2)
cor.test(-log10(heightGeneP2), coregulationBrain2)

wilcox.test(hpoPredictions2[,hpoTerm] ~ as.factor(hpoMatrix2[,hpoTerm]))
wilcox.test(-log10(heightGeneP2) ~ as.factor(hpoMatrix2[,hpoTerm]))
x<- wilcox.test(heightCoregulation2 ~ as.factor(hpoMatrix2[,hpoTerm]))

plot(hpoPredictions2[hpoMatrix2[,hpoTerm] == 0 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 0], xlab = "GADO score tall stature", ylab = "DEPICT2 prioritization height GWAS", main = "known tall stature genes only")
plot(hpoPredictions2[hpoMatrix2[,hpoTerm] == 1 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 1], xlab = "GADO score tall stature", ylab = "DEPICT2 prioritization height GWAS", main = "known tall stature genes only")
cor.test(hpoPredictions2[hpoMatrix2[,hpoTerm] == 1 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 1])


plot(hpoPredictions2[hpoMatrix2[,hpoTerm] == 1 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 1], xlab = "GADO score tall stature", ylab = "DEPICT2 prioritization height GWAS", main = "known tall stature genes only")


#wilcox.test((heightCoregulation + hpoPredictions2[,"HP:0000098"]) ~ as.factor(hpoMatrix2[,"HP:0000098"]))

str(x)

str(hpoMatrix)
