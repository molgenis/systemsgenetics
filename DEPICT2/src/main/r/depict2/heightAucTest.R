setwd("C:\\UMCG\\Genetica\\Projects\\Depict2Pgs")




library(readr)


table_tmp <- read_delim("../GeneNetwork/Data31995Genes05-12-2017/PCA_01_02_2018/PathwayMatrix/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt.gz", delim = "\t", quote = "")
hpoMatrix <- as.matrix(table_tmp[,c("HP:0000098", "HP:0002721", "HP:0002665", "HP:0004313", "HP:0004332", "HP:0002846", "HP:0002715", "HP:0100763", "HP:0001743")])
rownames(hpoMatrix) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("Height_v33/heightHpoPredictions.txt", delim = "\t", quote = "")
hpoPredictions <- as.matrix(table_tmp[,-1])
rownames(hpoPredictions) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./Height_v33/Height_genePvalues.txt", delim = "\t", quote = "")
heightGeneP <- as.matrix(table_tmp[,-1])
rownames(heightGeneP) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("./Height_v33//Height_Coregulation_Enrichment_zscoreExHla.txt", delim = "\t", quote = "")
heightCoregulation <- as.matrix(table_tmp[,-1])
rownames(heightCoregulation) <- table_tmp[,1][[1]]
rm(table_tmp)



all(rownames(hpoMatrix) ==rownames(hpoPredictions))


genes <- intersect(rownames(heightGeneP)[!is.na(heightGeneP)], rownames(heightCoregulation))
str(genes)

all(genes %in% rownames(heightCoregulation))

hpoMatrix2 <- hpoMatrix[match(genes, rownames(hpoMatrix)),]
hpoPredictions2 <- hpoPredictions[match(genes, rownames(hpoPredictions)),]
heightGeneP2 <- heightGeneP[match(genes, rownames(heightGeneP)),]
heightCoregulation2 <- heightCoregulation[match(genes, rownames(heightCoregulation)),]

all(rownames(hpoMatrix2) ==names(heightCoregulation2))
all(row.names(hpoPredictions2) ==names(heightCoregulation2))

hpoTerm = "HP:0000098"

layout(matrix(1:3, nrow = 1))
boxplot(hpoPredictions2[,hpoTerm] ~ as.factor(hpoMatrix2[,hpoTerm]), main = "GADO Tall stature")
boxplot(-log10(heightGeneP2) ~ as.factor(hpoMatrix2[,hpoTerm]), main = "-log10 Height gene p-value")
boxplot(heightCoregulation2 ~ as.factor(hpoMatrix2[,hpoTerm]), main = "Height gene co-regulation Z-score")
layout(1)

library(pROC)
gado <- roc(as.factor(hpoMatrix2[,hpoTerm]), hpoPredictions2[,hpoTerm])
geneP <- roc(as.factor(hpoMatrix2[,hpoTerm]), -log10(heightGeneP2))
geneCoReg <- roc(as.factor(hpoMatrix2[,hpoTerm]), heightCoregulation2)
#geneGadoAndCoReg <- roc(as.factor(hpoMatrix2[,hpoTerm]), (heightCoregulation2 + hpoPredictions2[,hpoTerm]))

roc.test(gado, geneP)
roc.test(geneP, geneCoReg)
roc.test(gado, geneCoReg)

plot.roc(gado, col = "springgreen2", main = "Prediction of Mendelian genes")
lines.roc(geneP, col = "goldenrod2")
lines.roc(geneCoReg, col = "dodgerblue3")
#lines.roc(geneReconstruction, col = "magenta4")
legend("bottomright", legend=c("GADO", "GWAS gene p-values", "GWAS DEPICT2 gene prioritization"), col=c("springgreen2", "goldenrod2", "dodgerblue3"), lwd=2)



wilcox.test(hpoPredictions2[,hpoTerm] ~ as.factor(hpoMatrix2[,hpoTerm]))
wilcox.test(-log10(heightGeneP2) ~ as.factor(hpoMatrix2[,hpoTerm]))
x<- wilcox.test(heightCoregulation2 ~ as.factor(hpoMatrix2[,hpoTerm]))

plot(hpoPredictions2[hpoMatrix2[,hpoTerm] == 0 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 0], xlab = "GADO score tall stature", ylab = "DEPICT2 prioritization height GWAS", main = "known tall stature genes only")
plot(hpoPredictions2[hpoMatrix2[,hpoTerm] == 1 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 1], xlab = "GADO score tall stature", ylab = "DEPICT2 prioritization height GWAS", main = "known tall stature genes only")
cor.test(hpoPredictions2[hpoMatrix2[,hpoTerm] == 1 ,hpoTerm], heightCoregulation2[hpoMatrix2[,hpoTerm] == 1])



#wilcox.test((heightCoregulation + hpoPredictions2[,"HP:0000098"]) ~ as.factor(hpoMatrix2[,"HP:0000098"]))

str(x)

str(hpoMatrix)
