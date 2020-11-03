x <- read.table("~/Desktop/depict2/output/maf_filtered/v51/height_2018_30124842_hg19_genePvalues.txt", header=T, row.names=1)
y <- read.table("~/Desktop/depict2/output/height_paper/height_2018_30124842_hg19_genePvalues.txt", header=T, row.names=1)
x[,1] <- -log10(x[,1])
y[,1] <- -log10(y[,1])

coreg.x <- read.table("~/Desktop/depict2/output/specific_coreg_extractions/height_FNDC3B_YAP1_skewness.txt", row.names=1, header=T)
coreg.y <- read.table("~/Desktop/depict2/output/specific_coreg_extractions/height_FNDC3B_YAP1_protein_coding.txt", row.names=1, header=T)

ol <- intersect()

par(mfrow=c(2,2))
m <- lm(coreg.x[rownames(x),"ENSG00000075420"] ~ x[,1])
plot(x[,1], coreg.x[rownames(x),"ENSG00000075420"],
     main=paste0("FNDC3B OLS model P: ", format(summary(m)$coefficients[2,4], scientific=T, digits=2)),
     xlab="-log10(gene P)",
     ylab="Coregulation skewness")
abline(m, col="blue")

m <- lm(coreg.y[rownames(y),"ENSG00000075420"] ~ y[,1])
plot(y[,1], coreg.y[rownames(y),"ENSG00000075420"],
     main=paste0("FNDC3B OLS model P: ", format(summary(m)$coefficients[2,4], scientific=T, digits=2)),
     xlab="-log10(gene P)",
     ylab="Coregulation protein coding")
abline(m, col="blue")


m <- lm(coreg.x[rownames(x),"ENSG00000137693"] ~ x[,1])
plot(x[,1], coreg.x[rownames(x),"ENSG00000137693"],
     main=paste0("YAP1 OLS model P: ", format(summary(m)$coefficients[2,4], scientific=T, digits=2)),
     xlab="-log10(gene P)",
     ylab="Coregulation skewness")
abline(m, col="blue")

m <- lm(coreg.y[rownames(y),"ENSG00000137693"] ~ y[,1])
plot(y[,1], coreg.y[rownames(y),"ENSG00000137693"],
     main=paste0("YAP1 OLS model P: ", format(summary(m)$coefficients[2,4], scientific=T, digits=2)),
     xlab="-log10(gene P)",
     ylab="Coregulation protein coding")
abline(m, col="blue")

