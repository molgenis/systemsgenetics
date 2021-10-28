library(data.table)

simple.qq.plot <- function (observedPValues) {
  plot(-log10(1:length(observedPValues)/length(observedPValues)), 
       -log10(sort(observedPValues)))
  abline(0, 1, col = "red")
}

genes <- fread("../../../ensgR75_protein_coding.txt", data.table=F)
expr <- fread("zcat ../exp.tmm.log2.covCorrect.freeze2.txt.gz")

expr <- as.data.frame(expr)
rownames(expr) <- expr[,1]
expr <- expr[,-1]


sample.cor <- cor(expr)

eigenvectors.bios    <- read.table("../bios_559_eigenvectors_protein_coding.txt")
coreg                <- cor(t(scale(eigenvectors.bios[, 1:100])))
coreg.559            <- cor(t(scale(eigenvectors.bios[, 1:559])))

#coreg.559           <- fread("../bios_559_eigenvectors_protein_coding_downstreamer_calculated.txt", data.table=F)
#rownames(coreg.559) <- coreg.559[,1]
#coreg.559           <- coreg.559[,-1]


gene      <- 3640
gene.n    <- "ENSG00000109320"
gm.output <- read.table("output/ENSG00000109320.assoc.txt", stringsAsFactors=F, row.names=1, header=T)
gm.output[gene, "beta"] <- 0

gene1.cor   <- cor(as.numeric(sc.expr[gene,]), t(sc.expr))
gene1.cor[gene1.cor == 1] <- 0

gene1.coreg <- coreg[,gene.n]
gene1.coreg[gene1.coreg==1] <- 0

gene1.coreg.559 <- coreg.559[,gene.n]
gene1.coreg.559[gene1.coreg.559 ==1] <- 0
names(gene1.coreg.559) <- rownames(coreg.559)

pdf(width=9, height=9, file=paste0("plots_gene_", gene, ".pdf"))
plot(gm.output$beta, gene1.cor , xlab="LMM beta", ylab="Pearson cor", main=paste0("BIOS coregulation for gene ", gene))
abline(a=0, b=1)


simple.qq.plot(gm.output$p)


ol <- intersect(rownames(gm.output), rownames(coreg))
plot(gm.output[ol,]$beta, gene1.coreg[ol] , xlab="LMM beta", ylab="Coregulation R 100pcs", main=paste0("BIOS coregulation for gene ", gene))
abline(a=0, b=1)


ol <- intersect(rownames(gm.output), rownames(coreg.559))
plot(gm.output[ol,]$beta, gene1.coreg.559[ol] , xlab="LMM beta", ylab="Coregulation R 559", main=paste0("BIOS coregulation for gene ", gene))
abline(a=0, b=1)


ol <- intersect(colnames(gene1.cor), rownames(coreg))
plot(gene1.cor[,ol], gene1.coreg[ol] , xlab="Pearson cor", ylab="Coregulation R 100pcs", main=paste0("BIOS coregulation for gene ", gene))
abline(a=0, b=1)

ol <- intersect(colnames(gene1.cor), rownames(coreg.559))
plot(gene1.cor[,ol], gene1.coreg.559[ol] , xlab="Pearson cor", ylab="Coregulation R 559pcs", main=paste0("BIOS coregulation for gene ", gene))
abline(a=0, b=1)


ol <- intersect(colnames(coreg.559), rownames(coreg))
plot(gene1.coreg.559[ol], gene1.coreg[ol] , xlab="Coregulation R 559pcs", ylab="Coregulation R 100pcs", main=paste0("BIOS coregulation for gene ", gene))
abline(a=0, b=1)


hist(gene1.cor, main="Pearson R", breaks=100, freq=F)
lines(density(rnorm(10000,mean=mean(gene1.cor), sd=sd(gene1.cor))), col="blue", lwd=2)


hist(gm.output$beta, main="LMM beta", breaks=100, freq=F)
lines(density(rnorm(10000,mean=mean(gm.output$beta), sd=sd(gm.output$beta))), col="blue", lwd=2)


hist(gene1.coreg, main="Coregulation R 100pcs", breaks=100, freq=F)
lines(density(rnorm(10000,mean=mean(gene1.coreg), sd=sd(gene1.coreg))), col="blue", lwd=2)


hist(gene1.coreg.559, main="Coregulation 559pcs", breaks=100, freq=F)
lines(density(rnorm(10000,mean=mean(gene1.coreg.559), sd=sd(gene1.coreg.559))), col="blue", lwd=2)

dev.off()