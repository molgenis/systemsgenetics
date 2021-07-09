library(data.table)
library(matrixStats)
library(DescTools)

# Read expression data from bios
exp.bios           <- fread("zcat exp.tmm.log2.covCorrect.freeze2.txt.gz", data.table=F)
rownames(exp.bios) <- exp.bios[,1]
exp.bios           <- exp.bios[,-1]
cn                 <- colnames(fread("zcat exp.tmm.log2.covCorrect.freeze2.txt.gz | head -n 1", data.table=F))
colnames(exp.bios) <- cn
exp.bios           <- as.matrix(exp.bios)

# Ensembl genes 
ensembl <- read.table("/groups/umcg-wijmenga/tmp04/projects/depict2/depict2_bundle/reference_datasets/human_b37/ensgR75_protein_coding.txt", row.names=1, header=T)

# Filter any lowly expressed and non varying genes
vars       <- rowVars(exp.bios)


# Histograms used to eyeball thresholds
pdf(width=5, height=5, file="bios_vars.pdf")
hist(vars, breaks=100)
hist(vars, breaks=1000, xlim=c(0,1))
hist(exp.bios[vars < 0.1,], breaks=1000)
hist(exp.bios[,1], breaks=1000)
hist(exp.bios[,1], breaks=1000, xlim=c(0,5))
dev.off()

exp.filter        <- rowSums(exp.bios > 0.5) > 0.75*ncol(exp.bios)
var.filter        <- vars > 0.01
filter            <- (var.filter + exp.filter) == 2
exp.bios.filtered <- exp.bios[filter,]
exp.bios.filtered <- exp.bios.filtered[intersect(rownames(exp.bios.filtered), rownames(ensembl)),]

# Calculate eigenvectors
cor.bios <- cor(t(exp.bios.filtered))
save(cor.bios, file="bios_protein_coding_gene_correlations_for_downstreamer")

eigenvectors.bios <- eigen(cor.bios)
save(eigenvectors.bios, file="bios_eigenvectors_for_downstreamer")

# Determine number of pcs for 75% variance
perc.var <- eigenvectors.bios$values / sum(eigenvectors.bios$values)

# first > 0/75 is at 559 components
sum(perc.var[1:559])

# Prepare eigenvectors
coreg.vecs           <- eigenvectors.bios$vectors[,1:559]
rownames(coreg.vecs) <- rownames(cor.bios)
write.table(coreg.vecs, file="bios_559_eigenvectors_protein_coding.txt", sep="\t", quote=F)
coreg.vecs           <- scale(coreg.vecs)

# Calculate coregulation zscores (for validation, way too slow for full comparision)
coreg <- apply(coreg.vecs[1:10,], 1, function(x){
	apply(coreg.vecs, 1, function(y)  {
		cur.cor <- cor.test(x, y)
		z       <- qnorm(cur.cor$p.value/2, lower.tail=F)
		z       <- sign(cur.cor$estimate) * z
		return(z)
	})
})

# Set diagonal to zero
diag(coreg) <- 0
save(coreg, file="bios_coregulation_zscores_protein_coding_for_downstreamer")

# Compare with downstreamer Z-scores
ds.coreg           <- fread("bios_559_eigenvectors_protein_coding_downstreamer_calculated.txt", data.table=F)
rownames(ds.coreg) <- ds.coreg[,1]
ds.coreg           <- ds.coreg[,-1]

ol <- intersect(colnames(coreg), colnames(ds.coreg))

png(width=500, height=500, file="bios_coreg_comparision_ds_vs_r.png")
plot(ds.coreg[rownames(coreg), ol[1]], coreg[,ol[1]])
dev.off()