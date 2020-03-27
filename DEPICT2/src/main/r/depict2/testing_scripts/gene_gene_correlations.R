
t2 <- fread("/groups/umcg-wijmenga/scr02/depict2/Height10/Height_genePvaluesNullGwas.txt", data.table=F)
rownames(t2) <- t2[,1]
t2 <- t2[,-1]

collumns2 <- fread("/groups/umcg-wijmenga/scr02/depict2/Height10/Height_geneVariantCount.txt", data.table=F)

t2 <- t2[collumns2[collumns2[,2] > 0, 1], ]



c3 <- cor(t(t2[,1:10000]))

c4 <- cor(t(t2[,10001:20000]))

png(file="c3_c4_scatterplot.png", width=500, height=500)
plot(c3[upper.tri(c3, diag=F)], c4[upper.tri(c4, diag=F)])
dev.off()


pdf(file="c3_hist.pdf")
hist(c3[upper.tri(c3, diag=F)], breaks=100)
dev.off()

pdf(file="c4_hist.pdf")
hist(c4[upper.tri(c4, diag=F)], breaks=100)
dev.off()



gene.pos <- fread("/groups/umcg-wijmenga/tmp04/umcg-obbakker/data/reference/ensembl/genes_94/ensembl_gene_position_export_unique.sorted.bed", data.table=F)
rownames(gene.pos) <- gene.pos$V5

ol <- intersect(rownames(c3), rownames(gene.pos))
gene.pos2 <- gene.pos[ol,]


per.chrom.dist <- sapply(1:22, function(chr){
	cur.pos <- gene.pos2[gene.pos2$V1 == chr,2]
	names(cur.pos) <- gene.pos2[gene.pos2$V1 == chr,4]
	return(outer(cur.pos, cur.pos, "-"))
})

per.chrom.plots <- sapply(1:22, function(chr){
	cur.pos <- gene.pos2[gene.pos2$V1 == chr,2]
	names(cur.pos) <- gene.pos2[gene.pos2$V1 == chr,4]

	chr.dist <- abs(outer(cur.pos, cur.pos, "-"))

	chr.cors <- c3[names(cur.pos), names(cur.pos)]

	selector <- upper.tri(chr.dist, diag=F)
	cur.cor <- cor.test(chr.cors[selector], chr.dist[selector])

	png(width=500, height=500, file=paste0("chr_", chr, "_abs_distance_vs_correlation.png"))

	plot(chr.cors[selector], chr.dist[selector], main=paste0("chr", chr, " r=", cur.cor$estimate, " p=", cur.cor$p.value))

	dev.off()

})

per.chrom.plots <- sapply(1:22, function(chr){
	cur.pos <- gene.pos2[gene.pos2$V1 == chr,2]
	names(cur.pos) <- gene.pos2[gene.pos2$V1 == chr,4]

	chr.dist <- abs(outer(cur.pos, cur.pos, "-"))

	chr.cors <- c3[names(cur.pos), names(cur.pos)]


	selector <- upper.tri(chr.dist, diag=F)

	chr.cors[abs(chr.cors) < 0.1] <- NA


	x <- chr.cors[selector]

	y <- chr.dist[selector]
	y <- y[!is.na(x)]
	x <- x[!is.na(x)]
	cur.cor <- cor.test(y, x)

	png(width=500, height=500, file=paste0("chr_", chr, "_abs_distance_vs_correlation_noRandom.png"))

	plot(x,y, main=paste0("chr", chr, " r=", cur.cor$estimate, " p=", cur.cor$p.value))

	dev.off()

})

