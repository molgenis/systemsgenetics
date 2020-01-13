load("~/Documents/data/gs_tcells/gene_expression/rpkm_Tcell_123samples.Rdata")

# Remove any genes with > 80% zeroes
tmp <- rpkm[rowSums(rpkm == 0) < 25,]
tmp <- na.omit(tmp)
dim(tmp)

write.table(tmp, quote=F, file="~/Documents/data/gs_tcells/gs_tcells_rpkm_80perc_zero_filtered.tsv", sep="\t")









# pc <- prcomp(t(tmp))


batch <- sapply(strsplit(rownames(tmp), "__TCC"), function(x){return (x[1])})

pca.plot(pc, fill=batch, size=4)
