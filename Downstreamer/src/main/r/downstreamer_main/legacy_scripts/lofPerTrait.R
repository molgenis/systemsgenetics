
source(paste0("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\Downstreamer\\src\\main\\r\\downstreamer_main/downstreamer_functions.r"))
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

traits <- read.delim("paperPheno.txt", stringsAsFactors = F)
ensg <- read.delim("ensgR75_protein_coding.txt", stringsAsFactors = F)



genePrioritizations <- lapply(paste0("final_paper/",traits$EnrichmentExcel), function(file){
  read.depict2(file)$GenePrioritization
})
names(genePrioritizations) <- traits$Name
str(genePrioritizations)

str(genePrioritizations[["height"]])

genePrioritizationsZscores <- lapply(genePrioritizations, function(x){
  return(x$Enrichment.Z.score[ match(genePrioritizations[["Height"]]$Gene.ID, x$Gene.ID)])
})
str(genePrioritizationsZscores)

genePrioritizationsBonfSig <- lapply(genePrioritizations, function(x){
  return(x$Bonferroni.significant[ match(genePrioritizations[["Height"]]$Gene.ID, x$Gene.ID)])
})
str(genePrioritizationsBonfSig)


genePascalPvalues <- lapply(genePrioritizations, function(x){
  return(x$GWAS.gene.P.value[ match(genePrioritizations[["Height"]]$Gene.ID, x$Gene.ID)])
})
str(genePascalPvalues)

pascalGeneP <- do.call(cbind, genePascalPvalues)
row.names(pascalGeneP) <- genePrioritizations[["Height"]]$Gene.ID

names(genePascalPvalues)[!names(genePascalPvalues) %in% colnames(pascalGeneP)]

dsGeneZ <- do.call(cbind, genePrioritizationsZscores)
row.names(dsGeneZ) <- genePrioritizations[["Height"]]$Gene.ID

dsGeneZ.max <- as.data.frame(dsGeneZ)
dsGeneZ.max$Max <- apply(dsGeneZ, 1, max, na.rm = T) 

pascalGeneP.max <- as.data.frame(pascalGeneP)
pascalGeneP.max$Max <- apply(pascalGeneP, 1, max, na.rm = T) 
pascalGeneP.max$Max[is.infinite(pascalGeneP.max$Max)] <- 0


library(heatmap3)
#heatmap3(dsGeneZ, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)

library(readr)
pli <- as.data.frame(read_delim("./gnomad.v2.1.1.lof_metrics.by_gene.txt.gz", delim = "\t", quote = ""))
row.names(pli) <- pli$gene_id
pli$lof_mis_z <- pmax(pli$lof_z, pli$mis_z, na.rm = T)
pli2 <- pli[row.names(dsGeneZ),]
str(pli2)


lofZCor <- cor(dsGeneZ.max, pli2[,"lof_z", drop= F], use = "pairwise.complete.obs")
str(lofZCor)
par(mar = c(15,5,3,1))
barplot(sort(lofZCor[,1]), las = 2, main = "Correlation DS-zscore vs gnomad mis_z")



misZCor <- cor(dsGeneZ.max, pli2[,"mis_z", drop= F], use = "pairwise.complete.obs")

LofMisZCor <- cor(dsGeneZ.max, pli2[,"lof_mis_z", drop= F], use = "pairwise.complete.obs")


LofMisZCorGeneP <- cor(pascalGeneP.max, pli2[,"lof_mis_z", drop= F], use = "pairwise.complete.obs")


synZCor <- cor(dsGeneZ.max, pli2[,"syn_z", drop= F], use = "pairwise.complete.obs")
par(mar = c(15,5,3,1))
barplot(sort(synZCor[,1]), las = 2, main = "Correlation DS-zscore vs gnomad syn_z")

barplot(sort(LofMisZCorGeneP[,1]), las = 2, main = "Correlation pascal -log10(p) vs gnomad max(lof_z, mis_z)")


barplot(sort(LofMisZCor[,1]), las = 2, main = "Correlation DS-zscore vs gnomad max(lof_z, mis_z)")

plot(LofMisZCorGeneP[,1], LofMisZCor[rownames(LofMisZCorGeneP),1], ylab = "Correlation DS-zscore vs gnomad max(lof_z, mis_z)", xlab = "Correlation pascal -log10(p) vs gnomad max(lof_z, mis_z)")
cor.test(LofMisZCorGeneP[,1], LofMisZCor[rownames(LofMisZCorGeneP),1])

plot(lofZCor, misZCor)
cor.test(lofZCor, misZCor)


plot(lofZCor, synZCor)
cor.test(lofZCor, synZCor)
