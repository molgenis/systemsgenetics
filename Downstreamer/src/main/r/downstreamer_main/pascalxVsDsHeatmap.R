load(file = "hpoEnrichmentSession4b.RData")


str(recount3Zscores)

pascalxZscores2 <- pascalxZscores
colnames(pascalxZscores2) <- gwas$Name[match(colnames(pascalxZscores2), gwas$machine_friendly_id)]


pascalDist <- as.dist(1 - cor(pascalxZscores2, use = "pa"))
pascalDist[is.nan(pascalDist)] <- 1
pascalDist[is.na(pascalDist)] <- 1
traitClustPascalX <- hclust(pascalDist, method = "ward.D2")



str(traitClust)



metaZscores2 <- metaZscores

colnames(metaZscores2) <- gwas$Name[match(colnames(metaZscores2), gwas$machine_friendly_id)]
str(recount3Zscores)
metaZDist <-  as.dist(1 - cor(metaZscores2, use = "pa"))
metaZDist[is.nan(metaZDist)] <- 1
metaZDist[is.na(metaZDist)] <- 1
traitClustMetaZ <- hclust(metaZDist, method = "ward.D2")
plot(traitClustMetaZ)


corMatrixMetaZ <- cor(metaZscores2, use = "pa")
corMatrixPascalX <- cor(pascalxZscores2, use = "pa")




corMatrixBoth <- corMatrixMetaZ[traitClustMetaZ$labels[traitClustMetaZ$order],traitClustMetaZ$labels[traitClustMetaZ$order]]
diag(corMatrixBoth) <- 0
corMatrixBoth[lower.tri(corMatrixBoth)] <- corMatrixPascalX[traitClustMetaZ$labels[traitClustMetaZ$order],traitClustMetaZ$labels[traitClustMetaZ$order]][lower.tri(corMatrixPascalX)]
corMatrixBoth[is.nan(corMatrixBoth)] <- 0
corMatrixBoth[is.na(corMatrixBoth)] <- 0


colHeatmap <- rev(c(colorRampPalette(c("#f03b20", "#feb24c", "#ffeda0"))(99), "white", colorRampPalette(c("#e0ecf4", "#9ebcda", "#8856a7"))(99)))
colBreaks <- c(seq(-1,-0.1,length.out= 100), seq(0.1,1,length.out= 100))


library(pheatmap)
pheatmap(corMatrixBoth, scale = "none", cluster_rows  = traitClustMetaZ, cluster_cols = traitClustMetaZ, col = colHeatmap, breaks = colBreaks)
pheatmap(corMatrixBoth, scale = "none", cluster_rows  = F, cluster_cols = F, col = colHeatmap, breaks = colBreaks)

pheatmap(corMatrixPascalX, scale = "none", cluster_rows  = traitClustMetaZ, cluster_cols = traitClustMetaZ, col = colHeatmap, breaks = colBreaks)
pheatmap(corMatrixMetaZ, scale = "none", cluster_rows  = traitClustMetaZ, cluster_cols = traitClustMetaZ, col = colHeatmap, breaks = colBreaks)

range(corMatrixPascalX)
range(corMatrixMetaZ, na.rm = T)


range(corMatrixBoth[lower.tri(corMatrixBoth)])
