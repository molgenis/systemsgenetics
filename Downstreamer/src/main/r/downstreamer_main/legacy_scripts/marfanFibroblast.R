setwd("C:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

library(GEOquery)
GSE8759 <- getGEO("GSE8759", GSEMatrix=TRUE, AnnotGPL=TRUE)

show(GSE8759)


if (length(GSE8759) > 1) idx <- grep("GPL538", attr(gset, "names")) else idx <- 1
GSE8759 <- GSE8759[[idx]]

str(GSE8759)


experimentData(GSE8759)

pheno <- pData(phenoData(GSE8759))

View(featureData(GSE8759)@data)

exp <- exprs(GSE8759[[1]])

dim(exp)
str(exp)
str(pheno)

pheno$sample <- as.factor(sapply(strsplit(as.character(pheno$title), " "), function(x){return(x[4])}))

expAgg <- t(apply(exp, 1, function(x){aggregate(x, by = list(pheno$sample), FUN = mean)$x}))
str(expAgg)

phenoAgg <- pheno[match(levels(pheno$sample), pheno$sample),]
rownames(expAgg) <- rownames(exp)
dimnames(expAgg)[[2]] <- levels(pheno$sample)

status <- as.factor(phenoAgg$`disease status:ch1`)




pcaRes <- prcomp(expAgg, center = T, scale. = F)
plot(pcaRes$rotation[,1], pcaRes$rotation[,2], col = status)

pairs(pcaRes$rotation[,1:5], col = status)

res <- apply(expAgg, 1, function(x){t.test(x ~ status)$p.value})


res2 <- data.frame(gene = featureData(GSE8759)@data$`Gene symbol`, pvalue = res, geneId = featureData(GSE8759)@data$`Gene ID`)

map <- read.delim("ensgNcbiId.txt", stringsAsFactors = F)

res2$ensg <- map$Gene.stable.ID[match(res2$geneId, map$NCBI.gene.ID)]
View(res2)

write.table(res2, file = "MarfanExpression/result.txt", quote = F, sep = "\t",row.names = F)



heightCoreGene <- read.delim("height_2018_30124842_hg19_48/height_2018_30124842_hg19_intermediates_Coregulation_Enrichment_zscoreExHla.txt")

heightCoreGene$marfan = res2$pvalue[match(heightCoreGene$X., res2$ensg)]

heightCoreGene2 <- heightCoreGene[!is.na(heightCoreGene$marfan),]

cor.test(heightCoreGene2$p, -log10(heightCoreGene2$marfan))
