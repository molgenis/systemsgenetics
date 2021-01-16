
source(paste0("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\Downstreamer\\src\\main\\r\\downstreamer_main/downstreamer_functions.r"))
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

mapping <- read.delim("ukbBurden/mapping.txt", stringsAsFactors = F)
ensg <- read.delim("ensgR75_protein_coding.txt", stringsAsFactors = F)


mapping$excel <- paste0("final_paper/", mapping$machine_friendly_id, "_enrichtments.xlsx")
mapping$excel[mapping$pheno != ""] <- paste0("final_paper/", mapping$machine_friendly_id, "_enrichtments_", mapping$pheno ,".xlsx")[mapping$pheno != ""]

mapping$trait2 <- mapping$trait
mapping$trait2[mapping$pheno != ""] <- mapping$pheno[mapping$pheno != ""]
mapping$trait2 <- gsub(".txt","",mapping$trait2)
mapping$trait2 <- gsub("_2019_ukbio_icbp_hg19","",mapping$trait2)
mapping$trait2 <- gsub("_2018_29507422_hg19","",mapping$trait2)

dsFiles <- unique(mapping$excel)
names(dsFiles) <- mapping$trait2[match(dsFiles, mapping$excel)]

genePrioritizations <- lapply(unique(dsFiles), function(file){
  read.depict2(file)$GenePrioritization
})
names(genePrioritizations) <- names(dsFiles)
str(genePrioritizations)

str(genePrioritizations[["height"]])

genePrioritizationsZscores <- lapply(genePrioritizations, function(x){
  return(x$Enrichment.Z.score[ match(genePrioritizations[["height"]]$Gene.ID, x$Gene.ID)])
})
str(genePrioritizationsZscores)

genePascalPvalues <- lapply(genePrioritizations, function(x){
  return(x$GWAS.gene.P.value[ match(genePrioritizations[["height"]]$Gene.ID, x$Gene.ID)])
})
str(genePascalPvalues)

pascalGeneP <- do.call(cbind, genePascalPvalues)
row.names(pascalGeneP) <- genePrioritizations[["height"]]$Gene.ID

dsGeneZ <- do.call(cbind, genePrioritizationsZscores)
row.names(dsGeneZ) <- genePrioritizations[["height"]]$Gene.ID


ukbCodes <- unique(mapping$UKB.phenotype.code)
names(ukbCodes) <- mapping$trait2[match(ukbCodes, mapping$UKB.phenotype.code)]

ukbCoding <- lapply(ukbCodes, function(code){
  read.delim(paste0("ukbBurden/allCoding/",code,".txt"), stringsAsFactors = F)
})

ukbLof<- lapply(ukbCodes, function(code){
  read.delim(paste0("ukbBurden/lof/",code,".txt"), stringsAsFactors = F)
})

ukbGenesCoding <- ukbCoding[[1]]$GENE

ukbGenesLof <- ukbLof[[1]]$GENE

all(sapply(ukbCoding, function(x){
  length(x$GENE) == length(ukbGenesCoding) & all(x$GENE==ukbGenesCoding)
}))
all(sapply(ukbLof, function(x){
  length(x$GENE) == length(ukbGenesLof) & all(x$GENE==ukbGenesLof)
}))

ukbCoding2 <- lapply(ukbCoding, function(x){
  pvalues <- x$P_BOLT_LMM_INF
})

ukbLof2 <- lapply(ukbLof, function(x){
  pvalues <- x$P_BOLT_LMM_INF
})

ukbCoding3 <- do.call(cbind, ukbCoding2)
row.names(ukbCoding3) <- ukbGenesCoding

ukbLof3 <- do.call(cbind, ukbLof2)
row.names(ukbLof3) <- ukbGenesLof

ukbCoding4 <- ukbCoding3[row.names(ukbCoding3) %in% ensg$Associated.Gene.Name,]
row.names(ukbCoding4) <- ensg$Ensembl.Gene.ID[match(row.names(ukbCoding4), ensg$Associated.Gene.Name)]

ukbLof4 <- ukbLof3[row.names(ukbLof3) %in% ensg$Associated.Gene.Name,]
row.names(ukbLof4) <- ensg$Ensembl.Gene.ID[match(row.names(ukbLof4), ensg$Associated.Gene.Name)]


str(ukbCoding4)


#pascalGeneP[is.na(pascalGeneP)] <- 1

ukbCodingCor <- cor(-log10(ukbCoding4))
ukbLofCor <- cor(-log10(ukbLof4))
dsGeneCor <- cor(dsGeneZ)
pascalGenePCor <- cor(-log10(pascalGeneP),use="pairwise.complete.obs")


pascalGenePNaCor <- cor(is.na(pascalGeneP))


corToZ <- function(r, df){
  t = sqrt(df) * abs(r) / sqrt(1 - (r*r))
  p = 2 * min(pt(t, df), pt(t, df, lower.tail=FALSE))
  z = qnorm((p/2))
  if(r > 0){
    z <- z * -1
  }
  if(z < -30){
    z <- -30
  } else if(z> 30){
    z <- 30
  }
  return(z)
}

ukbCodingCorZ <- apply(ukbCodingCor, 1:2, corToZ, df = nrow(ukbCoding4)-2)
ukbLofCorZ <- apply(ukbLofCor, 1:2, corToZ, df = nrow(ukbCoding4)-2)
dsGeneCorZ <- apply(dsGeneCor, 1:2, corToZ, df = nrow(dsGeneZ)-2)
#pascalGenePCorZ <- apply(pascalGenePCor, 1:2, corToZ, df = nrow(pascalGeneP)-2)


str(ukbCodingCor)
str(ukbCodingCorZ)


library(heatmap3)
heatmap3(ukbCodingCorZ, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)


heatmap3(ukbLofCorZ, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)


heatmap3(dsGeneCor, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)

heatmap3(pascalGenePCor, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)

heatmap3(pascalGenePNaCor, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)


geneOverlapCoding <- intersect(row.names(dsGeneZ), row.names(ukbCoding4))



dsGeneZ_c <- dsGeneZ[match(geneOverlapCoding, row.names(dsGeneZ)),]
pascalGeneP_c <- pascalGeneP[match(geneOverlapCoding, row.names(pascalGeneP)),]
ukbCoding4_c <- ukbCoding4[match(geneOverlapCoding, row.names(ukbCoding4)),]


dsVsCodingCor <- cor(dsGeneZ_c, -log10(ukbCoding4_c))
dsVsCodingCorZ <- apply(dsVsCodingCor, 1:2, corToZ, df = nrow(dsGeneZ_c)-2)


heatmap3(dsVsCodingCorZ, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)


geneOverlapLof <- intersect(row.names(dsGeneZ), row.names(ukbLof4))

dsGeneZ_l <- dsGeneZ[match(geneOverlapLof, row.names(dsGeneZ)),]
pascalGeneP_l <- pascalGeneP[match(geneOverlapLof, row.names(pascalGeneP)),]
ukbLof4_l <- ukbLof4[match(geneOverlapLof, row.names(ukbLof4)),]


dsVsLofCor <- cor(dsGeneZ_l, -log10(ukbLof4_l))
dsVsLofCorZ <- apply(dsVsLofCor, 1:2, corToZ, df = nrow(dsGeneZ_l)-2)


heatmap3(dsVsLofCorZ, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T, margins = c(10,10), Rowv = NA, Colv = NA)


plot(-log10(ukbLof4_l[,"NEU"]), dsGeneZ_l[,"NEU"])
cor.test(-log10(ukbLof4_l[,"NEU"]), dsGeneZ_l[,"NEU"])



load("gene_zcore_gene_prio_pearson_correlations.Rdata")

str(cor.m)
colnames(cor.m)


cor.m[,"height_2018_30124842"]

heatmap3(cor.m, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)




pascalVsCodingCor <- cor(-log10(pascalGeneP_c), -log10(ukbCoding4_c), use = "pairwise.complete.obs")
heatmap3(pascalVsCodingCor, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)

heatmap3(pascalVsCodingCor, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T, Rowv = NA, Colv= NA)





pascalVsCodingCor <- cor(-log10(pascalGeneP_c), -log10(ukbCoding4_c), use = "pairwise.complete.obs")
heatmap3(pascalVsCodingCor, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)

heatmap3(pascalVsCodingCor, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T, Rowv = NA, Colv= NA)



pascalVsLofCor <- cor(-log10(pascalGeneP_l), -log10(ukbLof4_l), use = "pairwise.complete.obs")
heatmap3(pascalVsCodingCor, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T)

heatmap3(pascalVsCodingCor, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T, Rowv = NA, Colv= NA)
