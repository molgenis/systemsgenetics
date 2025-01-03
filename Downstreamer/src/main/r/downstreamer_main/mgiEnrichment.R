

#remoter::server(port = 55556, sync = T)
#remoter::server(port = 54001, sync = T)

#ml R/4.2.1-foss-2022a-bare

remoter::client("localhost", port = 54001)#55556 55501   54104

setwd("/groups/umcg-fg/tmp02/projects/downstreamer/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")


library(readxl)


read.depict2 <- function(path, potential_traits=NULL) {
  if (is.null(potential_traits)) {
    potential_traits <- excel_sheets(path)
    potential_traits <- potential_traits[grep("Overview", potential_traits, invert=T)]
  }
  
  output <- list()
  for (sheet in potential_traits) {
    tmp <- tryCatch({data.frame(read_excel(path, sheet=sheet, col_types ="guess", trim_ws = T), stringsAsFactors=F)},
                    error=function(a){return(NA)},
                    warn=function(a){return(NA)})
    
    
    for (i in 1:ncol(tmp)) {
      if (class(tmp[,i]) == "character"){
        tmp[,i] <- type.convert(tmp[,i], as.is=T)
        
      }
    }
    
    rownames(tmp) <- tmp[,1]
    output[[sheet]] <- tmp
  }
  
  return(output)
}


geneInfo <- read.delim("depict2_bundle/reference_datasets/human_b37/genes_Ensembl94_protein_coding.txt", header = F)
str(geneInfo)

library(readr)


#table_tmp <- read_delim("depict2_bundle/reference_datasets/human_b37/hpo/phenotype_to_genes_V1268_OMIMandORPHA.txt_matrix.txt.gz", delim = "\t", quote = "", col_types = colTypes)

#change X1 in case of specified header for row names
colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/pathway_databases/HPO/2023_06_17/HPO_2023_06_17.txt.gz", delim = "\t", quote = "", col_types = colTypes)
hpoMatrix <- as.matrix(table_tmp[,-1])
rownames(hpoMatrix) <- table_tmp[,1][[1]]
rm(table_tmp) 

#hpoAnnotations <- as.matrix(read.delim("depict2_bundle/reference_datasets/human_b37/hpo/phenotype_to_genes_V1268_OMIMandORPHA.txt_matrix.colAnnotations.txt", row.names = 1))
hpoAnnotations <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/genenetwork/pathway_databases/HPO/2023_06_17/hpoToName.txt", row.names = 1))
hpoAnnotations <- hpoAnnotations[,1]


hpoPerGene <- apply(hpoMatrix, 1, function(x){sum(x>0)})
hpoMatrix <- hpoMatrix[hpoPerGene >= 1,]

genePerHpo <- apply(hpoMatrix, 2, function(x){sum(x>0)})
hpoMatrix <- hpoMatrix[,genePerHpo >= 10]
dim(hpoMatrix)



gwas <- read.delim("traits.txt", header = T)
traits <- gwas[,1]

traits2 <- c("Height_2022", "schizophrenia_2018_29483656_hg19", "IBD_deLange2017")

metaZscores <- matrix(NA, nrow = nrow(geneInfo), ncol = length(traits), dimnames = list(geneInfo$V1, traits))
recount3Zscores <- matrix(NA, nrow = nrow(geneInfo), ncol = length(traits), dimnames = list(geneInfo$V1, traits))
pascalxZscores <- matrix(NA, nrow = nrow(geneInfo), ncol = length(traits), dimnames = list(geneInfo$V1, traits))


trait <- traits[2]
sink <- lapply(traits, function(trait){
  
    
  cat(trait,"\n")

  
  
  pascalxRes <- read.delim(paste0("/groups/umcg-fg/tmp02/projects/downstreamer/PascalX_bundle/results_25kb/",trait, ".txt"), row.names =1 )
  #str(pascalxRes)
  pxg <- intersect(geneInfo$V1, row.names(pascalxRes))
  range(pascalxRes[pxg,"pvalue"])
  range(-qnorm(pascalxRes[pxg,"pvalue"]))
  pascalxZscores[pxg,trait] <<- -qnorm(pascalxRes[pxg,"pvalue"]/2)
  
  
  res <- read.depict2(paste0("depict2_bundle/output/ds2_B_25k/", trait,"/",trait,"_keygenes3_noCovCor_enrichtments.xlsx"))
  
  
  tissues <- names(res)[names(res) != "Recount3"]
  
  
  
  genes <- sapply(tissues, function(t){
    res[[t]]$Gene.ID
  })
  if(is.list(genes)){
    genes <- unique(do.call(c, genes))
  }
  
  
  tissueZscores <- matrix(NA, nrow = length(genes), ncol = length(tissues), dimnames = list(genes, tissues))
  
  sapply(tissues, function(t){
    resT <- res[[t]]
    resTz <- resT$Enrichment.Z.score
    if(all(resTz == 0)){
      tissueZscores[resT$Gene.ID, t] <<- rep(NA, length(resT$Gene.ID))
    } else {
      tissueZscores[resT$Gene.ID, t] <<- resTz
    }
    
    return(1)
  })
  
  
  metaZ <- apply(tissueZscores, 1, function(z){
    z <- z[!is.na(z)]
    
    return(sum(z) / sqrt(length(z)))
    
  })
  
  metaZ <- metaZ[!is.na(metaZ)]
  
  
  metaZscores[names(metaZ), trait] <<- metaZ
  
  
  resRecount3 <- res[["Recount3"]]

  recount3Zscores[resRecount3$Gene.ID, trait] <<- resRecount3$Enrichment.Z.score
  
  
  
})

#save(recount3Zscores, geneInfo, pascalxZscores, metaZscores, traits, hpoMatrix, hpoAnnotations, file = "depict2_bundle/hpoEnrichmentSession4.RData")


setwd("/groups/umcg-fg/tmp02/projects/downstreamer/")


load("hpoEnrichmentSession4.RData")


trait <- traits[1]
trait <- "IBD_deLange2017"
trait <- "Schizophrenia_Pardinas2018"
trait <- "ADHD_Demontis2018"
trait <- "Height_2022"

library(parallel)

clust <- makeCluster(10)


clusterExport(clust, "recount3Zscores")
clusterExport(clust, "metaZscores")
clusterExport(clust, "pascalxZscores")
clusterExport(clust, "hpoMatrix")
clusterExport(clust, "hpoAnnotations")



hpoEnrichmentList <- parLapply(clust, traits, function(trait){
  
  metaZ <- metaZscores[, trait] 
  metaZ <- metaZ[!is.na(metaZ)]
  length(metaZ)
  
  sum(is.na(metaZscores[, trait] ))
  
  
  sigThreshold <- -qnorm(0.05/length(metaZ)/2)
  sig <- metaZ >= sigThreshold
  
  
  overlapGenes <- intersect(names(sig), row.names(hpoMatrix))
  
  sig <- sig[overlapGenes]
  if(sum(sig) <= 10){ 
      tissueMetaHpoEnrichment <- NA
    }else {
      
      
    hpoMatrixS <- hpoMatrix[overlapGenes,]
    hpoMatrixS <- hpoMatrixS[,apply(hpoMatrixS, 2, function(x){sum(x>0)})>=10]
    hpoMatrixS <- hpoMatrixS[,apply(hpoMatrixS, 2, function(x){sum(x<1)})>=10]
    
    x <- hpoMatrixS[,1]
    tissueMetaHpoEnrichment <- t(apply(hpoMatrixS, 2, function(x){
      tryCatch( 
        {  
          
          f <- fisher.test(table(sig, x))
          return(c(f$p.value, f$estimate, as.vector(table(sig, x))))
          
          
        },
        error=function(cond) {
          print(table(sig, x)) 
        }
      )
    }))
    
    tissueMetaHpoEnrichment <- as.data.frame(tissueMetaHpoEnrichment)
    colnames(tissueMetaHpoEnrichment)[1] <- "Fisher-Pvalue"
    colnames(tissueMetaHpoEnrichment)[3:6] <- c("noKeyNoHpo", "keyNoHpo", "noKeyHpo", "keyHpo")
    tissueMetaHpoEnrichment$descip <- hpoAnnotations[rownames(tissueMetaHpoEnrichment)]
    
  }
  
  
  
  pascalZ <- pascalxZscores[, trait] 
  pascalZ <- pascalZ[!is.na(pascalZ)]
  length(pascalZ)
  str(pascalxZscores[,trait])
  
  sigThreshold <- -qnorm(0.05/length(pascalZ)/2)
  sig <- pascalZ >= sigThreshold
  
  
  overlapGenes <- intersect(names(sig), row.names(hpoMatrix))
  
  sig <- sig[overlapGenes]
  if(sum(sig) <= 10){ 
    pascalXHpoEnrichment <- NA
  }else {
    
    
    hpoMatrixS <- hpoMatrix[overlapGenes,]
    hpoMatrixS <- hpoMatrixS[,apply(hpoMatrixS, 2, function(x){sum(x>0)})>=10]
    hpoMatrixS <- hpoMatrixS[,apply(hpoMatrixS, 2, function(x){sum(x<1)})>=10]
    
    x <- hpoMatrixS[,1]
    pascalXHpoEnrichment <- t(apply(hpoMatrixS, 2, function(x){
      tryCatch( 
        {  
          
          f <- fisher.test(table(sig, x))
          return(c(f$p.value, f$estimate, as.vector(table(sig, x))))
          
          
        },
        error=function(cond) {
          print(table(sig, x)) 
        }
      )
    }))
    
    pascalXHpoEnrichment <- as.data.frame(pascalXHpoEnrichment)
    colnames(pascalXHpoEnrichment)[1] <- "Fisher-Pvalue"
    colnames(pascalXHpoEnrichment)[3:6] <- c("noKeyNoHpo", "keyNoHpo", "noKeyHpo", "keyHpo")
    pascalXHpoEnrichment$descip <- hpoAnnotations[rownames(pascalXHpoEnrichment)]
    
  }
  
  
  
  
  

  
  resRecount3 <- recount3Zscores[, trait] 
  resRecount3 <- resRecount3[!is.na(resRecount3)]
  
  
  
  sigThreshold <- -qnorm(0.05/length(resRecount3)/2)
  sig <- resRecount3 >= sigThreshold
  
  
  overlapGenes <- intersect(names(sig), row.names(hpoMatrix))
  
  sig <- sig[overlapGenes]
  if(sum(sig) <= 10){ 
    recount3HpoEnrichment <- NA
  }else {
    
    
    hpoMatrixS <- hpoMatrix[overlapGenes,]
    hpoMatrixS <- hpoMatrixS[,apply(hpoMatrixS, 2, function(x){sum(x>0)})>=10]
    hpoMatrixS <- hpoMatrixS[,apply(hpoMatrixS, 2, function(x){sum(x<1)})>=10]
    
    x <- hpoMatrixS[,1]
    recount3HpoEnrichment <- t(apply(hpoMatrixS, 2, function(x){
      tryCatch( 
        {  
          
          f <- fisher.test(table(sig, x))
          return(c(f$p.value, f$estimate, as.vector(table(sig, x))))
          
          
        },
        error=function(cond) {
          print(table(sig, x)) 
        }
      )
    }))
    
    recount3HpoEnrichment <- as.data.frame(recount3HpoEnrichment)
    colnames(recount3HpoEnrichment)[1] <- "Fisher-Pvalue"
    colnames(recount3HpoEnrichment)[3:6] <- c("noKeyNoHpo", "keyNoHpo", "noKeyHpo", "keyHpo")
    recount3HpoEnrichment$descip <- hpoAnnotations[rownames(recount3HpoEnrichment)]
    
  }
  
  
  return(list(recount3HpoEnrichment, tissueMetaHpoEnrichment, pascalXHpoEnrichment))

}) 

#save(recount3Zscores, metaZscores, pascalxZscores, traits, geneInfo, hpoEnrichmentList, hpoAnnotations, file = "depict2_bundle/hpoEnrichmentSession4bNewHpo.RData")
load(file = "hpoEnrichmentSession4bNewHpo.RData")

names(hpoEnrichmentList) <- traits

View(hpoEnrichmentList$Height[[2]])
str(hpoEnrichmentList$IBD_deLange2017[[3]])

enrichmentCompare <- merge(hpoEnrichmentList$IBD_deLange2017[[2]], hpoEnrichmentList$IBD_deLange2017[[3]], by = "descip", all = T, suffixes = c ("-DS", "-PX"))
plot(-log10(enrichmentCompare$`Fisher-Pvalue-PX`), -log10(enrichmentCompare$`Fisher-Pvalue-DS`), xlab = "PascalX enrichment", ylab = "Downtreamer enrichment", main = "IBD")
abline(coef = c(0,1))
dev.off()
write.table(enrichmentCompare, file = "ibd.txt", sep = "\t", quote = F, col.names = NA)



enrichmentCompare <- merge(hpoEnrichmentList$Height_2022[[2]], hpoEnrichmentList$Height_2022[[3]], by = "descip", all = T, suffixes = c ("-DS", "-PX"))
plot(-log10(enrichmentCompare$`Fisher-Pvalue-PX`), -log10(enrichmentCompare$`Fisher-Pvalue-DS`), xlab = "PascalX enrichment", ylab = "Downtreamer enrichment", main = "Height")
abline(coef = c(0,1))
dev.off()
write.table(enrichmentCompare, file = "height.txt", sep = "\t", quote = F, col.names = NA)




enrichmentCompare <- merge(hpoEnrichmentList$schizophrenia_2018_29483656_hg19[[2]], hpoEnrichmentList$schizophrenia_2018_29483656_hg19[[3]], by = "descip", all = T, suffixes = c ("-DS", "-PX"))
plot(-log10(enrichmentCompare$`Fisher-Pvalue-PX`), -log10(enrichmentCompare$`Fisher-Pvalue-DS`), xlab = "PascalX enrichment", ylab = "Downtreamer enrichment", main = "Schizophrenia")
abline(coef = c(0,1))
dev.off()
write.table(enrichmentCompare, file = "schizo.txt", sep = "\t", quote = F, col.names = NA)


#hpoEnrichmentListold <- hpoEnrichmentList
#hpoEnrichmentListNew <- hpoEnrichmentList

#hpoEnrichmentList50kb <- hpoEnrichmentList
#hpoEnrichmentList25kb <- hpoEnrichmentList



str(hpoEnrichmentList25kb$IBD_deLange2017[[2]])
str(hpoEnrichmentList50kb$IBD_deLange2017[[2]])

heigthEnrichmentCompare <- merge(hpoEnrichmentListold$IBD_deLange2017[[2]], hpoEnrichmentListNew$IBD_deLange2017[[2]], by = "descip", all = T, suffixes = c ("fdr05", "fdr25"))


plot(-log10(heigthEnrichmentCompare$`Fisher-Pvaluefdr05`), -log10(heigthEnrichmentCompare$`Fisher-Pvaluefdr25`), xlab = "Eigen FDR 0.05", ylab = "Eigen FDR 0.25", main = "IBD")
abline(coef = c(0,1))
dev.off()


heigthEnrichmentCompare <- merge(hpoEnrichmentListold$Height_2022[[2]], hpoEnrichmentListNew$Height_2022[[2]], by = "descip", all = T, suffixes = c ("fdr05", "fdr25"))


plot(-log10(heigthEnrichmentCompare$`Fisher-Pvaluefdr05`), -log10(heigthEnrichmentCompare$`Fisher-Pvaluefdr25`), xlab = "Eigen FDR 0.05", ylab = "Eigen FDR 0.25", main = "Height")
abline(coef = c(0,1))
dev.off()


heigthEnrichmentCompare <- merge(hpoEnrichmentListold$schizophrenia_2018_29483656_hg19[[2]], hpoEnrichmentListNew$schizophrenia_2018_29483656_hg19[[2]], by = "descip", all = T, suffixes = c ("fdr05", "fdr25"))
str(heigthEnrichmentCompare)

plot(-log10(heigthEnrichmentCompare$`Fisher-Pvaluefdr05`), -log10(heigthEnrichmentCompare$`Fisher-Pvaluefdr25`), xlab = "Eigen FDR 0.05", ylab = "Eigen FDR 0.25", main = "Schizophrenia")
abline(coef = c(0,1))
dev.off()






write.table(heigthEnrichmentCompare, file = "tmp.txt", sep = "\t", quote = F, col.names = NA)


plot(heigthEnrichmentCompare$`odds ratio50kb`, heigthEnrichmentCompare$`odds ratio25kb`)
dev.off()


str(pascalxZscores)
pascalxSig <- apply(pascalxZscores, 2, function(x){
  t <- -qnorm((0.05/sum(!is.na(x)))/2)
  
  x[!is.na(x)] <- x[!is.na(x)] >= t
  return(x)
})

write.table(pascalxSig, file = gzfile("signficantPascalx.txt.gz") , sep = "\t", quote = F, col.names = NA, na = "")
write.table(pascalxZscores, file = gzfile("pascalxZscores.txt.gz") , sep = "\t", quote = F, col.names = NA, na = "")
hist(pascalxZscores)
barplot(apply(pascalxSig, 2, sum, na.rm=T))
sort(apply(pascalxSig, 2, sum, na.rm=T))
metaZscores <- metaZscores[apply(metaZscores, 1, function(x){!all(is.na(x))}),]
str(metaZscores)

x <- metaZscores[,2]
metaZscoresSig <- apply(metaZscores, 2, function(x){
  t <- -qnorm((0.05/sum(!is.na(x)))/2)
  
  x[!is.na(x)] <- x[!is.na(x)] >= t
  return(x)
  })

write.table(metaZscoresSig, file = gzfile("signficantMetaKeyGene.txt.gz") , sep = "\t", quote = F, col.names = NA, na = "")
write.table(metaZscores, file = gzfile("metaKeyGene.txt.gz") , sep = "\t", quote = F, col.names = NA, na = "")

  str(metaZscoresSig)


apply(metaZscoresSig, 2, sum, na.rm=T)
apply(metaZscoresSig, 2, function(x){sum(is.na(x))})

table(metaZscoresSig, useNA = "a")

sum(metaZscores[,"Height_2022"] >= 4.66, na.rm = T)





trait <- "Alzheimers_Jansen2019"
trait <- "IBD_deLange2017"
trait <- "Height_2022"
trait <- traits[1]

str(hpoEnrichmentList)

x <- lapply(traits, function(trait){
  
  print(trait)
  tissueMetaHpoEnrichment <- hpoEnrichmentList[[trait]][[2]]
  
  if(length(tissueMetaHpoEnrichment)>1){
  
  
  
  bonThres <- 0.05/nrow(tissueMetaHpoEnrichment)
  
  
  tissueMetaHpoEnrichment$log2Or <- log2(tissueMetaHpoEnrichment$`odds ratio`)
  
  tissueMetaHpoEnrichmentFilter <- tissueMetaHpoEnrichment[tissueMetaHpoEnrichment$`Fisher-Pvalue` <= bonThres & tissueMetaHpoEnrichment$log2Or >= 0,] 

    
    
  #print(tissueMetaHpoEnrichmentFilter[order(tissueMetaHpoEnrichmentFilter$`Fisher-Pvalue`, decreasing=F)[1:10],], width = 500)
  print(tissueMetaHpoEnrichmentFilter[order(abs(tissueMetaHpoEnrichmentFilter$`log2Or`), decreasing=T)[1:10],], width = 500)
  } else {
    print(" ")
  }
})

str(hpoAnnotations)


hpoEnrichOrMatrix <- matrix(NA, nrow = length(hpoAnnotations), ncol = length(traits), dimnames = list(hpoAnnotations, traits))

trait <- traits[1]

x <- lapply(traits, function(trait){
  
  print(trait)
  tissueMetaHpoEnrichment <- hpoEnrichmentList[[trait]][[2]]
  
  if(length(tissueMetaHpoEnrichment)>1){
    
    bonThres <- 0.05/nrow(tissueMetaHpoEnrichment)
    
    tissueMetaHpoEnrichment[tissueMetaHpoEnrichment$descip=="Abnormality of the thorax",]
    tissueMetaHpoEnrichment$log2Or <- log2(tissueMetaHpoEnrichment$`odds ratio`)
    #View(tissueMetaHpoEnrichment)
    tissueMetaHpoEnrichmentFilter <- tissueMetaHpoEnrichment[tissueMetaHpoEnrichment$`Fisher-Pvalue` <= bonThres  & tissueMetaHpoEnrichment$log2Or >= 0,] 
    #View(tissueMetaHpoEnrichmentFilter)
    hpoEnrichOrMatrix[tissueMetaHpoEnrichmentFilter$descip,trait] <<- tissueMetaHpoEnrichmentFilter$log2Or
    
  }
  
})

hpoEnrichOrMatrix2 <- hpoEnrichOrMatrix[apply(hpoEnrichOrMatrix, 1, function(x){any(!is.na(x))}), apply(hpoEnrichOrMatrix, 2, function(x){any(!is.na(x))})]

str(hpoEnrichOrMatrix2)


colnames(hpoEnrichOrMatrix2) <- gwas$Name[match(colnames(hpoEnrichOrMatrix2), gwas$machine_friendly_id)]

#colnames(hpoMatrix)[!colnames(hpoMatrix) %in% names(hpoAnnotations)]

all(colnames(hpoMatrix) %in% names(hpoAnnotations))
hpoMatrix2 <- hpoMatrix
colnames(hpoMatrix2) <- hpoAnnotations[colnames(hpoMatrix)]
hpoMatrix2 <- hpoMatrix2[,rownames(hpoEnrichOrMatrix2)]


#hpoEnrichOrMatrix2ForCluster <- hpoEnrichOrMatrix2

#hpoEnrichOrMatrix2ForCluster[is.na(hpoEnrichOrMatrix2ForCluster)] <- 0

#hist(hpoEnrichOrMatrix2ForCluster)
#dev.off()


metaZscores2 <- metaZscores
colnames(metaZscores2) <- gwas$Name[match(colnames(metaZscores2), gwas$machine_friendly_id)]


traitDist <- as.dist(1 - cor(metaZscores2[,colnames(hpoEnrichOrMatrix2)], use = "pa"))
traitDist[is.nan(traitDist)] <- 1
traitDist[is.na(traitDist)] <- 1
traitClust <- hclust(traitDist, method = "ward.D2")


plot(traitClust)
dev.off()

maxLog2Or <- max(abs(hpoEnrichOrMatrix2[!is.infinite(hpoEnrichOrMatrix2)]), na.rm = T)

hpoEnrichOrMatrix2[is.na(hpoEnrichOrMatrix2)] <- 0



#hpoEnrichOrMatrix2[is.infinite(hpoEnrichOrMatrix2) & hpoEnrichOrMatrix2 < 0] <- -maxLog2Or
hpoEnrichOrMatrix3 <- hpoEnrichOrMatrix2
hpoEnrichOrMatrix3[is.infinite(hpoEnrichOrMatrix2) & hpoEnrichOrMatrix2 > 0] <- maxLog2Or


#colHeatmap <- rev(c(colorRampPalette(c("#f03b20", "#feb24c", "#ffeda0"))(99), "white", colorRampPalette(c("#e0ecf4", "#9ebcda", "#8856a7"))(99)))
#colBreaks <- c(seq(-maxLog2Or,-2,length.out= 100), seq(2,maxLog2Or,length.out= 100))

colHeatmap <- rev(c(colorRampPalette(c("#f03b20", "#feb24c", "#ffeda0"))(99), "white"))
#colBreaks <- c(0,seq(2,maxLog2Or,length.out= 100))
colBreaks <- c(0,seq(0.1,maxLog2Or,length.out= 100))

str(hpoMatrix2)
hpoDist <- dist(t(hpoMatrix2), method = "manhattan")
hpoClust <- hclust(hpoDist, method =  "ward.D2")
plot(hpoClust)


write.table(hpoEnrichOrMatrix2[hpoClust$labels[hpoClust$order],traitClust$labels[traitClust$order]], file = "hpoErichmentHeatmap25kNoCovCor.txt", sep = "\t", quote = F, col.names = NA, na = "")


unique(gwas$class)

gwasHeatmapAnnotation <- data.frame(Class = gwas$class, row.names = gwas$Name)
str(gwasHeatmapAnnotation)

library(pheatmap)


#rpng(width = 600, height = 1000)
pdf("hpoErichmentHeatmapNewHpo.pdf", height = 180, width = 20)#65
pheatmap(hpoEnrichOrMatrix3, cluster_rows  = hpoClust, cluster_cols = traitClust, col = colHeatmap, breaks = colBreaks, scale = "none", cellwidth = 10, cellheight = 10, treeheight_row = 300, treeheight_col = 100, annotation_col = gwasHeatmapAnnotation, legend_breaks = seq(0,6,2))
dev.off()

hpoExample <- read.delim("HpoExample.txt", header = F)$V1


hpoDist <- dist(t(hpoMatrix2[,hpoExample]), method = "manhattan")
hpoClust <- hclust(hpoDist)#, method =  "ward.D2"
plot(hpoClust)


pdf("hpoErichmentHeatmapNewHpoSelection.pdf", height = 11, width = 20)#65
pheatmap(hpoEnrichOrMatrix3[hpoExample,], cluster_rows  = hpoClust, cluster_cols = traitClust, col = colHeatmap, breaks = colBreaks, scale = "none", cellwidth = 10, cellheight = 10, treeheight_row = 300, treeheight_col = 100, annotation_col = gwasHeatmapAnnotation, legend_breaks = seq(0,6,2))
dev.off()



class = "cardiovascular"
class = "brain"
for(class in unique(gwas$class)){
  
  print(class)
  
  slectedTraits <- colnames(hpoEnrichOrMatrix3)[colnames(hpoEnrichOrMatrix3) %in% gwas$Name[gwas$class==class]]
  
  if(length(slectedTraits) > 1){
    traitDistClass <- as.dist(1 - cor(metaZscores2[,slectedTraits], use = "pa"))
    traitDistClass[is.nan(traitDistClass)] <- 1
    traitDistClass[is.na(traitDistClass)] <- 1
    traitClustClass <- hclust(traitDistClass, method = "ward.D2")
    
  } else {
    traitClustClass = FALSE
  }
  #needs new clustering per class
  
  hpoEnrichOrMatrixClass <- hpoEnrichOrMatrix3[,slectedTraits, drop = F]
  hpoEnrichOrMatrixClass <- hpoEnrichOrMatrixClass[apply(hpoEnrichOrMatrixClass,1,function(x){any(x>0)}),,drop = F]

  hpoClustClass <- hclust(dist(t(hpoMatrix2[,rownames(hpoEnrichOrMatrixClass),drop =F]), method = "manhattan"), method =  "ward.D2")
  pdf(paste0("hpoErichmentHeatmap_",class,"newHpo.pdf"), height = 75, width = 12)
  pheatmap(hpoEnrichOrMatrixClass, cluster_rows  = hpoClustClass, cluster_cols = traitClustClass, col = colHeatmap, breaks = colBreaks, scale = "none", cellwidth = 10, cellheight = 10, treeheight_row = 300, treeheight_col = 100, legend_breaks = seq(2,6,2), main = class)
  dev.off()
  
}


sum(is.na(hpoEnrichOrMatrix2ForCluster))
range(hpoEnrichOrMatrix2ForCluster)

sum(is.na(hpoEnrichOrMatrix2))
sum( is.infinite (hpoEnrichOrMatrix2))






#  pdf("ibdEnrichment.pdf")
#  par(xpd = NA)
#  plot(log2(tissueMetaHpoEnrichment[,2]), -log10(tissueMetaHpoEnrichment[,1]), pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), xlab = "log2(odds ratio)", ylab = "-log10(Fisher pvalue)", main = "Schizo tissue specific meta analyzed enrichment")
#  abline(h = -log10(0.05/nrow(tissueMetaHpoEnrichment)), lwd = 2, col = "darkred")
#  dev.off()

#  tissueMetaHpoEnrichment[order(tissueMetaHpoEnrichment[,1], decreasing = F)[1:20],]
#  tissueMetaHpoEnrichment[order(tissueMetaHpoEnrichment[,2], decreasing = T)[1:20],]


#  library(pROC)
#  table(sig, hpoMatrixS[,"HP:0002733"])
#  r <- roc(response = sig, predictor = hpoMatrixS[,"HP:0002733"])
#  plot.roc(r)


#  library(beeswarm)
#  beeswarm(metaZ[overlapGenes][hpoMatrixS[,"HP:0002583"]==1])
#  abline(h=sigThreshold)



#  library(vioplot)
#  vioplot(metaZ[overlapGenes], drawRect = F, ylab = "Schizo tissue meta keygene z-score", main = "Hyperactivity")
#  beeswarm(metaZ[overlapGenes][hpoMatrixS[,"HP:0000752"]==1], add = T, pch = 16, col = "orange2")
#  abline(h = sigThreshold, lwd = 2, col = "darkred")

#vioplot(resRecount3[overlapGenes, "Enrichment.Z.score"], drawRect = F, ylab = "IBD recount3 keygene z-score", main = "Colitis")
#beeswarm(resRecount3[overlapGenes, "Enrichment.Z.score"][hpoMatrixS[,"HP:0002583"]==1], add = T, pch = 16, col = "orange2")
#abline(h=min(abs(resRecount3$Enrichment.Z.score)[resRecount3$Bonferroni.significant]))
#dev.off()

#ovelap = intersect(rownames(recount3HpoEnrichment), rownames(tissueMetaHpoEnrichment))

# plot(-log10(tissueMetaHpoEnrichment[ovelap,1]), -log10(recount3HpoEnrichment[ovelap,1]), xlab = "Recount3 -log10(pvalue)", ylab = "Tissue meta analysis -log10(pvalue)", pch = 16, col=adjustcolor( ifelse(tissueMetaHpoEnrichment[ovelap,2] > 1, "dodgerblue2", "orange2"), alpha.f = 0.5))
#  plot(log2(tissueMetaHpoEnrichment[ovelap,2]), log2(recount3HpoEnrichment[ovelap,2]), xlab = "Recount3 log2(OR)", ylab = "Tissue meta analysis log2(OR)", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), xlim = c(-4,6), ylim = c(-4,6))

#head(sort(recount3HpoEnrichment), n = 25)

#table(hpoMatrixS[,5])

#sum(hpoMatrixS[,5] > 0)

#colnames(hpoMatrixS)[1]



#genes <- intersect(res$Recount3$Gene.ID, names(metaZ))
#zThreshold <- -log10(0.05/length(genes))

#table(res$Recount3[genes, "Enrichment.Z.score"] >= zThreshold, metaZ[genes] >= zThreshold)
#276* 100 / (54+276)

#pdf("ibdCompare.pdf")
#par(xpd = NA)
#plot(res$Recount3[genes, "Enrichment.Z.score"], metaZ[genes], cex = 0.8, pch = 16, col=adjustcolor(ifelse(res$Recount3[genes, "Enrichment.Z.score"] >= zThreshold, "darkred", "dodgerblue2"), alpha.f = 0.2),
#     ylab = "Meta analysis of tissue specific networks",
#     xlab = "Multi tissue network",
#     main = "IBD key-gene prediction scores")
#abline(h = zThreshold)
#abline(v = zThreshold)
#dev.off()
#cor.test(res$Recount3[genes, "Enrichment.Z.score"], metaZ[genes])