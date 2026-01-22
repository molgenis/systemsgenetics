

#remoter::server(verbose = T, port = 55556, sync = T)




remoter::client("localhost", port = 55556)#55556 55501

setwd("/groups/umcg-fg/tmp01/projects/downstreamer/depict2_bundle/")


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


geneInfo <- read.delim("/groups/umcg-fg/tmp01/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/genes_Ensembl94.txt")
str(geneInfo)

library(readr)

#change X1 in case of specified header for row names
colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("reference_datasets/human_b37/hpo/phenotype_to_genes_V1268_OMIMandORPHA.txt_matrix.txt.gz", delim = "\t", quote = "", col_types = colTypes)
hpoMatrix <- as.matrix(table_tmp[,-1])
rownames(hpoMatrix) <- table_tmp[,1][[1]]
rm(table_tmp)

hpoAnnotations <- as.matrix(read.delim("reference_datasets/human_b37/hpo/phenotype_to_genes_V1268_OMIMandORPHA.txt_matrix.colAnnotations.txt", row.names = 1))
hpoAnnotations <- hpoAnnotations[,1]


hpoPerGene <- apply(hpoMatrix, 1, function(x){sum(x>0)})
hpoMatrix <- hpoMatrix[hpoPerGene >= 1,]
str(hpoMatrix)

traits <- read.delim("/groups/umcg-fg/tmp01/projects/downstreamer/PascalX_bundle/runRealGWAS/pqtlList.txt", header = F)$V1

str(traits)

metaZscores <- matrix(NA, nrow = nrow(geneInfo), ncol = length(traits), dimnames = list(geneInfo$Gene.stable.ID, traits))
recount3Zscores <- matrix(NA, nrow = nrow(geneInfo), ncol = length(traits), dimnames = list(geneInfo$Gene.stable.ID, traits))



trait <- traits[1]
trait <- "IBD_deLange2017"
trait <- "Schizophrenia_Pardinas2018"
trait <- "ADHD_Demontis2018"

hpoEnrichmentList <- lapply(traits, function(trait){
  
  cat(trait,"\n")
  res <- read.depict2(paste0("output/ds2/", trait,"/",trait,"_keygenes_covCor_enrichtments.xlsx"))
  
  
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
  
  
 # dataForO <- as.data.frame(tissueZscores)
#  dataForO$meta <- metaZ
  
 # write.table(dataForO, file = "ibdData.txt", sep = "\t", quote = F, col.names = NA)
  
  sigThreshold <- -qnorm(0.05/length(metaZ)/2)
  sig <- metaZ >= sigThreshold
  
  
  overlapGenes <- intersect(names(sig), row.names(hpoMatrix))
  
  sig <- sig[overlapGenes]
  print("test3")
  # if(sum(sig) <= 10){ 
       tissueMetaHpoEnrichment <- NA
  #   }else {
  #     
  #     
  #   hpoMatrixS <- hpoMatrix[overlapGenes,]
  #   hpoMatrixS <- hpoMatrixS[,apply(hpoMatrixS, 2, function(x){sum(x>0)})>=10]
  #   hpoMatrixS <- hpoMatrixS[,apply(hpoMatrixS, 2, function(x){sum(x<1)})>=10]
  #   
  #   x <- hpoMatrixS[,1]
  #   tissueMetaHpoEnrichment <- t(apply(hpoMatrixS, 2, function(x){
  #     tryCatch( 
  #       {  
  #         
  #         f <- fisher.test(table(sig, x))
  #         return(c(f$p.value, f$estimate, as.vector(table(sig, x))))
  #         
  #         
  #       },
  #       error=function(cond) {
  #         print(table(sig, x)) 
  #       }
  #     )
  #   }))
  #   print("test5")
  #   
  #   #hist(-log10(tissueMetaHpoEnrichment[,1]))
  #   tissueMetaHpoEnrichment <- as.data.frame(tissueMetaHpoEnrichment)
  #   colnames(tissueMetaHpoEnrichment)[1] <- "Fisher-Pvalue"
  #   colnames(tissueMetaHpoEnrichment)[3:6] <- c("noKeyNoHpo", "keyNoHpo", "noKeyHpo", "keyHpo")
  #   tissueMetaHpoEnrichment$descip <- hpoAnnotations[rownames(tissueMetaHpoEnrichment)]
  #   
  # }
  print("test4")
  
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
  
  resRecount3 <- res[["Recount3"]]
  str(resRecount3)
  recount3Zscores[resRecount3$Gene.ID, trait] <<- resRecount3$Enrichment.Z.score
  
  
  overlapGenes <- intersect(resRecount3$Gene.ID, row.names(hpoMatrix))
  
  resRecount3 <- resRecount3[overlapGenes,]
  hpoMatrixS <- hpoMatrix[overlapGenes,]
  hpoMatrixS <- hpoMatrixS[,apply(hpoMatrixS, 2, function(x){sum(x>0)})>=10]
  
  hpoMatrixS <- hpoMatrixS[,apply(hpoMatrixS, 2, function(x){sum(x<1)})>=10]
  
  resRecount3$key <- resRecount3$Enrichment.Z.score > 0 & resRecount3$Bonferroni.significant
  table(resRecount3$key)
  print("test2")
  # if(sum(resRecount3$key) <= 10){
     recount3HpoEnrichment <- NA
  # } else {
  #   
  #   recount3HpoEnrichment <-t(apply(hpoMatrixS, 2, function(x){
  #     tryCatch( 
  #       {  
  #         
  #         f <- fisher.test(table(resRecount3$key, x))
  #         return(c(f$p.value, f$estimate, as.vector(table(resRecount3$key, x))))
  #         
  #         
  #       },
  #       error=function(cond) {
  #         print(table(resRecount3$key, x)) 
  #       }
  #     )
  #   }))
  #   recount3HpoEnrichment <- as.data.frame(recount3HpoEnrichment)
  #   colnames(recount3HpoEnrichment)[1] <- "Fisher-Pvalue"
  #   colnames(recount3HpoEnrichment)[3:6] <- c("noKeyNoHpo", "keyNoHpo", "noKeyHpo", "keyHpo")
  #   recount3HpoEnrichment$descip <- hpoAnnotations[rownames(recount3HpoEnrichment)]
  # }
  
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
  print("test")
  return(list(recount3HpoEnrichment, tissueMetaHpoEnrichment))

}) 

names(hpoEnrichmentList) <- traits

hist(metaZscores)
dev.off()

str(metaZscores)








metaZscores <- metaZscores[apply(metaZscores, 1, function(x){any(!is.na(x))}),]
recount3Zscores <- recount3Zscores[apply(recount3Zscores, 1, function(x){any(!is.na(x))}),]


str(metaZscores)


#save(recount3Zscores, metaZscores, traits, hpoEnrichmentList, hpoAnnotations, file = "hpoEnrichmentSession2.RData")

load(file = "hpoEnrichmentSession2.RData")


str(metaZscores)


metaZscoresSig <- apply(metaZscores, 2, function(x){
  t <- -qnorm((0.05/sum(!is.na(x)))/2)
  
  x[!is.na(x)] <- x[!is.na(x)] >= t
  
  return(x)
  })




#write.table(metaZscoresSig, file = gzfile("signficantMetaKeyGene.txt.gz") , sep = "\t", quote = F, col.names = NA)
#write.table(metaZscores, file = gzfile("metaKeyGene.txt.gz") , sep = "\t", quote = F, col.names = NA)

str(metaZscoresSig)


apply(metaZscoresSig, 2, sum, na.rm=T)
apply(metaZscoresSig, 2, function(x){sum(is.na(x))})

table(metaZscoresSig, useNA = "a")


metaZscoresSig2 <- metaZscoresSig[apply(metaZscoresSig,1,function(x){any(x>0,na.rm=T)}),apply(metaZscoresSig,2,function(x){any(x>0,na.rm=T)})]

distGenes <- dist(metaZscoresSig2, method = "manhattan")
distPqtl <- dist(t(metaZscoresSig2), method = "manhattan")

plot(hclust(distGenes))
dev.off()

library(heatmap3)
heatmap3(metaZscoresSig2,Rowv = as.dendrogram(hclust(distGenes)), Colv = as.dendrogram(hclust(distPqtl)), scale = "none")
dev.off()


str(metaZscoresSig)

metaZscoresSig["ENSG00000185811",]

table(metaZscoresSig2)

trait <- "Alzheimers_Jansen2019"
trait <- "IBD_deLange2017"

trait <- traits[1]

hpoE <- hpoEnrichmentList[[1]]

x <- lapply(traits, function(trait){
  
  print(trait)
  tissueMetaHpoEnrichment <- hpoEnrichmentList[[trait]][[1]]
  
  if(length(tissueMetaHpoEnrichment)>1){
  
  
  
  bonThres <- 0.05/nrow(tissueMetaHpoEnrichment)
  
  
  tissueMetaHpoEnrichment$log2Or <- log2(tissueMetaHpoEnrichment$`odds ratio`)
  
  tissueMetaHpoEnrichmentFilter <- tissueMetaHpoEnrichment[tissueMetaHpoEnrichment$`Fisher-Pvalue` <= bonThres & abs(tissueMetaHpoEnrichment$log2Or) >= 2,] 

    
    
  print(tissueMetaHpoEnrichmentFilter[order(tissueMetaHpoEnrichmentFilter$`Fisher-Pvalue`, decreasing=F)[1:10],], width = 500)
  #print(tissueMetaHpoEnrichmentFilter[order(abs(tissueMetaHpoEnrichmentFilter$`log2Or`), decreasing=T)[1:10],], width = 500)
  } else {
    print(" ")
  }
})




hpoEnrichOrMatrix <- matrix(NA, nrow = length(hpoAnnotations), ncol = length(traits), dimnames = list(hpoAnnotations, traits))


trait <- traits[1]

x <- lapply(traits, function(trait){
  
  print(trait)
  tissueMetaHpoEnrichment <- hpoEnrichmentList[[trait]][[2]]
  
  if(length(tissueMetaHpoEnrichment)>1){
    
    bonThres <- 0.05/nrow(tissueMetaHpoEnrichment)
    
    
    tissueMetaHpoEnrichment$log2Or <- log2(tissueMetaHpoEnrichment$`odds ratio`)
    
    tissueMetaHpoEnrichmentFilter <- tissueMetaHpoEnrichment[tissueMetaHpoEnrichment$`Fisher-Pvalue` <= bonThres  & abs(tissueMetaHpoEnrichment$log2Or) >= 2,] 
    
    hpoEnrichOrMatrix[tissueMetaHpoEnrichmentFilter$descip,trait] <<- tissueMetaHpoEnrichmentFilter$log2Or
    
  }
  
})

hpoEnrichOrMatrix2 <- hpoEnrichOrMatrix[apply(hpoEnrichOrMatrix, 1, function(x){any(!is.na(x))}), apply(hpoEnrichOrMatrix, 2, function(x){any(!is.na(x))})]
str(hpoEnrichOrMatrix2)






all(colnames(hpoMatrix) %in% names(hpoAnnotations))
hpoMatrix2 <- hpoMatrix
colnames(hpoMatrix2) <- hpoAnnotations[colnames(hpoMatrix)]
hpoMatrix2 <- hpoMatrix2[,rownames(hpoEnrichOrMatrix2)]

hpoEnrichOrMatrix2ForCluster <- hpoEnrichOrMatrix2

hpoEnrichOrMatrix2ForCluster[is.na(hpoEnrichOrMatrix2ForCluster)] <- 0

hist(hpoEnrichOrMatrix2ForCluster)
dev.off()


traitDist <- as.dist(1 - cor(metaZscores[,colnames(hpoEnrichOrMatrix2)], use = "pa"))
traitDist[is.nan(traitDist)] <- 1
traitDist[is.na(traitDist)] <- 1
str(traitDist)

traitClust <- hclust(traitDist)
plot(traitClust)
dev.off()

maxLog2Or <- max(abs(hpoEnrichOrMatrix2[!is.infinite(hpoEnrichOrMatrix2)]), na.rm = T)

hpoEnrichOrMatrix2[is.na(hpoEnrichOrMatrix2)] <- 0


write.table(hpoEnrichOrMatrix2, file = "hpoErichmentHeatmap2.txt", sep = "\t", quote = F, col.names = NA)


hpoEnrichOrMatrix2[is.infinite(hpoEnrichOrMatrix2) & hpoEnrichOrMatrix2 < 0] <- -maxLog2Or
hpoEnrichOrMatrix2[is.infinite(hpoEnrichOrMatrix2) & hpoEnrichOrMatrix2 > 0] <- maxLog2Or



colHeatmap <- rev(c(colorRampPalette(c("#f03b20", "#feb24c", "#ffeda0"))(99), "white", colorRampPalette(c("#e0ecf4", "#9ebcda", "#8856a7"))(99)))
colBreaks <- c(seq(-maxLog2Or,-2,length.out= 100), seq(2,maxLog2Or,length.out= 100))

str(hpoMatrix2)
hpoDist <- dist(t(hpoMatrix2), method = "manhattan")
hpoClust <- hclust(hpoDist, method =  "ward.D")
str(hpoClust)
library(pheatmap)
#rpng(width = 600, height = 1000)
pdf("hpoErichmentHeatmap2.pdf", height = 80, width = 15)
pheatmap(hpoEnrichOrMatrix2, Rowv  = as.dendrogram(hpoClust), Colv = as.dendrogram(traitClust), col = colHeatmap, breaks = colBreaks, scale = "none", cellwidth = 10, cellheight = 10)
dev.off()
sum(is.na(hpoEnrichOrMatrix2ForCluster))
range(hpoEnrichOrMatrix2ForCluster)





pascalResDir <- "/groups/umcg-fg/tmp01/projects/downstreamer/PascalX_bundle/results/"





recount3Zscores[resRecount3$Gene.ID, trait] <<- resRecount3$Enrichment.Z.score



pascalZscores <- matrix(NA, nrow = nrow(geneInfo), ncol = length(traits), dimnames = list(geneInfo$Gene.stable.ID, traits))



trait <- traits[1]
trait <- "IBD_deLange2017"
trait <- "Schizophrenia_Pardinas2018"
trait <- "ADHD_Demontis2018"

hpoEnrichmentList <- lapply(traits, function(trait){
  
  res <- read.delim(paste0(pascalResDir, trait, ".txt"))
  
  res$zscore <- qnorm(res$pvalue/2)
  

  pascalZscores[res$gene, trait] <<- res$zscore
  
  
})

pascalZscores <- pascalZscores[apply(pascalZscores, 1, function(x){any(!is.na(x))}),]


pcGenes <- read.delim("reference_datasets/human_b37/genes_Ensembl94_protein_coding.txt", header = F)$V1

pascalZscores <- pascalZscores[rownames(pascalZscores) %in% pcGenes,]

qnorm((0.05/nrow(pascalZscores))/2)

write.table(pascalZscores, file = gzfile("pascalGeneZ.txt.gz") , sep = "\t", quote = F, col.names = NA)
