
#remoter::server(verbose = T, port = 55001, sync = T)


remoter::client("localhost", port = 55001)#55501  55556

setwd("/groups/umcg-fg/tmp02/projects/downstreamer/")

traits <- read.delim("traits.txt", header = T)
str(traits)

geneInfo <- read.delim("/groups/umcg-fg/tmp02/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/genes_Ensembl94_protein_coding.txt", header = F)
str(geneInfo)

combinedResults <- matrix(NA, nrow = nrow(geneInfo), ncol = length(traits$machine_friendly_id), dimnames = list(geneInfo$V1, traits$machine_friendly_id))

t <- traits$machine_friendly_id[1]
drain <- lapply(traits$machine_friendly_id, function(t){
  cat(t)
  tRes <- read.delim(paste0("/groups/umcg-fg/tmp02/projects/downstreamer/PascalX_bundle/results_25kb/", t, ".txt"))
  
  tRes <- tRes[tRes$gene %in% geneInfo$V1,]
  
  tRes$Zscore <- -qnorm(tRes$pvalue/2)
  
  #combinedResults[tRes$gene,t] <<- tRes$Zscore
  combinedResults[tRes$gene,t] <<- qnorm((rank(tRes$Zscore,na.last="keep")-0.5)/sum(!is.na(tRes$Zscore)))
})

str(combinedResults)

combinedResults <- combinedResults[apply(combinedResults, 1, function(x){any(!is.na(x))}),]
str(combinedResults)

diseaseClasses <- table(traits$class)
#apply(combinedResults, 2, mean, na.rm = T)

boxplot(combinedResults)
dev.off()

c<-unique(traits$class)[1]

meanPerClass <- sapply(unique(traits$class), function(c){
  print(str(c))
  traitsC <- traits$machine_friendly_id[traits$class==c]

  return(apply(combinedResults[,traitsC], 1, mean, na.rm = T))

})
str(meanPerClass)

boxplot(meanPerClass)
dev.off()


pairs(meanPerClass, upper.panel = NULL)
dev.off()

medianPerClass <- sapply(unique(traits$class), function(c){
  print(str(c))
  traitsC <- traits$machine_friendly_id[traits$class==c]
  
  return(apply(combinedResults[,traitsC], 1, median, na.rm = T))
  
})
str(meanPerClass)

meanSignal <- apply(meanPerClass, 1, mean, na.rm = T)

meanSignal2 <- data.frame(meanSignal)
boxplot(meanSignal2)
dev.off()

sum(is.na(meanSignal2))

write.table(meanSignal2, file = "/groups/umcg-fg/tmp02/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/mean_gene_pvalues_88_traits_per_class_calculated_pascalX_forceNorm_25kb.txt", sep = "\t", col.names= NA, quote = FALSE)

mean50k <- read.table("/groups/umcg-fg/tmp02/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/mean_gene_pvalues_88_traits_per_class_calculated_pascalX_forceNorm.txt")
str(mean50k)
str(meanSignal2)

g <- intersect(rownames(meanSignal2),rownames(mean50k))
plot(mean50k[g,1],meanSignal2[g,1])
dev.off()

medianSignal <- apply(meanPerClass, 1, median, na.rm = T)
str(meanSignal)
str(meanPerClass[,1])
plot(meanSignal, meanPerClass[,1], xlab = "New mean of means", ylab = "mean of brain")
plot(meanSignal, meanPerClass[,3], xlab = "New mean of means", ylab = "mean of immune")
plot(meanSignal, meanPerClass[,6], xlab = "New mean of means", ylab = "mean of skelettal")
dev.off()
cor(meanSignal, medianSignal)


colnames(meanPerClass)

medianSignal <- apply(combinedResults, 1, median, na.tm = T)

str(medianSignal)
library(beeswarm)
rpng(width = 600 , height = 600)
beeswarm(cor(medianSignal, combinedResults, use = "complete.obs"), main = "Correlation between median GWAS signal and traits")
dev.off()

medianSignal2 <- data.frame(medianGwasSignal =medianSignal[!is.na(medianSignal)])

write.table(medianSignal2, file = "/groups/umcg-fg/tmp01/projects/downstreamer/depict2_bundle/medianGwasSignal.txt", sep = "\t", col.names= NA, quote = FALSE)

combinedResultsAndMedian <- cbind(combinedResults, medianSignal)


library(heatmap3)
rpng(width = 1000 , height =1000)
pdf("/groups/umcg-fg/tmp01/projects/downstreamer/depict2_bundle/traitCor.pdf", width = 15, height = 15)
heatmap3(cor(combinedResultsAndMedian, use = "complete.obs"), scale = "none", balanceColor= T, main = "Traits correlation and cor to median")
dev.off()



ldScores <- read.delim("/groups/umcg-fg/tmp01/projects/genenetwork/onco_downstreamer/ld_score/mean_ldscore_baselineLF_v2.2.UKB_Ensembl94_per_gene_25kb_window.tsv", row.names = 1)

ldScores2 <- read.delim("/groups/umcg-fg/tmp01/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/mean_ldscore_eur_per_gene_25kb_window.txt", row.names = 1)
str(ldScores)
str(ldScores2)

combinedResultsAndMedianAndLd <- merge(combinedResultsAndMedian, ldScores, by = 0, all = T)
str(combinedResultsAndMedianAndLd)

plot(combinedResultsAndMedianAndLd$ldscore, combinedResultsAndMedianAndLd$medianSignal, xlab = "Mean LD score", ylab = "Median GWAS Z-score", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
cor.test(combinedResultsAndMedianAndLd$ldscore, combinedResultsAndMedianAndLd$medianSignal, use = "complete.obs")
dev.off()
#


oldMean <- read.delim("/groups/umcg-fg/tmp01/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/mean_gene_pvalues_44_traits_per_class_calculated.txt", row.names = 1)
str(oldMean)



oldAndNewFn <- merge(oldMean, meanSignal2, by = 0 )
str(oldAndNewFn)
plot(oldAndNewFn$mean, oldAndNewFn$meanSignal, xlab = "old", ylab = "new")
dev.off()
cor(oldAndNewFn$mean, oldAndNewFn$meanSignal)

rpng(width = 1000, height = 500)
layout(matrix(1:2, nrow = 1))
oldAndLd <- merge(oldMean, ldScores2, by = 0 )
str(oldAndLd)
cor.test(oldAndLd$mean, oldAndLd$ldscore)
plot(oldAndLd$ldscore, oldAndLd$mean, xlab = "LD scores", ylab = "Old mean gwas", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
dev.off()


newAndLdFn <- merge(meanSignal2, ldScores, by = 0 )
str(newAndLdFn)
cor.test(log2(newAndLdFn$ldscore), newAndLdFn$meanSignal)
plot(log2(newAndLdFn$ldscore), newAndLdFn$meanSignal, xlab = "LD scores", ylab = "New mean gwas", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
dev.off()













oldMean <- read.delim("/groups/umcg-fg/tmp02/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/mean_gene_pvalues_44_traits_per_class_calculated.txt", row.names = 1)
colnames(oldMean) <- "oldMean"
str(oldMean)



ldScores2 <- read.delim("/groups/umcg-fg/tmp02/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/mean_ldscore_eur_per_gene_25kb_window.txt", row.names = 1)
colnames(ldScores2) <- "ldScoreMean"
str(ldScores2)


globalMedian <- read.delim("/groups/umcg-fg/tmp02/projects/downstreamer/depict2_bundle/medianGwasSignal.txt", row.names = 1)
colnames(globalMedian) <- "globalMedian"
str(globalMedian)



mean88fn <- read.delim("/groups/umcg-fg/tmp02/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/mean_gene_pvalues_88_traits_per_class_calculated_pascalX_forceNorm.txt", row.names = 1)
colnames(mean88fn) <- "mean88fn"
str(mean88fn)




mean88fn_25k <- read.delim("/groups/umcg-fg/tmp02/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/mean_gene_pvalues_88_traits_per_class_calculated_pascalX_forceNorm_25kb.txt", row.names = 1)
colnames(mean88fn_25k) <- "mean88fn"
str(mean88fn)




mean88 <- read.delim("/groups/umcg-fg/tmp02/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/mean_gene_pvalues_88_traits_per_class_calculated_pascalX.txt", row.names = 1)
colnames(mean88fn) <- "mean88"
str(mean88fn)


trait <- "Height_2022"
resMean88fn <- read.depict2(paste0("depict2_bundle/output/ds2_B/", trait,"/",trait,"_keygenes_enrichtments.xlsx"))
resMean88 <- read.depict2(paste0("depict2_bundle/output/ds2_B/", trait,"/",trait,"_keygenes_test_mean88_enrichtments.xlsx"))
resOldMean <- read.depict2(paste0("depict2_bundle/output/ds2_B/", trait,"/",trait,"_keygenes_test_oldMean_enrichtments.xlsx"))
resGlobalMedian <- read.depict2(paste0("depict2_bundle/output/ds2_B/", trait,"/",trait,"_keygenes_test_globalMedian57_enrichtments.xlsx"))
resLdscore <- read.depict2(paste0("depict2_bundle/output/ds2_B/", trait,"/",trait,"_keygenes_test_ldScore_enrichtments.xlsx"))



sharedGenes <- intersect(intersect(intersect(intersect(rownames(oldMean), rownames(ldScores2)), intersect(rownames(globalMedian), rownames(mean88fn))), rownames(mean88)),rownames(mean88fn_25k))
str(sharedGenes)

covariateMatrix <- cbind(ldScores2[sharedGenes,], oldMean[sharedGenes,], globalMedian[sharedGenes,], mean88fn[sharedGenes,], mean88[sharedGenes,], mean88fn_25k[sharedGenes,])
rownames(covariateMatrix) <- sharedGenes
colnames(covariateMatrix) <- c("ldScoreMean", "oldMean", "globalMedian", "mean88fn", "mean88", "mean88fn_25kb")
covariateMatrix[,"ldScoreMean"] <- log2(covariateMatrix[,"ldScoreMean"])


dsZscoreMatrix <- cbind(resLdscore$cartilage.tenosynovium[sharedGenes,]$Enrichment.Z.score, resOldMean$cartilage.tenosynovium[sharedGenes,]$Enrichment.Z.score, resGlobalMedian$cartilage.tenosynovium[sharedGenes,]$Enrichment.Z.score, resMean88fn$cartilage.tenosynovium[sharedGenes,]$Enrichment.Z.score , resMean88$cartilage.tenosynovium[sharedGenes,]$Enrichment.Z.score)
rownames(dsZscoreMatrix) <- sharedGenes
colnames(dsZscoreMatrix) <- c("DS_Zscore-ldScoreMean", "DS_Zscore-oldMean", "DS_Zscore-globalMedian", "DS_Zscore-mean88fn", "DS_Zscore-mean88")



panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 2)
}

rpng(width = 1000, height = 1000)
pairs(covariateMatrix,  upper.panel = panel.cor, pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
dev.off()

rpng(width = 1000, height = 1000)
pairs(dsZscoreMatrix,  upper.panel = panel.cor, pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
dev.off()
str(dsZscoreMatrix)



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




betasMean88fn <- read.delim(paste0("depict2_bundle/output/ds2_B/", trait,"/",trait,"_keygenes_intermediates/Recount3_eigenvectorBetas.txt"))
betasMean88 <- read.delim(paste0("depict2_bundle/output/ds2_B/", trait,"/",trait,"_keygenes_test_mean88_intermediates/Recount3_eigenvectorBetas.txt"))
betasOldMean <- read.delim(paste0("depict2_bundle/output/ds2_B/", trait,"/",trait,"_keygenes_test_oldMean_intermediates/Recount3_eigenvectorBetas.txt"))
betasGlobalMedian <- read.delim(paste0("depict2_bundle/output/ds2_B/", trait,"/",trait,"_keygenes_test_globalMedian57_intermediates/Recount3_eigenvectorBetas.txt"))
betasLdscore <- read.delim(paste0("depict2_bundle/output/ds2_B/", trait,"/",trait,"_keygenes_test_ldScore_intermediates/Recount3_eigenvectorBetas.txt"))

remoter::client("localhost", port = 54104)#55556 55501


str(betaMatrix)
betaMatrix <- cbind(betasLdscore[,2], betasOldMean[,2], betasGlobalMedian[,2], betasMean88fn[,2], betasMean88[,2])
rownames(betaMatrix) <- betasMean88fn[,1]
colnames(betaMatrix) <- c("ldScoreMean", "oldMean", "globalMedian", "mean88fn", "mean88")
str(betaMatrix)



rpng(width = 1000, height = 1000)
pairs(betaMatrix,  upper.panel = panel.cor, pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
dev.off()
