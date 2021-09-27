setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

library(readr)
eqtl <- read_delim("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz", delim = "\t", quote = "")




  #eqtlFull <- read_delim("2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz", delim = "\t", quote = "")
#saveRDS(eqtlFull, "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.rds")
#testedGenes <- unique(eqtlFull$Gene)
#write.table(testedGenes, file = "eqtlGenTestedGenes.txt", quote = F, row.names = F)
testedGenes <- read.delim("eqtlGenTestedGenes.txt")[,1]


maxEqtl <- aggregate(eqtl$Zscore, list(gene = eqtl$Gene), max)


rm(eqtlFull)

str(eqtlFull2)
hist(maxEqtlFull$x)

dim(maxEqtlFull)



eqtlFull2 <- eqtlFull[,c("Gene", "Zscore")]
maxEqtlFull <- aggregate(eqtlFull2$Zscore, list(gene = eqtlFull2$Gene), max)

load("keyCisTransRare/key_gene_overlap_files_for_patrick.RData")

qtlTestedDsTested <- intersect(testedGenes, row.names(gene.pvalues) )

str(key.genes)

allKey <- unique(do.call("c", key.genes))

str(maxEqtl)

maxEqtlTested <- maxEqtlFull[maxEqtlFull$gene %in% qtlTestedDsTested, ]

boxplot(maxEqtlTested$x ~ as.factor((maxEqtlTested$gene %in% key.genes$inflammatory_bowel_disease_2017_29906448)))
t.test(maxEqtlTested$x ~ as.factor((maxEqtlTested$gene %in% key.genes$inflammatory_bowel_disease_2017_29906448)))

boxplot(maxEqtlTested$x ~ as.factor((maxEqtlTested$gene %in% key.genes$height_2018_30124842)))
t.test(maxEqtlTested$x ~ as.factor((maxEqtlTested$gene %in% key.genes$height_2018_30124842)))

sigEqtlGenes <- unique(eqtl$Gene)

sigEqtlGenes2 <- sigEqtlGenes[sigEqtlGenes %in% qtlTestedDsTested]
length(sigEqtlGenes2)
#a = eqtl and key gene
#b = eqtl not key gene
#c = not eqtl and key gene
#d = not eqtl not key gene

keyGenes <- key.genes$inflammatory_bowel_disease_2017_29906448
keyGenes2 <- keyGenes[keyGenes %in% qtlTestedDsTested]

a <- sum(sigEqtlGenes2 %in% keyGenes2)
b <- length(sigEqtlGenes2) - a
c <- length(keyGenes2) - a
d <- length(qtlTestedDsTested) - a - b - c 

(m <- matrix(c(a,b,c,d), nrow = 2, byrow = T, dimnames = list(c("isCisQtl", "notCisQtl"), c("isKeyGene", "notKeyGene"))))
fisher.test(m)

keyGenes <- key.genes$height_2018_30124842
keyGenes2 <- keyGenes[keyGenes %in% qtlTestedDsTested]

a <- sum(sigEqtlGenes2 %in% keyGenes2)
b <- length(sigEqtlGenes2) - a
c <- length(keyGenes2) - a
d <- length(qtlTestedDsTested) - a - b - c 

(m <- matrix(c(a,b,c,d), nrow = 2, byrow = T, dimnames = list(c("isCisQtl", "notCisQtl"), c("isKeyGene", "notKeyGene"))))
fisher.test(m)

keyGenes <- key.genes$educational_attainment_2018_30038396
keyGenes2 <- keyGenes[keyGenes %in% qtlTestedDsTested]

a <- sum(sigEqtlGenes2 %in% keyGenes2)
b <- length(sigEqtlGenes2) - a
c <- length(keyGenes2) - a
d <- length(qtlTestedDsTested) - a - b - c 

(m <- matrix(c(a,b,c,d), nrow = 2, byrow = T, dimnames = list(c("isCisQtl", "notCisQtl"), c("isKeyGene", "notKeyGene"))))
str(fisher.test(m))





fisherPerTrait <- sapply(names(key.genes), function(trait){


  keyGenes <- key.genes[[trait]]
  keyGenes2 <- keyGenes[keyGenes %in% qtlTestedDsTested]
  
  if(length(keyGenes2) < 10){
    return(c(NA,NA,NA,NA))
  }
  
  a <- sum(sigEqtlGenes2 %in% keyGenes2)
  b <- length(sigEqtlGenes2) - a
  c <- length(keyGenes2) - a
  d <- length(qtlTestedDsTested) - a - b - c 
  
  
  (m <- matrix(c(a,b,c,d), nrow = 2, byrow = T, dimnames = list(c("isCisQtl", "notCisQtl"), c("isKeyGene", "notKeyGene"))))
  res <- fisher.test(m)
  return(c(res$p.value, res$estimate, res$conf.int))
})
fisherPerTrait
str(fisherPerTrait)


fisherPerTrait2 <- fisherPerTrait[,!is.na(fisherPerTrait[2,])]
str(fisherPerTrait2)

par(mar = c(15,4,2,1))
barplot(sort(log2(fisherPerTrait2[2,])), las = 2)


fisherPerTrait3 <- as.data.frame(t(fisherPerTrait[,!is.na(fisherPerTrait[2,])]))
fisherPerTrait3$log2or <- log2(fisherPerTrait3$`odds ratio`)
fisherPerTrait3$log2or_low <- log2(fisherPerTrait3$V3)
fisherPerTrait3$log2or_high <- log2(fisherPerTrait3$V4)
colnames(fisherPerTrait3) <- c("pvalue", "or", "or_low", "or_high", "log2or", "log2or_low", "log2or_high")
fisherPerTrait3$bonf <- fisherPerTrait3$pvalue < (0.05 / nrow(fisherPerTrait3))
write.table(fisherPerTrait3, file = "cisQtlEnrichment.txt", sep = "\t", quote = FALSE, col.names = NA)











eqtlBrain <- read_delim("metaBrainQtl/eQTLProbesFDR0.05-ProbeLevel.txt.gz", delim = "\t", quote = "")

brainTested <- read_delim("metaBrainQtl/testedgenes.txt", delim = "\t", quote = "", col_names = F)$X1
length(brainTested)
View(eqtlBrain)


brainTested <- gsub("\\..*$", "", brainTested)
eqtlBrain$ProbeName <- gsub("\\..*$", "", eqtlBrain$ProbeName)

qtlBrainTestedDsTested <- intersect(brainTested, row.names(gene.pvalues) )

sigEqtlGenesBrain2 <- unique(eqtlBrain$ProbeName[eqtlBrain$ProbeName %in% qtlBrainTestedDsTested])
str(sigEqtlGenesBrain2)

keyGenes <- key.genes$inflammatory_bowel_disease_2017_29906448
keyGenes2 <- keyGenes[keyGenes %in% qtlBrainTestedDsTested]

a <- sum(sigEqtlGenesBrain2 %in% keyGenes2)
b <- length(sigEqtlGenesBrain2) - a
c <- length(keyGenes2) - a
d <- length(qtlBrainTestedDsTested) - a - b - c 

(m <- matrix(c(a,b,c,d), nrow = 2, byrow = T, dimnames = list(c("isCisQtl", "notCisQtl"), c("isKeyGene", "notKeyGene"))))
fisher.test(m)

sum(keyGenes2 %in% sigEqtlGenesBrain2 )

keyGenes <- key.genes$educational_attainment_2018_30038396
keyGenes2 <- keyGenes[keyGenes %in% qtlBrainTestedDsTested]

a <- sum(sigEqtlGenesBrain2 %in% keyGenes2)
b <- length(sigEqtlGenesBrain2) - a
c <- length(keyGenes2) - a
d <- length(qtlBrainTestedDsTested) - a - b - c 

(m <- matrix(c(a,b,c,d), nrow = 2, byrow = T, dimnames = list(c("isCisQtl", "notCisQtl"), c("isKeyGene", "notKeyGene"))))
fisher.test(m)


keyGenes <- key.genes$height_2018_30124842
keyGenes2 <- keyGenes[keyGenes %in% qtlBrainTestedDsTested]

a <- sum(sigEqtlGenesBrain2 %in% keyGenes2)
b <- length(sigEqtlGenesBrain2) - a
c <- length(keyGenes2) - a
d <- length(qtlBrainTestedDsTested) - a - b - c 

(m <- matrix(c(a,b,c,d), nrow = 2, byrow = T, dimnames = list(c("isCisQtl", "notCisQtl"), c("isKeyGene", "notKeyGene"))))
fisher.test(m)



fisherPerTrait <- sapply(names(key.genes), function(trait){
  
  
  keyGenes <- key.genes[[trait]]
  keyGenes2 <- keyGenes[keyGenes %in% qtlBrainTestedDsTested]
  
  if(length(keyGenes2) < 10){
    return(c(NA,NA,NA,NA))
  }
  
  a <- sum(sigEqtlGenesBrain2 %in% keyGenes2)
  b <- length(sigEqtlGenesBrain2) - a
  c <- length(keyGenes2) - a
  d <- length(qtlBrainTestedDsTested) - a - b - c 
  
  
  (m <- matrix(c(a,b,c,d), nrow = 2, byrow = T, dimnames = list(c("isCisQtl", "notCisQtl"), c("isKeyGene", "notKeyGene"))))
  res <- fisher.test(m)
  return(c(res$p.value, res$estimate, res$conf.int))
})
fisherPerTrait
str(fisherPerTrait)


fisherPerTrait2 <- fisherPerTrait[,!is.na(fisherPerTrait[2,])]


par(mar = c(15,4,2,1))
barplot(sort(log2(fisherPerTrait2[2,])), las = 2)


fisherPerTrait3Brain <- as.data.frame(t(fisherPerTrait[,!is.na(fisherPerTrait[2,])]))
fisherPerTrait3Brain$log2or <- log2(fisherPerTrait3Brain$`odds ratio`)
fisherPerTrait3Brain$log2or_low <- log2(fisherPerTrait3Brain$V3)
fisherPerTrait3Brain$log2or_high <- log2(fisherPerTrait3Brain$V4)
colnames(fisherPerTrait3Brain) <- c("pvalue", "or", "or_low", "or_high", "log2or", "log2or_low", "log2or_high")
fisherPerTrait3Brain$bonf <- fisherPerTrait3Brain$pvalue < (0.05 / nrow(fisherPerTrait3Brain))
write.table(fisherPerTrait3Brain, file = "cisQtlEnrichmentBrain.txt", sep = "\t", quote = FALSE, col.names = NA)



dim(fisherPerTrait3Brain)
