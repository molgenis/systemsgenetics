
projectDir <- getwd()
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

source(paste0(projectDir, "/depict2_functions.r"))

heightDs <- read.depict2("final_paper/inflammatory_bowel_disease_2017_29906448_hg19_enrichtments_exHla.xlsx")
str(heightDs$Coregulation)



library(readr)
#table_tmp <- read_delim('phenotype_to_genes_V1268_OMIMandORPHA.txt_matrix.txt.gz', delim = "\t", quote = "")
#hpoGenes <- as.matrix(table_tmp[,-1])
#rownames(hpoGenes) <- table_tmp[,1][[1]]
#rm(table_tmp)

#saveRDS(hpoGenes, file = "phenotype_to_genes_V1268_OMIMandORPHA.txt_matrix.rds")

hpoGenes <- readRDS("phenotype_to_genes_V1268_OMIMandORPHA.txt_matrix.rds")


hpoGenes <- hpoGenes[,apply(hpoGenes, 2, sum) >= 10]

#table_tmp <- read_delim('../GeneNetwork/MGI/2020_10_20/MGI_2020_10_20_PhenoGenoMP_matrix.txt.gz', delim = "\t", quote = "")
#mgiGenes <- as.matrix(table_tmp[,-1])
#rownames(mgiGenes) <- table_tmp[,1][[1]]
#rm(table_tmp)

#saveRDS(mgiGenes, file = "../GeneNetwork/MGI/2020_10_20/MGI_2020_10_20_PhenoGenoMP_matrix.rds")

mgiGenes <- readRDS("../GeneNetwork/MGI/2020_10_20/MGI_2020_10_20_PhenoGenoMP_matrix.rds")

mgiGenes <- mgiGenes[,apply(mgiGenes, 2, sum) >= 10]

heightDsCoreg <- heightDs$Coregulation

all(heightDsCoreg$Gene.set %in% rownames(hpoGenes))
all(heightDsCoreg$Gene.set %in% rownames(mgiGenes))

signGenes <- heightDsCoreg$Gene.set[heightDsCoreg$Bonferroni.significant & heightDsCoreg$Enrichment.Z.score > 0]
otherGenes <- heightDsCoreg$Gene.set[!(heightDsCoreg$Bonferroni.significant & heightDsCoreg$Enrichment.Z.score > 0)]

signGenesInMgi <- signGenes[signGenes %in% rownames(mgiGenes)]
otherGenesInMgi <- otherGenes[otherGenes %in% rownames(mgiGenes)]



str(signGenes)


ft <- rbind(
  table(factor(hpoGenes[otherGenes,"HP:0002733"], levels = c("0", "1"))),
  table(factor(hpoGenes[signGenes,"HP:0002733"], levels = c("0", "1")))
)
(res <- fisher.test(ft))
str(res)

pvalues <- apply(hpoGenes, 2, function(x){
  ft <- rbind(
    table(factor(x[otherGenes], levels = c("0", "1"))),
    table(factor(x[signGenes], levels = c("0", "1")))
    
  )
  return(fisher.test(ft)$p.value)
})


odds <- apply(hpoGenes, 2, function(x){
  ft <- rbind(
    table(factor(x[otherGenes], levels = c("0", "1"))),
    table(factor(x[signGenes], levels = c("0", "1")))
  )
  return(fisher.test(ft)$estimate)
})


pvaluesMgi <- apply(mgiGenes, 2, function(x){
  ft <- rbind(
    table(factor(x[signGenesInMgi], levels = c("0", "1"))),
    table(factor(x[otherGenesInMgi], levels = c("0", "1")))
    
  )
  return(fisher.test(ft)$p.value)
})


oddsMgi <- apply(mgiGenes, 2, function(x){
  ft <- rbind(
    table(factor(x[otherGenesInMgi], levels = c("0", "1"))),
    table(factor(x[signGenesInMgi], levels = c("0", "1")))
  )
  return(fisher.test(ft)$estimate)
})


x <- matrix(c(18030, 1078,386,22), ncol = 2)
fisher.test(x, alternative = "greater")$p.value
fisher.test(x, alternative = "less")$p.value

1078/18030
22/386


n11 <- 19010
n12 <- 387
n21 <- 98
n22 <- 21


fisher.test(matrix(c(n11,n12,n21,22), ncol = 2), alternative = "greater")$p.value
fisher.test(matrix(c(n11,n12,n21,22), ncol = 2), alternative = "less")$p.value

387/19010
21/98




heightOdds <- data.frame(hpo = colnames(hpoGenes), pvalue = pvalues, or = odds)
heightOddsMgi <- data.frame(mp = colnames(mgiGenes), pvalue = pvaluesMgi, or = oddsMgi)

hpoDesc <- read.delim("TermToDesc.txt", stringsAsFactors = F, header = F)

mgiDesc <- read.delim("../GeneNetwork/MGI/2020_10_20/VOC_MammalianPhenotype.rpt", stringsAsFactors = F, header = F)


heightOdds$hpoDesc <- hpoDesc$V2[match(heightOdds$hpo, hpoDesc$V1)]


heightOddsMgi$mpDesc <- mgiDesc$V2[match(heightOddsMgi$mp, mgiDesc$V1)]
