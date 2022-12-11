
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\IBD_test")

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


dsRes <- read.depict2("inflammatory_bowel_disease_2017_29906448_hg19_enrichtments_5.xlsx")
str(dsRes)



library(readr)


colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("Colon_eigenVec_ForceNorm.txt.gz", delim = "\t", quote = "", col_types = colTypes)
colonEigen <- as.matrix(table_tmp[,-1])
rownames(colonEigen) <- table_tmp[,1][[1]]


colonBetas <- read.delim(gzfile("Colon_Enrichment_betasExHla.txt.gz"))
colnames(colonBetas)[2] <- "Beta"
colon <- merge(dsRes$Colon, colonBetas, by.x = "Gene.set", by.y = "X.")


colonSelected <- colon[colon$FDR.5..significant,c("Gene.set", "Beta")]
#colonSelected <- colon[colon$Enrichment.P.value <= 0.05,c("Gene.set", "Beta")]

colonGenePrio <- (colonEigen[,colonSelected$Gene.set] %*% colonSelected$Beta)[,1]
str(colonGenePrio)

oldGenePrio <- dsRes$GenePrioritization
rownames(oldGenePrio) <- oldGenePrio$Gene.ID
str(oldGenePrio)

geneOverlap <- intersect(names(colonGenePrio), oldGenePrio$Gene.ID)
plot(colonGenePrio[geneOverlap], oldGenePrio[geneOverlap, "Enrichment.Z.score"], xlab = "Colon prioritization using significant components", ylab = "Coregulation gene network 165 comps", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.5)
cor.test(colonGenePrio[geneOverlap], oldGenePrio[geneOverlap, "Enrichment.Z.score"])
str(colonGenePrio[geneOverlap])


cat(names(tail(sort(colonGenePrio[geneOverlap]), n = 100)), sep = "\n")

cat(names(tail(sort(as.matrix(oldGenePrio[geneOverlap, "Enrichment.Z.score", drop =F])[,1]), n = 100)), sep = "\n")





colonSelected2 <- colon[colon$Enrichment.P.value <= 0.05,c("Gene.set", "Beta")]
colonGenePrio2 <- (colonEigen[,colonSelected$Gene.set] %*% colonSelected$Beta)[,1]

plot(colonGenePrio[geneOverlap], colonGenePrio2[geneOverlap], xlab = "Colon prioritization using FDR significant components", ylab = "Colon prioritization using nominal significant components", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.5)
cor.test(colonGenePrio, colonGenePrio2)











colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("T.Cells_eigenVec_ForceNorm.txt.gz", delim = "\t", quote = "", col_types = colTypes)
tcellsEigen <- as.matrix(table_tmp[,-1])
rownames(tcellsEigen) <- table_tmp[,1][[1]]

tcellsBetas <- read.delim(gzfile("T_Cells_Enrichment_betasExHla.txt.gz"))
colnames(tcellsBetas)[2] <- "Beta"
tcells <- merge(dsRes$T_Cells, tcellsBetas, by.x = "Gene.set", by.y = "X.")
str(tcells)

tcellsSelected <- tcells[tcells$FDR.5..significant,c("Gene.set", "Beta")]

tcellsGenePrio <- (tcellsEigen[,tcellsSelected$Gene.set, drop = F] %*% tcellsSelected$Beta)[,1]
str(tcellsGenePrio)
tail(sort(tcellsGenePrio))


oldGenePrio <- dsRes$GenePrioritization
str(oldGenePrio)

geneOverlap <- intersect(names(tcellsGenePrio), oldGenePrio$Gene.ID)
plot(tcellsGenePrio[geneOverlap], oldGenePrio[geneOverlap, "Enrichment.Z.score"], xlab = "T-cell prioritization using significant components", ylab = "Coregulation gene network 165 comps", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.5)
cor.test(tcellsGenePrio[geneOverlap], oldGenePrio[geneOverlap, "Enrichment.Z.score"])








colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("Whole.Blood_eigenVec_ForceNorm.txt.gz", delim = "\t", quote = "", col_types = colTypes)
wholeBloodEigen <- as.matrix(table_tmp[,-1])
rownames(wholeBloodEigen) <- table_tmp[,1][[1]]

wholeBloodBetas <- read.delim(gzfile("Whole_Blood_Enrichment_betasExHla.txt.gz"))
colnames(wholeBloodBetas)[2] <- "Beta"
wholeBlood <- merge(dsRes$Whole_Blood, wholeBloodBetas, by.x = "Gene.set", by.y = "X.")
str(wholeBlood)

wholeBloodSelected <- wholeBlood[wholeBlood$FDR.5..significant,c("Gene.set", "Beta")]

wholeBloodGenePrio <- (wholeBloodEigen[,wholeBloodSelected$Gene.set, drop = F] %*% wholeBloodSelected$Beta)[,1]
str(wholeBloodGenePrio)
tail(sort(wholeBloodGenePrio))


oldGenePrio <- dsRes$GenePrioritization
str(oldGenePrio)

geneOverlap <- intersect(names(wholeBloodGenePrio), oldGenePrio$Gene.ID)
plot(wholeBloodGenePrio[geneOverlap], oldGenePrio[geneOverlap, "Enrichment.Z.score"], xlab = "Whole blood prioritization using significant components", ylab = "Coregulation gene network 165 comps", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.5)
cor.test(wholeBloodGenePrio[geneOverlap], oldGenePrio[geneOverlap, "Enrichment.Z.score"])





