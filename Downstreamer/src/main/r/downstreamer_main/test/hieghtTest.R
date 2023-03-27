
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Height_test")

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


dsRes <- read.depict2("height_2018_30124842_hg19_enrichtments_9.xlsx")
dsRes2 <- read.depict2("height_2018_30124842_hg19_enrichtments_10.xlsx")
str(dsRes)



library(readr)


colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("Myotube_Enrichment_normalizedPathwayScores_ExHla.txt.gz", delim = "\t", quote = "", col_types = colTypes)
myotubeEigen <- as.matrix(table_tmp[,-1])
rownames(myotubeEigen) <- table_tmp[,1][[1]]

myotubeBetas <- read.delim(gzfile("Myotube_Enrichment_betasExHla.txt.gz"))
colnames(myotubeBetas)[2] <- "Beta"
myotube <- merge(dsRes$Myotube, myotubeBetas, by.x = "Gene.set", by.y = "X.")

str(myotube)

myotubeSelected <- myotube[myotube$FDR.5..significant,c("Gene.set", "Beta")]

myotubeGenePrio <- (myotubeEigen[,myotubeSelected$Gene.set] %*% myotubeSelected$Beta)[,1]
str(myotubeGenePrio)

hist(myotubeGenePrio)




colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("cartilage_tenosynovium_Enrichment_normalizedPathwayScores_ExHla.txt.gz", delim = "\t", quote = "", col_types = colTypes)
cartilage_tenosynoviumEigen <- as.matrix(table_tmp[,-1])
rownames(cartilage_tenosynoviumEigen) <- table_tmp[,1][[1]]

cartilage_tenosynoviumBetas <- read.delim(gzfile("cartilage_tenosynovium_Enrichment_betasExHla.txt.gz"))
colnames(cartilage_tenosynoviumBetas)[2] <- "Beta"
cartilage_tenosynovium <- merge(dsRes$cartilage_tenosynovium, cartilage_tenosynoviumBetas, by.x = "Gene.set", by.y = "X.")


cartilage_tenosynoviumSelected <- cartilage_tenosynovium[cartilage_tenosynovium$FDR.5..significant,c("Gene.set", "Beta")]

cartilage_tenosynoviumGenePrio <- (cartilage_tenosynoviumEigen[,cartilage_tenosynoviumSelected$Gene.set] %*% cartilage_tenosynoviumSelected$Beta)[,1]
str(cartilage_tenosynoviumGenePrio)

hist(cartilage_tenosynoviumGenePrio)
cat(tail(names(sort(cartilage_tenosynoviumGenePrio)), n = 100),sep = "\n")






colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("Fibro_Enrichment_normalizedPathwayScores_ExHla.txt.gz", delim = "\t", quote = "", col_types = colTypes)
FibroEigen <- as.matrix(table_tmp[,-1])
rownames(FibroEigen) <- table_tmp[,1][[1]]

FibroBetas <- read.delim(gzfile("Fibro_Enrichment_betasExHla.txt.gz"))
colnames(FibroBetas)[2] <- "Beta"
Fibro <- merge(dsRes$Fibro, FibroBetas, by.x = "Gene.set", by.y = "X.")


FibroSelected <- Fibro[Fibro$FDR.5..significant,c("Gene.set", "Beta")]

FibroGenePrio <- (FibroEigen[,FibroSelected$Gene.set] %*% FibroSelected$Beta)[,1]
str(FibroGenePrio)

hist(FibroGenePrio)
cat(tail(names(sort(FibroGenePrio)), n = 100),sep = "\n")



overlap <- intersect(names(FibroGenePrio), names(cartilage_tenosynoviumGenePrio))
plot(FibroGenePrio[overlap], cartilage_tenosynoviumGenePrio[overlap],pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.5, xlab = "Fibroblast key gene scores", ylab = "Cartilage key gene scores")
cor.test(FibroGenePrio[overlap], cartilage_tenosynoviumGenePrio[overlap])







colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("Prostate_Enrichment_normalizedPathwayScores_ExHla.txt.gz", delim = "\t", quote = "", col_types = colTypes)
ProstateEigen <- as.matrix(table_tmp[,-1])
rownames(ProstateEigen) <- table_tmp[,1][[1]]

ProstateBetas <- read.delim(gzfile("Prostate_Enrichment_betasExHla.txt.gz"))
colnames(ProstateBetas)[2] <- "Beta"
Prostate <- merge(dsRes$Prostate, ProstateBetas, by.x = "Gene.set", by.y = "X.")


ProstateSelected <- Prostate[Prostate$FDR.5..significant,c("Gene.set", "Beta")]

ProstateGenePrio <- (ProstateEigen[,ProstateSelected$Gene.set] %*% ProstateSelected$Beta)[,1]
str(ProstateGenePrio)

hist(ProstateGenePrio)
cat(tail(names(sort(ProstateGenePrio)), n = 100),sep = "\n")



overlap <- intersect(names(ProstateGenePrio), names(cartilage_tenosynoviumGenePrio))
plot(ProstateGenePrio[overlap], cartilage_tenosynoviumGenePrio[overlap],pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.5, xlab = "Prostateblast key gene scores", ylab = "Cartilage key gene scores")
cor.test(ProstateGenePrio[overlap], cartilage_tenosynoviumGenePrio[overlap])




colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("Full_Enrichment_normalizedPathwayScores_ExHla.txt.gz", delim = "\t", quote = "", col_types = colTypes)
FullEigen <- as.matrix(table_tmp[,-1])
rownames(FullEigen) <- table_tmp[,1][[1]]

FullBetas <- read.delim(gzfile("Full_Enrichment_betasExHla.txt.gz"))
colnames(FullBetas)[2] <- "Beta"
Full <- merge(dsRes2$Full, FullBetas, by.x = "Gene.set", by.y = "X.")


FullSelected <- Full[Full$FDR.5..significant,c("Gene.set", "Beta")]

FullGenePrio <- (FullEigen[,FullSelected$Gene.set] %*% FullSelected$Beta)[,1]
str(FullGenePrio)

hist(FullGenePrio)
cat(tail(names(sort(FullGenePrio)), n = 100),sep = "\n")



overlap <- intersect(names(FullGenePrio), names(cartilage_tenosynoviumGenePrio))
plot(FullGenePrio[overlap], cartilage_tenosynoviumGenePrio[overlap],pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.5, xlab = "Fullblast key gene scores", ylab = "Cartilage key gene scores")
cor.test(FullGenePrio[overlap], cartilage_tenosynoviumGenePrio[overlap])

