
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\LudeEigen")

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



trait <- "inflammatory_bowel_disease_2017_29906448_hg19"
trait_i <- "8"


dsRes <- read.depict2(paste0(trait,"_enrichtments_",trait_i,".xlsx"))
str(dsRes)



library(readr)


colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim(paste0(trait,"_intermediates/EigenLude_Enrichment_normalizedPathwayScores_ExHla.txt.gz"), delim = "\t", quote = "", col_types = colTypes)
eigenLude <- as.matrix(table_tmp[,-1])
rownames(eigenLude) <- table_tmp[,1][[1]]

betasLude <- read.delim(gzfile(paste0(trait,"_intermediates/EigenLude_Enrichment_betasExHla.txt.gz")))
colnames(betasLude)[2] <- "Betabb"

lude <- merge(dsRes$EigenLude, betasLude, by.x = "Gene.set", by.y = "X.")

ludeSelected <- lude[lude$FDR.5..significant,c("Gene.set", "Beta")]

ludeGenePrio <- (eigenLude[,ludeSelected$Gene.set] %*% ludeSelected$Beta)[,1]

write.table(ludeGenePrio, file = gzfile(paste0(trait,"_ludeEigen.txt.gz")), sep = "\t", quote = F, col.names = F)






colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim(paste0(trait,"_intermediates/Full_Enrichment_normalizedPathwayScores_ExHla.txt.gz"), delim = "\t", quote = "", col_types = colTypes)
eigenrecount3 <- as.matrix(table_tmp[,-1])
rownames(eigenrecount3) <- table_tmp[,1][[1]]

betasrecount3 <- read.delim(gzfile(paste0(trait,"_intermediates/Full_Enrichment_betasExHla.txt.gz")))
colnames(betasrecount3)[2] <- "Beta"

recount3 <- merge(dsRes$Full, betasrecount3, by.x = "Gene.set", by.y = "X.")

recount3Selected <- recount3[recount3$FDR.5..significant,c("Gene.set", "Beta")]

recount3GenePrio <- (eigenrecount3[,recount3Selected$Gene.set] %*% recount3Selected$Beta)[,1]

write.table(recount3GenePrio, file = gzfile(paste0(trait,"_recount3.txt.gz")), sep = "\t", quote = F, col.names = F)

tail(sort(recount3GenePrio))
