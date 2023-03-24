setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\output\\ds2\\")

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


traits <- read.delim("../../traitList.txt", header = F)$V1

minZfdrSig <- Inf

x <- traits[1]
resultList <- lapply(traits, function(x){
  res <- read.depict2(paste0(x,"/",x,"_tissueEnrichment_enrichtments.xlsx"))[["tissueMedianCentered"]]
  tissues <<- res$Gene.set[order(res$Gene.set)]
  
  if(any(res$FDR.5..significant)){
    z <- min(abs(res$Enrichment.Z.score[res$FDR.5..significant]))
    if(z < minZfdrSig){
      minZfdrSig <<- z
    }
  }
  
  return(res[order(res$Gene.set),"Enrichment.Z.score"])
}) 

tissueEnrichments <- do.call("cbind", resultList)
rownames(tissueEnrichments) <- tissues
colnames(tissueEnrichments) <- traits


library(pheatmap)

sigZ <- -qnorm(0.05/nrow(tissueEnrichments)/2)


load("../../combinedHealthyTissue_medianPerTissueCenteredCluster.RData", verbose = T)


#"#f03b20"
colHeatmap <- c(colorRampPalette(c("#feb24c", "#ffeda0"))(99), "white", colorRampPalette(c("#e0ecf4", "#9ebcda", "#8856a7"))(99))
colBreaks <- c(seq(min(tissueEnrichments),-sigZ,length.out= 100), seq(sigZ,max(tissueEnrichments),length.out= 100))


pdf("plot.pdf", width = 20, height =20)
pheatmap(tissueEnrichments[medianPerTissueCenteredCluster$labels,], scale = "none", col = colHeatmap, breaks = colBreaks,cluster_rows = medianPerTissueCenteredCluster)
dev.off()
  
str(tissueEnrichments)

str(medianPerTissueCenteredCluster$labels)
str(tissueEnrichments)
