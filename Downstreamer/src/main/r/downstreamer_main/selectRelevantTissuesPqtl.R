#remoter::server(verbose = T, port = 55556, sync = T)




remoter::client("localhost", port = 55001)#55556 55501

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


#traits <- read.delim("/groups/umcg-fg/tmp01/projects/downstreamer/PascalX_bundle/traitList.txt", header = F)$V1
traits <- read.delim("/groups/umcg-fg/tmp01/projects/downstreamer/PascalX_bundle/runRealGWAS/pqtlList.txt", header = F)$V1

eigenFiles <- read.delim("/groups/umcg-fg/tmp01/projects/downstreamer/depict2_bundle/scripts/recount3Eigen.txt", header =F)
str(eigenFiles)


x <- traits[1]
resultList <- lapply(traits, function(x){
  res <- read.depict2(paste0("output/ds2/",x,"/",x,"_TissueEnrichment_enrichtments.xlsx"))[["tissue"]]
  
  sigTissues <- res$Gene.set[res$Enrichment.Z.score >0 & res$Bonferroni.significant]
  
  eigenSelection <- eigenFiles[eigenFiles$V1 %in% sigTissues,]
  str(eigenSelection)
  
  write.table(eigenSelection, file = paste0("output/ds2/",x,"/",x,"_eigenFilesSignificant.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  
}) 



