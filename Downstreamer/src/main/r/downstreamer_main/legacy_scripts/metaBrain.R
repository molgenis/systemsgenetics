#remoter::server(verbose = T, port = 55557, password = "laberkak", sync = T)

remoter::client("localhost", port = 55557, password = "laberkak")

setwd("/groups/umcg-wijmenga/tmp04/projects/depict2/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

source(paste0("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\Downstreamer\\src\\main\\r\\downstreamer_main/downstreamer_functions.r"))


library("readr")

if(FALSE){
  #don't run by accident
  table_tmp <- read_delim("/groups/umcg-biogen/tmp04/downstreamer/reference_datasets/pathway_databases/coregulation_all_samples/brain_coreg_zscores_cis_removed.txt.txt", delim = "\t", quote = "")
  metaBrainCoExp <- as.matrix(table_tmp[,-1])
  rownames(metaBrainCoExp) <- table_tmp[,1][[1]]
  rm(table_tmp)
  saveRDS(metaBrainCoExp, "/groups/umcg-biogen/tmp04/downstreamer/reference_datasets/pathway_databases/coregulation_all_samples/brain_coreg_zscores_cis_removed.rds")
}
metaBrainCoExp <- readRDS("MetaBrain//brain_coreg_zscores_cis_removed.rds")


str(metaBrainCoExp)



enrichFiles <- list.files(path="MetaBrain//", pattern="*_enrichtments", full.names=TRUE, recursive=FALSE)
traits <- gsub(".*/", "", enrichFiles)
traits <- gsub("\\_harm_jan_enrichtments.xlsx", "", traits)
traits <- gsub("\\_enrichtments.xlsx", "", traits)
traits
names(enrichFiles) <- traits



genePrioritizations <- lapply(enrichFiles, function(file){
  read.depict2(file)$GenePrioritization_MetaBrain
})

str(genePrioritizations)

prioGenes <- lapply(genePrioritizations, function(x){
  x$Gene.ID[x$Bonferroni.significant ]
})
prioGenes <- unique(do.call(c, prioGenes))

cisGenes <- lapply(genePrioritizations, function(x){
  x$Gene.ID[!is.na(x$GWAS.gene.P.value) & x$GWAS.gene.P.value <= 2.570694e-06 ]
})
cisGenes <- unique(do.call(c, cisGenes))
str(cisGenes)

library(heatmap3)
heatmap3(metaBrainCoExp[prioGenes,cisGenes], scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T, margins = c(12,12))

schizophrenia_ripke2014






pdf("MetaBrain/coregulationGwasVsPrio2.pdf")

for(trait in traits){

#trait = "schizophrenia_ripke2014"
  res <- genePrioritizations[[trait]]
  prioGenes <- res$Gene.ID[res$Bonferroni.significant]
  cisGenes <- res$Gene.ID[!is.na(res$GWAS.gene.P.value) & res$GWAS.gene.P.value <=  2.570694e-06 & res$Enrichment.Z.score >= 1.96]
  
  if(length(prioGenes) >= 2 & length(cisGenes) >= 2){
    coExp <- metaBrainCoExp[prioGenes, cisGenes]
    
    x <- sapply(prioGenes, function(pg){
      cisGenes[coExp[pg,] > 3]
    })
    
    prioDrivers <- stack(x)
    colnames(prioDrivers)  <- c("gwasGene", "DS_gene")
    write.table(prioDrivers, paste0("MetaBrain/",trait,"prioDrivers.txt"), row.names = F, quote = F, sep = "\t")
    
    heatmap3(coExp, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T, margins = c(10,10), xlab = "GWAS signifcant genes", ylab = "Prioritized genes", main = trait)
  }
}

dev.off()



genePrioritizationsZscores <- lapply(genePrioritizations, function(x){
  return(x$Enrichment.Z.score[ match(genePrioritizations[["schizophrenia_ripke2014"]]$Gene.ID, x$Gene.ID)])
})
dsGeneZ <- do.call(cbind, genePrioritizationsZscores)
row.names(dsGeneZ) <- genePrioritizations[["schizophrenia_ripke2014"]]$Gene.ID

pdf("Correlation of prioritization z-scores.pdf")
heatmap3(cor(dsGeneZ), scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T, margins = c(15,15), main = "Correlation of prioritization z-scores")
dev.off()