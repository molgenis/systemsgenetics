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


if(FALSE){
  #don't run by accident
  table_tmp <- read_delim("MetaBrain/MetaBrain.eigenvectors.1000_eigenvectors.txt.gz", delim = "\t", quote = "")
  metaBrainEigen <- as.matrix(table_tmp[,-1])
  rownames(metaBrainEigen) <- table_tmp[,1][[1]]
  rm(table_tmp)
  saveRDS(metaBrainEigen, "MetaBrain/MetaBrain.eigenvectors.1000_eigenvectors.rds")
}
metaBrainEigen <- readRDS("MetaBrain/MetaBrain.eigenvectors.1000_eigenvectors.rds")




str(metaBrainCoExp)
str(metaBrainEigen)


enrichFiles <- list.files(path="MetaBrain//", pattern="^[^~].*_enrichtments", full.names=TRUE, recursive=FALSE)
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
  x$Gene.ID[x$Bonferroni.significant & x$Enrichment.Z.score > 0]
})
prioGenes <- unique(do.call(c, prioGenes))

cisGenes <- lapply(genePrioritizations, function(x){
  x$Gene.ID[!is.na(x$GWAS.gene.P.value) & x$GWAS.gene.P.value <= 2.570694e-06 ]
})
cisGenes <- unique(do.call(c, cisGenes))
str(cisGenes)

library(heatmap3)
heatmap3(metaBrainCoExp[prioGenes,cisGenes], scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T, margins = c(12,12))


trait = "schizophrenia_ripke2014"
res <- genePrioritizations[[trait]]
prioGenes <- res$Gene.ID[res$Bonferroni.significant & res$Enrichment.Z.score > 0]

schizoHpo <- read.delim("MetaBrain/schizoTerms.txt", stringsAsFactors = F, row.names = 1)
schizoReactome <- read.delim("MetaBrain/schizoTermsReactome.txt", stringsAsFactors = F, row.names = 1)
pdf("MetaBrain/schizophrenia_ripke2014_prioGenes_Hpo.pdf")
heatmap3(cbind(schizoHpo[prioGenes,c("HP.0001249", "HP.0001263")], schizoReactome[prioGenes,]), scale= "none", distfun = function(x){dist(x,method = "manhattan")}, balanceColor = T, method = "ward.D2", keep.dendro = T, margins = c(10,10), xlab = "HPO", ylab = "Prioritized genes", main = trait)
dev.off()

metaBrainEigenTrait <- metaBrainEigen[prioGenes,]

library(uwot)
geneUmap <- umap(metaBrainEigenTrait, n_threads = 24, n_epochs = 10000, init = metaBrainEigenTrait[,1:2],  n_neighbors = 20, min_dist = 0.5, init_sdev = 1e-4, learning_rate = 100, spread = 1 ,scale = "none", nn_method = "fnn")
rownames(geneUmap) <- rownames(metaBrainEigenTrait)
colnames(geneUmap) <- colnames(c("UMAP1", "UMAP2"))
plot(geneUmap, col = schizoReactome[rownames(geneUmap), "R.HSA.4839726"] + 1, pch = 16)


geneUmap <- umap(metaBrainCoExp[prioGenes,prioGenes], 
                 n_threads = 24, 
                 n_epochs = 10000, 
                 init = "pca",  
                 n_neighbors = 10, 
                 min_dist = 0, 
                 init_sdev = 1e-4, 
                 learning_rate = 100, 
                 spread = 10,
                 scale = "none",
                 nn_method = "fnn")
rownames(geneUmap) <- rownames(metaBrainEigenTrait)
colnames(geneUmap) <- colnames(c("UMAP1", "UMAP2"))
layout(matrix(1:3,ncol  = 3))
plot(geneUmap, col = schizoReactome[rownames(geneUmap), "R.HSA.4839726"] + 1, pch = 16, main = "Chromatin organization")
plot(geneUmap, col = schizoHpo[rownames(geneUmap), "HP.0001263"] + 1, pch = 16, main = "Global developmental delay")
plot(geneUmap, col = schizoHpo[rownames(geneUmap), "HP.0001249"] + 1, pch = 16, main = "Intellectual disability ")


#write.table(sampleUmap, file = paste0"MetaBrain/.txt", sep = "\t", quote = F)
#sampleUmap <- read.delim("Umap/umap2.txt")

heatmap3(cbind(schizoHpo[prioGenes,], schizoReactome[prioGenes,]), scale= "none", distfun = function(x){dist(x,method = "manhattan")}, balanceColor = T, method = "ward.D2", keep.dendro = T, margins = c(10,10), xlab = "HPO", ylab = "Prioritized genes", main = trait)


sum(schizoHpo[prioGenes,"HP.0001249"]==1 | schizoHpo[prioGenes, "HP.0001263"]==1)
sum((schizoHpo[prioGenes,"HP.0001249"]==1 | schizoHpo[prioGenes, "HP.0001263"]==1) & schizoReactome[prioGenes,])

cat(prioGenes[schizoHpo[prioGenes,"HP.0001249"]==1 | schizoHpo[prioGenes, "HP.0001263"]==1],sep="\n")

cat(cisGenes[schizoHpo[cisGenes,"HP.0001249"]==1 | schizoHpo[cisGenes, "HP.0001263"]==1],sep="\n")

cisAndPrio <- unique(prioGenes, cisGenes)

cat(cisAndPrio[schizoHpo[cisAndPrio,"HP.0001249"]==1 | schizoHpo[cisAndPrio, "HP.0001263"]==1], sep ="\n")

sum(schizoReactome[cisAndPrio,])

sum(schizoReactome[cisAndPrio,] == 1 & (schizoHpo[cisAndPrio,"HP.0001249"]==1 | schizoHpo[cisAndPrio, "HP.0001263"]==1))

heatmap3(cbind(schizoHpo[cisAndPrio,c("HP.0001249", "HP.0001263")], schizoReactome[cisAndPrio,]), scale= "none", distfun = function(x){dist(x,method = "manhattan")}, balanceColor = T, method = "ward.D2", keep.dendro = T, margins = c(10,10), xlab = "HPO", ylab = "Prioritized genes", main = trait)


pdf("MetaBrain/coregulationGwasVsPrio2.pdf")

for(trait in traits){

#trait = "schizophrenia_ripke2014"
#trait = "ALS_sumstats_EUR_ASIA"
  
  res <- genePrioritizations[[trait]]
  prioGenes <- res$Gene.ID[res$Bonferroni.significant & res$Enrichment.Z.score > 0]
  cisGenes <- res$Gene.ID[!is.na(res$GWAS.gene.P.value) & res$GWAS.gene.P.value <=  2.570694e-06 & res$Enrichment.Z.score >= 1.96]
  
  if(length(prioGenes) >= 2 & length(cisGenes) >= 2){
    
    coExp <- metaBrainCoExp[prioGenes, cisGenes, drop = F]
    
    x <- sapply(prioGenes, function(pg){
      cisGenes[coExp[pg,] > 3]
    })
    
    prioDrivers <- stack(x)
    colnames(prioDrivers)  <- c("gwasGene", "DS_gene")
    prioDrivers$coregZ <- apply(prioDrivers, 1, function(pair){
      coExp[pair[2], pair[1]]
    })
    
    prioDrivers$pascalPGwasGene <- res[prioDrivers$gwasGene,"GWAS.gene.P.value"]
    prioDrivers$pascalPDsGene <- res[prioDrivers$DS_gene,"GWAS.gene.P.value"]
    
    prioDrivers$PrioZGwasGene <- res[prioDrivers$gwasGene,"Enrichment.Z.score"]
    prioDrivers$PrioZDsGene <- res[prioDrivers$DS_gene,"Enrichment.Z.score"]
    
    prioDrivers$gwasGeneSymbol <- res[prioDrivers$gwasGene,"Gene.symbol"]
    prioDrivers$DS_geneSymbol <- res[prioDrivers$DS_gene,"Gene.symbol"]
    
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












trait = "ALS_sumstats_EUR_ASIA"

res <- genePrioritizations[[trait]]
prioGenes <- res$Gene.ID[res$Enrichment.P.value <= 1E-4 & res$Enrichment.Z.score > 0]
cisGenes <- res$Gene.ID[!is.na(res$GWAS.gene.P.value) & res$GWAS.gene.P.value <=  1e-4 & res$Enrichment.Z.score >= 0]

coExp <- metaBrainCoExp[prioGenes, cisGenes, drop = F]

x <- sapply(prioGenes, function(pg){
  cisGenes[coExp[pg,] > 3]
})

prioDrivers <- stack(x)
colnames(prioDrivers)  <- c("gwasGene", "DS_gene")
prioDrivers$coregZ <- apply(prioDrivers, 1, function(pair){
  coExp[pair[2], pair[1]]
})

prioDrivers$pascalPGwasGene <- res[prioDrivers$gwasGene,"GWAS.gene.P.value"]
prioDrivers$pascalPDsGene <- res[prioDrivers$DS_gene,"GWAS.gene.P.value"]

prioDrivers$PrioZGwasGene <- res[prioDrivers$gwasGene,"Enrichment.Z.score"]
prioDrivers$PrioZDsGene <- res[prioDrivers$DS_gene,"Enrichment.Z.score"]

prioDrivers$gwasGeneSymbol <- res[prioDrivers$gwasGene,"Gene.symbol"]
prioDrivers$DS_geneSymbol <- res[prioDrivers$DS_gene,"Gene.symbol"]


write.table(prioDrivers, paste0("MetaBrain/",trait,"prioDrivers_special.txt"), row.names = F, quote = F, sep = "\t")
pdf(paste0("MetaBrain/",trait,"prioDrivers_special.pdf"))
heatmap3(coExp, scale= "none", balanceColor = T, method = "ward.D2", keep.dendro = T, margins = c(10,10), xlab = "GWAS signifcant genes", ylab = "Prioritized genes", main = trait)
dev.off()
