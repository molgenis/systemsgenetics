#remoter::server(verbose = T, port = 55557, password = "laberkak", sync = T)

remoter::client("localhost", port = 55580, password = "laberkak")

setwd("/groups/umcg-wijmenga/tmp04/projects/depict2/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

source(paste0("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\Downstreamer\\src\\main\\r\\downstreamer_main/downstreamer_functions.r"))


library("readr")


if(FALSE){
  #don't run by accident
  table_tmp <- read_delim("depict2_bundle/reference_datasets/human_b37/pathway_databases/gene_network_v2_1/eigenvectors_165.txt", delim = "\t", quote = "")
  eigen <- as.matrix(table_tmp[,-1])
  rownames(eigen) <- table_tmp[,1][[1]]
  rm(table_tmp)
  saveRDS(eigen, "depict2_bundle/reference_datasets/human_b37/pathway_databases/gene_network_v2_1/eigenvectors_165.rds")
}
eigen <- readRDS("depict2_bundle/reference_datasets/human_b37/pathway_databases/gene_network_v2_1/eigenvectors_165.rds")


str(eigen)

genes <- read.delim("depict2_bundle/reference_datasets/human_b37/ensgR75_protein_coding.txt", stringsAsFactors = F)

eigen2 <- eigen[row.names(eigen) %in% genes$Ensembl.Gene.ID,]
str(eigen2)
library(uwot)

geneUmap <- umap(eigen2, 
                 n_threads = 24, 
                 n_epochs = 1000, 
                 init = "pca",  
                 n_neighbors = 100, 
                 min_dist = 1, 
                 init_sdev = 1e-4, 
                 learning_rate = 100, 
                 spread = 100,
                 scale = "none",
                 nn_method = "fnn")

rownames(geneUmap) <- rownames(eigen2)
colnames(geneUmap) <- colnames(c("UMAP1", "UMAP2"))

write.table(geneUmap, file = "Umap/geneUmap.txt", sep = "\t", quote = F)


geneUmap <- read.delim("Umap/geneUmap.txt", row.names = 1)
str(geneUmap)

plot(geneUmap, pch = 16, cex = 0.5, col = adjustcolor("dodgerblue2", alpha.f = 0.3))
dev.off()
  
colnames(dsBonSig)

trait = "Inflammatory bowel disease"

for(trait in colnames(dsBonSig)){


colVec <- rep(adjustcolor("dodgerblue2", alpha.f = 0.2), times = nrow(geneUmap))
colVec[dsBonSig[rownames(geneUmap), trait] & dsGeneZ[rownames(geneUmap), trait] > 0] <- adjustcolor("firebrick", alpha.f = 0.5)

plotOrder <- rep(0, times = nrow(geneUmap))
plotOrder[dsBonSig[rownames(geneUmap), trait] & dsGeneZ[rownames(geneUmap), trait] > 0] <- 1
plotOrder <- order(plotOrder,decreasing = F)

png(paste0("Umap/geneUmap/", trait,".png"), width = 1000, height = 1000)
plot(geneUmap[plotOrder,], pch = 16, cex = 0.5, col = colVec[plotOrder], main = trait)
  dev.off()
}


table_tmp <- read_delim("Umap/geneUmap/gtexV8Expression.txt", delim = "\t", quote = "")
gtex <- as.matrix(table_tmp[,-1])
rownames(gtex) <- table_tmp[,1][[1]]
rm(table_tmp)

table_tmp <- read_delim("Umap/geneUmap/single_cell_atlas_2020_11_14_geneZscores.txt", delim = "\t", quote = "")
sca <- as.matrix(table_tmp[,-1])
rownames(sca) <- table_tmp[,1][[1]]
rm(table_tmp)


colfunc<-colorRampPalette(c("lightblue1","gray90", "firebrick"))
colMapTissues <- colfunc(20)
colMapTissues2 <- adjustcolor(colMapTissues, alpha.f = 0.3)

tissue <- "Brain - Cortex"
for(tissue in colnames(gtex)){
  sharedGenes <- intersect(rownames(geneUmap), rownames(gtex))
  expression <- gtex[sharedGenes,tissue]
  m <- max(abs(expression))
  
  geneUmapShared <- geneUmap[sharedGenes,]
  
  colVec <- colMapTissues2[cut(expression, breaks = seq(-m,m,length.out = 20))]
  png(paste0("Umap/geneUmap/gtex/", tissue,".png"), width = 1000, height = 1000)
  plot(geneUmapShared, pch = 16, cex = 1, col = colVec, main = tissue)
  dev.off()
}


for(tissue in colnames(sca)){
  sharedGenes <- intersect(rownames(geneUmap), rownames(sca))
  expression <- sca[sharedGenes,tissue]
  m <- max(abs(expression))
  
  geneUmapShared <- geneUmap[sharedGenes,]
  
  colVec <- colMapTissues2[cut(expression, breaks = seq(-m,m,length.out = 20))]
  png(paste0("Umap/geneUmap/sca/", tissue,".png"), width = 1000, height = 1000)
  plot(geneUmapShared, pch = 16, cex = 1, col = colVec, main = tissue)
  dev.off()
}




reactome <- readRDS("pathwayDatabases/reactome_2020_07_18_raw.rds")
str(reactome)
colnames(reactome)
sharedGenes <- intersect(rownames(geneUmap), rownames(reactome))
reactome <- reactome[sharedGenes,]

pathway <- "R-HSA-112316"

colVec <- rep(adjustcolor("dodgerblue2", alpha.f = 0.2), times = nrow(reactome))
colVec[reactome[,pathway] == 1] <- adjustcolor("firebrick", alpha.f = 0.5)

plotOrder <- rep(0, times = nrow(geneUmap))
plotOrder[reactome[,pathway] == 1] <- 1
plotOrder <- order(plotOrder,decreasing = F)

geneUmapShared <- geneUmap[sharedGenes,]

plot(geneUmapShared[plotOrder ,], pch = 16, cex = 0.5, col = colVec[plotOrder], main = pathway)
