#remoter::server(verbose = T, port = 55557, password = "laberkak", sync = T)

remoter::client("localhost", port = 55557, password = "laberkak")

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
                 n_neighbors = 200, 
                 min_dist = 1, 
                 init_sdev = 1e-4, 
                 learning_rate = 100, 
                 spread = 100,
                 scale = "none",
                 nn_method = "fnn")

rownames(geneUmap) <- rownames(eigen)
colnames(geneUmap) <- colnames(c("UMAP1", "UMAP2"))

write.table(geneUmap, file = "Umap/geneUmap.txt", sep = "\t", quote = F)


geneUmap <- read.delim("Umap/geneUmap.txt")


plot(geneUmap)
dev.off()
