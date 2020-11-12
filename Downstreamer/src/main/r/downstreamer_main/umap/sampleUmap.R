#remoter::server(verbose = T, port = 55557, password = "laberkak", sync = T)

remoter::client("localhost", port = 55557, password = "laberkak")

setwd("/groups/umcg-wijmenga/tmp04/projects/depict2/")

library(readr)


#if(FALSE){
#  #don't run by accident
#  table_tmp <- read_delim("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.txt.gz", delim = "\t", quote = "")
#  gnExp <- as.matrix(table_tmp[,-1])
#  rownames(gnExp) <- table_tmp[,1][[1]]
#  rm(table_tmp)
#  saveRDS(gnExp, "/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.rds")
#}
#gnExp <- readRDS("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.rds")

ensg75pc <- read.delim("ensgR75_protein_coding.txt", stringsAsFactors = F)

table_tmp <- read_delim("/groups/umcg-wijmenga/prm02/data_projects/Gado/GeneNetwork_V2_01-02-2018/Covariates/PCA/pc-scores250.txt.gz", delim = "\t", quote = "")
pcs <- as.matrix(table_tmp[,-1])
rownames(pcs) <- table_tmp[,1][[1]]
rm(table_tmp)
pcs <- pcs[,1:100]
str(pcs)

library(uwot)

init <- pcs[,c(2,1)]
init <- init * -1

tissueColMap <- read.delim("Umap/tissueCol.txt", stringsAsFactors = F)
sampleAnnotation <- read.delim("Umap/sampleAnnotations.txt", stringsAsFactors = F, row.names = 1 )
str(sampleAnnotation)

all(row.names(sampleAnnotation)  %in% row.names(pcs))

sampleAnnotation <- sampleAnnotation[match(row.names(pcs), rownames(sampleAnnotation) ),]

all(row.names(sampleAnnotation) == row.names(pcs))


sampleAnnotation$pch = 1
sampleAnnotation$pch[sampleAnnotation$CellLine != ""] = 5
sampleAnnotation$pch[sampleAnnotation$TissueType != ""] = 21

sampleAnnotation$col = adjustcolor("grey", alpha.f = 0.1)
sampleAnnotation$bg = NA

for(i in 1:nrow(tissueColMap)){
  sampleAnnotation$col[sampleAnnotation$TissueType == tissueColMap$Tissue[i]] <- adjustcolor(tissueColMap$col[i], alpha.f = 0.1)
  sampleAnnotation$bg[sampleAnnotation$TissueType == tissueColMap$Tissue[i]] <- adjustcolor(tissueColMap$col[i], alpha.f = 0.3)
}


 
#sampleUmap <- umap(pcs, n_threads = 24, n_epochs = 100000, init = init,  n_neighbors = 100, min_dist = 5, init_sdev = 1e-4, learning_rate = 1, spread = 10)

rpng( width = 1200, height = 1000)
layout(matrix(2:1, nrow =1),widths = c(5,1))
par(mar = c(0,0,0,0))
plot.new()
plot.window(0:1,0:1)
legend(x = "center", legend = tissueColMap$Tissues, fill = tissueColMap$col)
par(mar = c(3,3,1,0))
plot(sampleUmap, col = sampleAnnotation$col, pch = sampleAnnotation$pch, bg = sampleAnnotation$bg, cex = 0.5)
dev.off()
  