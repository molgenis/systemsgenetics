#remoter::server(verbose = T)

remoter :: client("localhost")


library(readr)


if(FALSE){
  #don't run by accident
  table_tmp <- read_delim("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.txt.gz", delim = "\t", quote = "")
  gnExp <- as.matrix(table_tmp[,-1])
  rownames(gnExp) <- table_tmp[,1][[1]]
  rm(table_tmp)
  saveRDS(gnExp, "/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.rds")
}

gnExp <- readRDS("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.rds")

ensg75pc <- read.delim("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/ensgR75_protein_coding.txt", stringsAsFactors = F)

gnExp <- gnExp[row.names(gnExp) %in% ensg75pc$Ensembl.Gene.ID,]

str(gnExp)

library(umap)

custom.settings = umap.defaults
custom.settings$n_neighbors = 15
custom.settings$min_dist = 0.01
custom.settings$n_epochs = 100000
gnExpUmap <- umap(t(gnExp), config = custom.settings)

str(gnExpUmap)
#pdf("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/plots/sampleUmap.pdf")

rpng( width = 1000, height = 1000)
plot(gnExpUmap$layout)
dev.off()


sampleAnnotation <- read.delim("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved_normalizedGenes.colAnnotations.txt", stringsAsFactors = F, row.names = 1 )

all(row.names(sampleAnnotation) == row.names(gnExpUmap$layout))

table(sampleAnnotation$TissueType)
sort(table(sampleAnnotation$CellType))
sort(table(sampleAnnotation$CellLine))

sampleAnnotation$col <- "grey"
sampleAnnotation$col[sampleAnnotation$TissueType == "blood" | sampleAnnotation$TissueType == "heparinised blood"] <- "red"
sampleAnnotation$col[sampleAnnotation$TissueType == "heart"] <- "orange"
sampleAnnotation$col[sampleAnnotation$TissueType == "muscle"] <- "purple4"
sampleAnnotation$col[sampleAnnotation$TissueType == "brain"] <- "blue"
sampleAnnotation$col[sampleAnnotation$TissueType == "liver"] <- "green"
sampleAnnotation$col[sampleAnnotation$TissueType == "pancreas"] <- "magenta"


rpng( width = 1000, height = 1000)
plot(gnExpUmap$layout, col = adjustcolor(sampleAnnotation$col, alpha.f = 0.1))
dev.off()
