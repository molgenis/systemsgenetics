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

gnExpUmap <- umap(t(gnExp))

str(gnExpUmap)
#pdf("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/plots/sampleUmap.pdf")

rpng( width = 1000, height = 1000)
plot(gnExpUmap$layout)
dev.off()

