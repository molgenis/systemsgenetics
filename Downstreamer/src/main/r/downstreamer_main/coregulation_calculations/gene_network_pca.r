library(data.table)

gn.expression          <- fread("zcat /groups/umcg-wijmenga/tmp04/projects/depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.txt.gz", data.table=F)

rownames(gn.expression) <- gn.expression[,1]
gn.expression           <- gn.expression[,-1]

protein.coding.genes    <- read.table("/groups/umcg-wijmenga/tmp04/projects/depict2/depict2_bundle/reference_datasets/human_b37/ensgR75_protein_coding.txt", stringsAsFactors=F, sep="\t", header=T)
cellline.samples        <- read.table("annotated_and_predicted_celllines_165_component_logisitic_regression.tsv", stringsAsFactors=F)[,1]


row.ol                  <- intersect(protein.coding.genes$Ensembl.Gene.ID, rownames(gn.expression))
col.ol                  <- colnames(gn.expression)[!colnames(gn.expression) %in% cellline.samples]
gn.expression.f         <- gn.expression[row.ol, col.ol]

write.table(gn.expression.f, file="gene_network_expression_data_no_celllines_deseq_normalized_protein_coding_ensgR75_samples_as_cols.txt", sep="\t", quote=F)
write.table(t(gn.expression.f), file="gene_network_expression_data_no_celllines_deseq_normalized_protein_coding_ensgR75_samples_as_cols.txt", sep="\t", quote=F)

# Use Downstreamer to calculate correlation matrix

# Gene-gene correlation matrix
#cor.m                    <- cor(t(gn.expression.f))
cor.m           <- fread("gene_gene_correlation_matrix_gene_network_expression_data_no_celllines_deseq_normalized_protein_coding_ensgR75.txt", data.table=F)
rownames(cor.m) <- cor.m[,1]
cor.m           <- cor.m[,-1]
diag(cor.m)     <- 1

# Eigen decomposition on correlation matrix
eigens                   <- eigen(cor.m)

# Number of components to save
components               <- 5000
rownames(eigens$vectors) <- rownames(cor.m)

# Write files
write.table(eigen$values, file=paste0("eigenvalues_gene_network_no_celllines_", components, "_components.tsv"), sep="\t")
write.table(eigen$vectors, file=paste0("eigenvectors_gene_network_no_celllines_", components, "_components.tsv"), sep="\t")
