# Cluster script
library(data.table)
n.causal <- 100
window   <- 50000

ensembl   <- fread("../depict2_bundle/reference_datasets/human_b37/ensgR75_protein_coding.txt", data.table=F)
ensembl   <- ensembl[ensembl[,2] !="MT",]
sim.pheno <- sapply(1:110, function(x){rnorm(nrow(ensembl))})
pheno.ids <- 11:20
  
for (i in 1:10) {
  causal.genes <- ensembl[order(sim.pheno[,i], decreasing=T)[1:n.causal],]
  write.table(causal.genes, file=paste0("causal_genes_pheno", pheno.ids[i], ".txt"), quote=F)
  causal.genes[,3] <- causal.genes[,3] - window
  causal.genes[,4] <- causal.genes[,4] + window
  write.table(causal.genes[,c(2,3,4)], file=paste0("causal_genes_50k_window_pheno", pheno.ids[i], ".bed"), quote=F, col.names=F, row.names=F, sep="\t")
}
rownames(sim.pheno) <- ensembl[,1]
colnames(sim.pheno) <- paste0("randPheno100Genes", 1:ncol(sim.pheno))
write.table(sim.pheno, file="random_pathway_database_100_genes.txt", quote=F, sep="\t")
write.table(ensembl, file="random_pathway_database_100_genes.colAnnotations", quote=F, sep="\t")
