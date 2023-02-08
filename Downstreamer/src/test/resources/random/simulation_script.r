library(MASS)
library(cowplot)
library(data.table)

source("/nfs/users/nfs_o/ob7/tools/rscripts/plotting_functions.r")
#source("downstreamer_functions.r")

output.prefix <- "../analysis/simulated_data/sim_data_300_genes"
dir.create(output.prefix, recursive = T)

# ------------------------------------------------------------------------------
gene.info           <- fread("../../data/ensembl_108_genes_b38.txt", data.table=F)
rownames(gene.info) <- make.names(gene.info$`Gene stable ID`, unique=T)
gene.info$arm       <- sapply(gene.info$`Karyotype band`, substr, 1,1)
gene.info$arm       <- paste0(gene.info$`Chromosome/scaffold name`, gene.info$arm)
gene.info           <- gene.info[order(gene.info$`Chromosome/scaffold name`,
                                       gene.info$`Gene start (bp)`),]
# ------------------------------------------------------------------------------
# Take gene gene correlation matrix and remove correlations between chromosome
# arms
# This is just a gene gene correlation matrix, can be any correlation matrix
# in principle, doesnt really matter as it is our ground truth, but decided
# to simulate data using a real correlation matrix from the old IBD permutations
# to be as representatitve as possible.
load("../data/test_data_ibd/cache_sigma.RData")

ol                  <- intersect(rownames(gene.info), colnames(sigma))
gene.info           <- gene.info[ol,]
arms                <- unique(gene.info$arm)

sigma.hat           <- matrix(0, ncol=length(ol), nrow=length(ol))
rownames(sigma.hat) <- ol
colnames(sigma.hat) <- ol

# Set correlation between arms to zero
for (arm in arms) {
  cur.genes <- gene.info[gene.info$arm == arm, ]
  cur.genes <- cur.genes[order(cur.genes$`Gene start (bp)`), 1]
  
  cur.genes <- intersect(cur.genes, ol)
  
  # Maintain only the correlations on a single chromosome arm
  sigma.hat[cur.genes, cur.genes] <- sigma[cur.genes, cur.genes]
}

# ------------------------------------------------------------------------------
# Subsample the matrix with only genes with lots of correlation to get a 
# relevant simulation.

# Copy correlation matrix
tmp             <- sigma.hat

# Ignore low correlation genes for this goal
tmp[tmp < 0.05] <- 0

# Select genes that have at least 40 correlations > 0.05
keep            <- rowSums(tmp > 0)
keep            <- keep[keep > 40]

# Sample 300 genes from this set
keep <- sample(names(keep), 300, replace=F)

# Order by genomic position
keep <- keep[order(gene.info[keep, "Chromosome/scaffold name"],
                   gene.info[keep, "Gene start (bp)"])]

final.gene.order <- keep
cor.ld           <- sigma.hat[final.gene.order, final.gene.order]
gene.info        <- gene.info[final.gene.order,]

# Convert to the nearest PD matrix so we can compare using inverse of cholesky
cor.ld           <- nearPD(cor.ld)$mat

# Sanity check correlations
simple.hm(cor.ld, cellsize = 0.5, show_rownames=F, show_colnames=F, cluster=F)

# Cleanup
rm(sigma.hat, sigma, tmp, arm, arms, cur.genes, keep, ol)
gc()

# This final matrix is the ground truth for the simulations
# Again, how this is arrived at doesn't matter that much.
write.table(as.matrix(cor.ld), file=paste0(output.prefix, "/300G_true_correlation_matrix.txt"), sep="\t", col.names = NA, quote=F)
write.table(final.gene.order, file=paste0(output.prefix, "/300G_gene_order.txt"), sep="\t", quote=F, row.names=F, col.names=F)


# ------------------------------------------------------------------------------
# Gene annotation file
ds.gene.info <- gene.info[final.gene.order, c("Gene stable ID", 
                                  "Chromosome/scaffold name",
                                  "Gene start (bp)" ,
                                  "Gene end (bp)",
                                  "Gene type",
                                  "Karyotype band",
                                  "Gene name"
)]

write.table(ds.gene.info, file=paste0(output.prefix, "/300G_gene_info.txt"), sep="\t", quote=F, row.names=F, col.names=T)

# ------------------------------------------------------------------------------
# Simulate the datasets using the correlation matrix

# Gene pvalues
# Formally the Sigma argument is a covariance matrix, but as we are producing data
# with mean zero and sd 1, this is also valid.
data.y    <- mvrnorm(n=1, mu=rep(0, ncol(cor.ld)), Sigma=cor.ld)
data.perm <- mvrnorm(n=20100, mu=rep(0, ncol(cor.ld)), Sigma=cor.ld)

# Pathway
data.x    <- t(mvrnorm(n=1000, mu=rep(0, ncol(cor.ld)), Sigma=cor.ld))

write.table(data.y, file=paste0(output.prefix, "/300G_y_gwas.txt"), sep="\t", col.names = NA, quote=F)
write.table(t(data.perm), file=paste0(output.prefix, "/300G_20100_random_gwass.txt"), sep="\t", col.names = NA, quote=F)
write.table(data.x, file=paste0(output.prefix, "/300G_x_pathways.txt"), sep="\t", col.names = NA, quote=F)

# ------------------------------------------------------------------------------
# Eigen decompositions
# ------------------------------------------------------------------------------

# Block diagonal decomposition
# ------------------------------------------------------------------------------
blocks            <- data.frame(matrix(NA, nrow=nrow(gene.info), ncol=0))
blocks$idx        <- 1:nrow(blocks)
blocks$new_idx    <- 1:nrow(blocks)
blocks$gene       <- final.gene.order
blocks$block      <- gene.info[final.gene.order, "arm"]

write.table(blocks, file=paste0(output.prefix, "/300G_gene_blocks.txt"), sep="\t", quote=F, row.names=F)

values            <- c()
vectors           <- matrix(0, nrow=nrow(cor.ld), ncol(cor.ld))

master.index <- 1
for (block in unique(blocks$block)) {
  genes           <- blocks[blocks$block == block, "gene"]
  idx             <- blocks[blocks$block == block,"idx"] 
  
  e               <- eigen(as.matrix(cor.ld[idx, idx]))
  v               <- e$vectors
  
  for (j in 1:ncol(v)) {
    for (i in 1:nrow(v)) {
      vectors[idx[i], master.index] <- v[i,j]
    }
    
    master.index <- master.index + 1
  }
  
  values                    <- c(values, e$values)
}

rownames(vectors) <- final.gene.order

# Sort by eigenvalue
blocked.vectors <- vectors[,order(values, decreasing=T)]
blocked.values  <- sort(values, decreasing=T)

write.table(blocked.values, file=paste0(output.prefix, "/300G_r_calculated_blocked_eigenvalues.txt"), sep="\t", col.names = NA, quote=F)
write.table(blocked.vectors, file=paste0(output.prefix, "/300G_r_calculated_blocked_eigenvectors.txt"), sep="\t", col.names = NA, quote=F)

rm(vectors, values)
# ------------------------------------------------------------------------------

full.e       <- eigen(cor.ld)
full.vectors <- full.e$vectors
full.values  <- full.e$values

rownames(full.vectors) <- rownames(cor.ld)

# Eigenvalues are well correlated, nearly identical
plot(blocked.values, full.values)

# Eigenvectors are well correlated (nearly identical except for very small eigenvalues)
plot(blocked.vectors[,1], full.vectors[,1])
plot(blocked.vectors[,200], full.vectors[,200])
plot(blocked.vectors[,280], full.vectors[,280])
plot(blocked.vectors[,290], full.vectors[,290])

# It falls appart a bit with very small eigenvalues. I think this is not suprising
plot(blocked.vectors[,299], full.vectors[,299])

write.table(full.values, file=paste0(output.prefix, "/300G_r_calculated_full_decompose_eigenvalues.txt"), sep="\t", col.names = NA, quote=F)
write.table(full.vectors, file=paste0(output.prefix, "/300G_r_calculated_full_decompose_eigenvectors.txt"), sep="\t", col.names = NA, quote=F)

# ------------------------------------------------------------------------------
# invert some eigenvectors

subset           <- sample(1:300, 100)

full.vectors.inv <- full.vectors
full.vectors.inv[, subset] <-  full.vectors.inv[, subset] * -1
write.table(full.vectors.inv, file=paste0(output.prefix, "/300G_r_calculated_full_decompose_eigenvectors_random_33perc_inverted.txt"), sep="\t", col.names = NA, 
quote=F)

blocked.vectors.inv <- blocked.vectors
blocked.vectors.inv[, subset] <-  blocked.vectors.inv[, subset] * -1
write.table(blocked.vectors.inv, file=paste0(output.prefix, "/300G_r_calculated_blocked_eigenvectors_random_33perc_inverted.txt"), sep="\t", col.names = NA, 
quote=F)



