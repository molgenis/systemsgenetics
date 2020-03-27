
# For ensembl readnames
ensembl <- read.table("~/Documents/data/reference/ensembl/ensembl_gene_position_export.txt", sep="\t", header=T, stringsAsFactors = F)
ensembl <- unique(ensembl[,c(1, 5)])
rownames(ensembl) <- ensembl[,1]

# Read coregulation matrix
cur.df <- read.table("~/Desktop/gene_network_niersteen_filtered_for_patrick_2020-02-07.txt", row.names=1, header=T)
cur.df <- read.table("~/Desktop/kidney_network_gene_gene_cor_zcores.txt", row.names=1, header=T)

edge.threshold <- -100

# Make node matrix
nodes           <- as.data.frame(c(colnames(cur.df), rownames(cur.df)), stringsAsFactors = F)
nodes$annot     <- rep("", nrow(nodes))
colnames(nodes) <- c("gene_id", "annot")
nodes$annot     <- as.character(nodes$annot)
nodes$gene_name <- as.character(ensembl[nodes$gene_id, 2])
nodes           <- nodes[!duplicated(nodes$gene_id),]
rownames(nodes) <- nodes$gene_id

# Construct edges based on zscore threshold
tmp.edges       <- cur.df >= edge.threshold
edges           <- as.data.frame(matrix(nrow=1, ncol=3))
colnames(edges) <- c("from", "to", "zscore")

# Remove duplicated edges for nodes which appear in both columns and rows
ol <- intersect(rownames(cur.df), colnames(cur.df))
tmp.edges[ol, ol][upper.tri(tmp.edges[ol, ol])] <- F
diag(tmp.edges) <- F

for (row in 1:nrow(tmp.edges)) {
  for(col in 1:ncol(tmp.edges)) {
    if (tmp.edges[row, col]) {
      #edges <- rbind(edges, c((ncol(tmp.edges) + row), col, cur.df[row, col]))
      rowgene <- rownames(cur.df)[row]
      colgene <- colnames(cur.df)[col]
      edges   <- rbind(edges, c(
        which(nodes$gene_id == rowgene),
        which(nodes$gene_id == colgene),
        cur.df[row, col]))
    }
  }
}

# Replace edge indices with gene id's
edges      <- na.omit(edges)
edges$from <- rownames(nodes)[edges$from]
edges$to   <- rownames(nodes)[edges$to]


write.table(nodes, file="~/Desktop/depict2/niersteen_kidney_network_nodes.tsv", sep="\t", row.names=F, quote=F)
write.table(edges, file="~/Desktop/depict2/niersteen_kidney_network_edges.tsv", sep="\t", quote=F, row.names=F)

write.table(nodes, file="~/Desktop/depict2/niersteen_gene_network_nodes.tsv", sep="\t", row.names=F, quote=F)
write.table(edges, file="~/Desktop/depict2/niersteen_gene_network_edges.tsv", sep="\t", quote=F, row.names=F)

