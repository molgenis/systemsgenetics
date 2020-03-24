#library(tidygraph)
#library(ggraph)
library(data.table)

# IO args
args           <- commandArgs(trailingOnly=TRUE)
bundle.dir     <- args[4]
input.path     <- args[1]
phenotype      <- args[2]
output.path    <- args[3]
input.prefix   <- paste0(input.path, "/", phenotype)

# Args
edge.threshold <- 2
edges.within.cols.and.rows <- T

cat("Following args have been supplied:\n")
cat(paste0("\tbundle dir:\t",     bundle.dir, "\n"))
cat(paste0("\tinput path:\t",     input.path, "\n"))
cat(paste0("\tphenotype:\t",      phenotype, "\n"))
cat(paste0("\toutput path:\t",    output.path, "\n"))
cat(paste0("\tedge threshold:\t", edge.threshold, "\n"))
cat(paste0("\tmake all edges:\t", edges.within.cols.and.rows, "\n"))

#------------------------------------------------------------------------------
# Read data
#------------------------------------------------------------------------------
# Gene names
ensembl <- read.table(paste0(bundle.dir, "/reference_datasets/human_b37/ensembl_gene_position_export.txt"), sep="\t", header=T, stringsAsFactors = F)
ensembl <- unique(ensembl[,c(1, 5)])
rownames(ensembl) <- ensembl[,1]

# Read genes
rows <- as.character(read.table(paste0(input.prefix, ".rows.txt"), stringsAsFactors = F)[,1])
cols <- as.character(read.table(paste0(input.prefix, ".cols.txt"), stringsAsFactors = F)[,1])

# Read coreg matrix
coreg           <- fread(paste0(input.path, "/", phenotype, "_coreg_genes.txt"), data.table = F)
rownames(coreg) <- coreg[,1]
cur.data        <- list(full.coregulation.matrix=coreg, rows=rows, cols=cols)

#------------------------------------------------------------------------------
# Function definitions
#------------------------------------------------------------------------------
# Determines the type of the edge, returns string
determine.edge.type <- function(gene1, gene2, rows, cols) {
  # Not very elegant code, but it does the job
  if (gene1 %in% rows & gene2 %in% rows & !gene1 %in% cols & !gene2 %in% cols) {
    return("within_row")
  }
  
  if (gene1 %in% cols & gene2 %in% cols & !gene1 %in% rows & !gene2 %in% rows) {
    return("within_col")
  }
  
  if (gene1 %in% rows & gene2 %in% cols & !gene1 %in% cols) {
    return("row_to_col")
  }
  
  if (gene2 %in% rows & gene1 %in% cols & !gene2 %in% cols) {
    return("row_to_col")
  }
  
  if (gene1 %in% rows & gene2 %in% cols & gene1 %in% cols) {
    return("row_and_col")
  }
  
  if (gene2 %in% rows & gene1 %in% cols & gene2 %in% cols) {
    return("row_and_col")
  }
  
  return("NA")
}

# Constructs network, returns list with nodes df and edges df
depict.to.network <- function(cur.data, edge.threshold=2, edges.within.cols.and.rows=T) {
  rows  <- cur.data$rows
  cols  <- cur.data$cols
  coreg <- cur.data$full.coregulation.matrix
 
  if (edges.within.cols.and.rows) {
    # If so X is a square matrix
    coreg <- coreg[unique(c(rows, cols)), unique(c(rows, cols))]
  } else {
    coreg <- coreg[rows, cols]
  }
  
  # Construct the nodes dataframe
  nodes                            <- as.data.frame(unique(c(rows, cols)), stringsAsFactors = F)
  # Add the source of the node
  nodes$annot                      <- rep("", nrow(nodes))
  nodes$annot[nodes[,1] %in% rows] <- "rows"
  nodes$annot[nodes[,1] %in% cols] <- "cols" 
  nodes$annot[nodes[,1] %in% rows & nodes[,1] %in% cols] <- "both"
  
  # Add other annotations and make sure characters are not factors
  colnames(nodes)                  <- c("gene_id", "annot")
  rownames(nodes)                  <- nodes$gene_id
  nodes$annot                      <- as.character(nodes$annot)
  nodes$gene_name                  <- as.character(ensembl[nodes$gene_id, 2])
  #nodes$zscore                     <- cur.data$coregulation[nodes$gene_id]
  #nodes$gene_p
  nodes                            <- nodes[order(nodes$annot, decreasing = T),]

  # Construct edges at zscore threshold
  tmp.edges       <- abs(coreg) >= edge.threshold
  edges           <- as.data.frame(matrix(nrow=1, ncol=4))
  colnames(edges) <- c("from", "to", "zscore", "edge_type")
  
  # Remove duplicated edges for nodes which appear in both columns and rows
  ol <- intersect(rownames(tmp.edges), colnames(tmp.edges))
  tmp.edges[ol, ol][upper.tri(tmp.edges[ol, ol])] <- F
  
  # Construct the edges dataframe
  for (row in 1:nrow(tmp.edges)) {
    for(col in 1:ncol(tmp.edges)) {
      if (tmp.edges[row, col]) {
        rowgene <- rownames(coreg)[row]
        colgene <- colnames(coreg)[col]
        edges   <- rbind(edges, c(
          which(nodes$gene_id == rowgene),
          which(nodes$gene_id == colgene),
          coreg[row, col],
          determine.edge.type(rowgene, colgene, rows, cols)))
      }
    }
  }
  edges            <- na.omit(edges)
  edges$from       <- as.numeric(edges$from)
  edges$to         <- as.numeric(edges$to)
  
  named.edges      <- edges
  named.edges$from <- rownames(nodes)[edges$from]
  named.edges$to   <- rownames(nodes)[edges$to]

  return(list(nodes=nodes, edges=edges, named.edges=named.edges))
}

#------------------------------------------------------------------------------
# Run depict.to.network
#------------------------------------------------------------------------------
network <- depict.to.network(cur.data, edge.threshold=edge.threshold, edges.within.cols.and.rows=edges.within.cols.and.rows)

# Write nodes and edges files
write.table(network$nodes,file=paste0(output.path, "/", phenotype, "_nodes.tsv"), sep="\t", row.names=F, quote=F)
write.table(network$edges, file=paste0(output.path, "/", phenotype, "_edges.tsv"), sep="\t", quote=F, row.names=F)
write.table(network$named.edges, file=paste0(output.path, "/", phenotype, "_named_edges.tsv"), sep="\t", quote=F, row.names=F)

