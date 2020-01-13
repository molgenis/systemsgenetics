library(tidygraph)
library(ggraph)
library(data.table)


#path           <- "/home/work/Desktop/depict2/coregulation_data/HP:0000991"
#phenotype      <- "30690_raw"

#path           <- "/home/work/Desktop/depict2/coregulation_data/HP:0003124/"
#phenotype      <- "30690_raw"

path           <- "/home/work/Desktop/depict2/coregulation_data/HP:0011153/"
phenotype      <- "educational_attainment_2018_30038396_hg19"


path           <- "/home/work/Desktop/depict2/coregulation_data/alzheimers_top10/"
phenotype      <- "alzheimers_disease_2018_29777097_hg19"

ngene          <- 100
edge.threshold <- 2

ensembl <- read.table("~/Documents/data/reference/ensembl/ensembl_gene_position_export.txt", sep="\t", header=T, stringsAsFactors = F)
ensembl <- unique(ensembl[,c(1, 5)])
rownames(ensembl) <- ensembl[,1]

read.coreg.data <- function(path, phenotype, ngene=100) {
  
  out <- list()
  genep.tmp       <- fread(paste0(path, "/", phenotype, "_gene_pvalues.txt"), data.table = F)
  genep           <- genep.tmp[,2]
  names(genep)    <- genep.tmp[,1]
  genep           <- genep[order(genep)]
  
  tmp <- sum(genep == min(genep, na.rm=T), na.rm=T)
  if (tmp >= ngene) {
    print(paste0("[WARN] Detected ", tmp, " ties"))
    ngene <- tmp
  }

  genes.to.keep   <- names(genep[1:ngene])
  coreg           <- fread(paste0(path, "/", phenotype, "_coreg_hpo_genes.txt"), data.table = F)
  rownames(coreg) <- coreg[,1]
  coreg.out       <- coreg[genes.to.keep,-1]
  
  coreg.tmp       <- fread(paste0(path, "/", phenotype, "_coregulation.txt"), data.table = F)
  crg             <- coreg.tmp[,-1]
  names(crg)      <- coreg.tmp[,1]

  out$gene.pvalues <- genep
  out$coregulation <- crg
  out$data         <- coreg.out
  out$full.matrix  <- coreg
  
  return(out)
}

depict.to.network <- function(cur.data, edge.threshold=2) {
  out             <- list()
  cur.df          <- cur.data$data
  nodes           <- as.data.frame(c(colnames(cur.df), rownames(cur.df)), stringsAsFactors = F)
  nodes           <- cbind(nodes, c(rep("Mendelian gene", ncol(cur.df)), rep("GWAS gene", nrow(cur.df))))
  colnames(nodes) <- c("gene_id", "annot")
  nodes$annot     <- as.character(nodes$annot)
  nodes$gene_name <- as.character(c(ensembl[colnames(cur.df), 2], rep(NA, nrow(cur.df))))
  nodes$zscore    <- cur.data$coregulation[nodes$gene_id]

  # Set negative coregulation to zero
  nodes[nodes$zscore < 0, ]$zscore <- 0
  
  # Remove duplicated nodes from node list, keeping the first (mendelian gene)
  nodes           <- nodes[!duplicated(nodes$gene_id),]
  rownames(nodes) <- nodes$gene_id
  nodes[intersect(colnames(cur.df), rownames(cur.df)), "annot"] <- "Both"
  nodes           <- nodes[order(nodes$annot, decreasing = T),]
  
  # Construct edges at zscore threshold
  tmp.edges       <- cur.df >= edge.threshold
  edges           <- as.data.frame(matrix(nrow=1, ncol=3))
  colnames(edges) <- c("from", "to", "zscore")
  
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
  
  edges <- na.omit(edges)
  #nodes <- nodes[unique(c(edges[,1], edges[,2])), ,drop=F]
  
  top.zscore.genes <- unique(nodes[edges[order(edges$zscore, decreasing = T),1][1:10], "gene_id"])
  nodes[top.zscore.genes, "gene_name"] <- ensembl[top.zscore.genes, 2]
  
  
  tbl   <- tbl_graph(nodes=nodes, edges=edges, directed=F)
  return(list(graph=tbl, nodes=nodes, edges=edges))
}

cur.data <- read.coreg.data(path, phenotype, ngene=ngene)
network  <- depict.to.network(cur.data, edge.threshold=edge.threshold)

# Make plots
pdf(width=12, height=10, useDingbats = F, file=paste0("~/Desktop/depict2/", phenotype, "_edges_", edge.threshold, "_network_plot.pdf"))
par(xpd=NA)
p <- ggraph(network$graph,layout="linear", circular=T) +
  geom_edge_arc(aes(width = zscore, alpha=zscore))+
  geom_node_point(aes(colour=annot, size=zscore)) +
  theme_graph(base_family = 'sans') +
  scale_size(range=c(1,5), name="Depict2 (zscore)") +
  scale_edge_width(range = c(0.5, 2), name="Coregulation (Zscore)") +
  scale_edge_alpha(name="Coregulation (Zscore)") +
  scale_color_manual(name="Gene type", values=c(`Both`="#FFBF40",
                                                `GWAS gene`="#33CCCC",
                                                `Mendelian gene`="#FF4040")) +
  guides(colour = guide_legend(override.aes = list(size=5)))

p
p + geom_node_label(aes(label=gene_name, fill=annot),
                    colour="white",
                    show.legend = F,
                    label.size = 0,
                    repel=T, segment.colour="black") +
  scale_fill_manual(values=c(`Both`="#FFBF40",
                             `GWAS gene`="#33CCCC",
                             `Mendelian gene`="#FF4040"))

# Remove nodes without an edge
sub_mygraph <- to_subgraph(network$graph, gene_id %in% network$nodes[unique(c(network$edges[,1], network$edges[,2])), 1], subset_by = "nodes")$subgraph

p2 <- ggraph(sub_mygraph, layout="linear", circular=T)  +
  geom_edge_arc(aes(width = zscore, alpha=zscore)) +
  geom_node_point(aes(colour=annot, size=zscore)) +
  theme_graph(base_family = 'sans') +
  scale_size(range=c(1,5), name="Depict2 (zscore)") +
  scale_edge_width(range = c(0.5, 2), name="Coregulation (Zscore)") +
  scale_edge_alpha(name="Coregulation (Zscore)") +
  scale_color_manual(name="Gene type", values=c(`Both`="#FFBF40",
                                                  `GWAS gene`="#33CCCC",
                                                  `Mendelian gene`="#FF4040")) +
  guides(colour = guide_legend(override.aes = list(size=5)))

p2
p2 + geom_node_label(aes(label=gene_name, fill=annot),
                     colour="white",
                     show.legend = F,
                     label.size = 0,
                     repel=F) +
  scale_fill_manual(values=c(`Both`="#FFBF40",
                              `GWAS gene`="#33CCCC",
                              `Mendelian gene`="#FF4040"))

p3 <- ggraph(sub_mygraph, layout="gem") +
  geom_edge_link(aes(width = zscore, alpha=zscore)) +
  geom_node_point(aes(colour=annot, size=zscore)) +
  theme_graph(base_family = 'sans') +
  scale_size(range=c(1,5), name="Depict2 (zscore)") +
  scale_edge_width(range = c(0.5, 2), name="Coregulation (Zscore)") +
  scale_edge_alpha(name="Coregulation (Zscore)") +
  scale_color_manual(name="Gene type", values=c(`Both`="#FFBF40",
                                                  `GWAS gene`="#33CCCC",
                                                  `Mendelian gene`="#FF4040")) +
  guides(colour = guide_legend(override.aes = list(size=5)))


p3
p3 + geom_node_label(aes(label=gene_name, fill=annot),
                     colour="white",
                     show.legend = F,
                     label.size = 0,
                     repel=F) +
  scale_fill_manual(values=c(`Both`="#FFBF40",
                             `GWAS gene`="#33CCCC",
                             `Mendelian gene`="#FF4040"))


dev.off()

