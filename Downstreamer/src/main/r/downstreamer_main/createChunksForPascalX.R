
#remoter::server(verbose = T, port = 55556, sync = T)

library(pbdZMQ)
sessionInfo()

remoter::client("localhost", port = 55556)#55501


genes <- read.delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/genes_Ensembl94.txt.gz", check.names = F)
str(genes)



#Ensembl Gene ID          Chromosome Name          Gene Start (bp)          Gene End (bp)            Gene Biotype             Band

genes <- genes[,c("Gene stable ID", "Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)", "Gene type", "Karyotype band")]

str(genes[101:200,])

write.table(genes, gzfile("/groups/umcg-fg/tmp01/projects/downstreamer/PascalX_bundle/genes_Ensembl94.txt.gz"), sep = "\t", quote = FALSE, row.names = F)


chr = "1"
chunk = 1

genesPerChunk <- 100

chrChunks <- sapply(as.character(1:22), function(chr){
  genesChr <- genes[genes$`Chromosome/scaffold name` == chr,]
  chunks <- ceiling(nrow(genesChr) / genesPerChunk)
  
  for(chunk in 1:chunks){
    startChunk <- ((chunk -1) * 100) + 1
    endChunk <- startChunk + 99
    if(endChunk > nrow(genesChr)){
      endChunk <- nrow(genesChr)
    }
    
    write.table(genesChr[startChunk:endChunk,], paste0("/groups/umcg-fg/tmp01/projects/downstreamer/PascalX_bundle/Ensembl94Chunks/chr",chr,"_chunk",chunk,".txt"), sep = "\t", quote = FALSE, row.names = F)
    

  }
  
  return(c(chr, chunk))
  
})
chrChunks <- t(chrChunks)
colnames(chrChunks) <- c("Chr", "Chunks")

write.table(chrChunks, "/groups/umcg-fg/tmp01/projects/downstreamer/PascalX_bundle/genes_Ensembl94_chunks.txt", sep = "\t", quote = FALSE, row.names = F)
