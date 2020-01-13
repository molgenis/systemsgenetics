
patrick.old <- list()
patrick.old$data <- read.table("~/Desktop/depict2/topGeneP.txt")
tmp <- read.table("~/Desktop/coreg_edu_att_maf_filtred.csv", sep="\t", header=T)
patrick.old$coregulation <- tmp$Enrichment.Z.score
names(patrick.old$coregulation) <- tmp$Gene.set
network <- depict.to.network(patrick.old)