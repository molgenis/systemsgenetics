#remoter::server(verbose = T, port = 55557, password = "laberkak", sync = T)

remoter::client("localhost", port = 55557, password = "laberkak")

setwd("/groups/umcg-wijmenga/tmp04/projects/depict2/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

library(readr)


#if(FALSE){
#  #don't run by accident
#  table_tmp <- read_delim("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.txt.gz", delim = "\t", quote = "")
#  gnExp <- as.matrix(table_tmp[,-1])
#  rownames(gnExp) <- table_tmp[,1][[1]]
#  rm(table_tmp)
#  saveRDS(gnExp, "/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.rds")
#}
#gnExp <- readRDS("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.rds")

ensg75pc <- read.delim("ensgR75_protein_coding.txt", stringsAsFactors = F)

table_tmp <- read_delim("/groups/umcg-wijmenga/prm02/data_projects/Gado/GeneNetwork_V2_01-02-2018/Covariates/PCA/pc-scores250.txt.gz", delim = "\t", quote = "")
pcs <- as.matrix(table_tmp[,-1])
rownames(pcs) <- table_tmp[,1][[1]]
rm(table_tmp)
pcs <- pcs[,1:165]
str(pcs)

library(uwot)

init <- pcs[,c(2,1)]
init <- init * -1



 
#sampleUmap <- umap(pcs, n_threads = 24, n_epochs = 1000, init = init,  n_neighbors = 500, min_dist = 2, init_sdev = 1e-4, learning_rate = 1, spread = 40 ,scale = "none", nn_method = "fnn")
#sampleUmap <- umap(pcs, n_threads = 24, n_epochs = 10000, init = init,  n_neighbors = 500, min_dist = 2, init_sdev = 1e-4, learning_rate = 2, spread = 40 ,scale = "none", nn_method = "fnn")


#rownames(sampleUmap) <- rownames(pcs)
#colnames(sampleUmap) <- colnames(c("UMAP1", "UMAP2"))

#write.table(sampleUmap, file = "Umap/umap2.txt", sep = "\t", quote = F)


sampleUmap <- read.delim("Umap/umap2.txt")




tissueColMap <- read.delim("Umap/tissueCol2.txt", stringsAsFactors = F)
sampleAnnotation <- read.delim("Umap/sampleAnnotations.txt", stringsAsFactors = F, row.names = 1 )
str(sampleAnnotation)

all(row.names(sampleAnnotation)  %in% row.names(sampleUmap))

sampleAnnotation <- sampleAnnotation[match(row.names(sampleUmap), rownames(sampleAnnotation) ),]

all(row.names(sampleAnnotation) == row.names(sampleUmap))


sampleAnnotation$pch = 1
sampleAnnotation$pch[sampleAnnotation$CellLine != ""] = 5
sampleAnnotation$pch[sampleAnnotation$TissueType != ""] = 21

sampleAnnotation$col = adjustcolor("grey", alpha.f = 0.1)
#sampleAnnotation$col = NA
sampleAnnotation$bg = NA

for(i in 1:nrow(tissueColMap)){
  sampleAnnotation$col[sampleAnnotation$PlotClass == tissueColMap$PlotClass[i]] <- adjustcolor(tissueColMap$col[i], alpha.f = 0.1)
  sampleAnnotation$bg[sampleAnnotation$PlotClass == tissueColMap$PlotClass[i]] <- adjustcolor(tissueColMap$col[i], alpha.f = 0.3)
}

#sampleAnnotation$bg[sampleAnnotation$PlotClass == "pancreas"] <- adjustcolor("turquoise1", alpha.f = 0.5)
#sampleAnnotation$bg[sampleAnnotation$PlotClass == "osteoblastic cell"] <- adjustcolor("red", alpha.f = 0.5)

#rpng( width = 1200, height = 1000)
layout(matrix(2:1, nrow =1),widths = c(5,1))
par(mar = c(0,0,0,0))
plot.new()
plot.window(0:1,0:1)
legend(x = "center", legend = tissueColMap$PlotClass, fill = tissueColMap$col)
par(mar = c(3,3,1,0))
plot(sampleUmap, xlim = c(-370,390), ylim = c(-410,340), col = sampleAnnotation$col, pch = 21, bg = sampleAnnotation$bg, cex = 0.5)
#dev.off()
  
abline(h=-70)
abline(h=-100)
abline(v=-40)
abline(v=-60)

library(RColorBrewer)
display.brewer.all()
cat(RColorBrewer::brewer.pal(12, "Set3"), sep = "\n")



colfunc <- colorRampPalette(c("khaki2", "orangered", "brown4"))


cols = c("orangered1", "orchid1", "magenta", "indianred1", "hotpink2", "darkorange1", "darkorchid3", "deeppink4", "brown3", "darkred", "khaki2", "goldenrod1", "orange", "mediumorchid4")
length(cols)
plot(rep(1,14),col=cols,pch=19,cex=3)


cols = c("", "", "", "", "", "", "", "", "", "", "", "", "")


source(paste0("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\Downstreamer\\src\\main\\r\\downstreamer_main/downstreamer_functions.r"))
ced <- read.depict2("final_paper/multiple_sclerosis_2013_24076602_hg19_enrichtments_exHla.xlsx")

ced <- read.depict2("multiple_sclerosis_patsopoulos_harm_jan_enrichtments_exHla_1.xlsx")
names(ced)

cedExp <- ced$ 
cedExp <- cedExp[match(row.names(sampleUmap), cedExp$Gene.set),]

all(cedExp$Sample == row.names(sampleUmap))

str(cedExp)
#colCed <- brewer.pal(9, "GnBu")[as.numeric(cut(cedExp$Enrichment.Z.score[match(row.names(sampleUmap), cedExp$Sample)],breaks = 9))]

colCed <- rep("grey90", nrow(sampleUmap)) #lightblue1
#colCed[cedExp$Enrichment.Z.score > 0 & cedExp$FDR.5..significant] <- "orange"
#colCed[cedExp$Enrichment.Z.score > 0 & cedExp$Bonferroni.significant] <- "red4"

colfunc<-colorRampPalette(c("gold2","orange","orangered","red3"))
colMap <- colfunc(20)
colCed[cedExp$Enrichment.Z.score > 0 & cedExp$FDR.5..significant] <- colMap[as.numeric(cut(cedExp$Enrichment.Z.score[cedExp$Enrichment.Z.score > 0 & cedExp$FDR.5..significant],20))]


colfunc<-colorRampPalette(c("lightblue1","dodgerblue4"))
colMap <- colfunc(20)
colCed[cedExp$Enrichment.Z.score < 0 & cedExp$FDR.5..significant] <- colMap[as.numeric(cut(cedExp$Enrichment.Z.score[cedExp$Enrichment.Z.score < 0 & cedExp$FDR.5..significant],20))]


plot(sampleUmap, xlim = c(-370,390), ylim = c(-410,340), col = adjustcolor(colCed, alpha.f = 0.5), pch = 21, bg = adjustcolor(colCed, alpha.f = 0.1), cex = 0.5)


