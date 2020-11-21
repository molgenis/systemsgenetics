#remoter::server(verbose = T, port = 55557, password = "laberkak", sync = T)

remoter::client("localhost", port = 55557, password = "laberkak")

setwd("/groups/umcg-wijmenga/tmp04/projects/depict2/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

library(readr)
library(gatepoints)

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




tissueColMap <- read.delim("Umap/tissueCol5.txt", stringsAsFactors = F)
sampleAnnotation <- read.delim("Umap/sampleAnnotations.txt", stringsAsFactors = F, row.names = 1 )
str(sampleAnnotation)

all(row.names(sampleAnnotation)  %in% row.names(sampleUmap))

sampleAnnotation <- sampleAnnotation[match(row.names(sampleUmap), rownames(sampleAnnotation) ),]

all(row.names(sampleAnnotation) == row.names(sampleUmap))


#sampleAnnotation$pch = 1
#sampleAnnotation$pch[sampleAnnotation$CellLine != ""] = 5
#sampleAnnotation$pch[sampleAnnotation$TissueType != ""] = 21

defaultCol <- "gray90"

sampleAnnotation$col = adjustcolor(defaultCol, alpha.f = 0.1)
#sampleAnnotation$col = NA
#sampleAnnotation$bg = NA



meanX <- sapply(tissueColMap$PlotClass, function(plotClass){
  mean(sampleUmap[sampleAnnotation$PlotClass == plotClass,1])
})

meanY <- sapply(tissueColMap$PlotClass, function(plotClass){
  mean(sampleUmap[sampleAnnotation$PlotClass == plotClass,2])
})

#clusterLabels <- data.frame("PlotClass" = tissueColMap$PlotClass, "centerX" = meanX, "centerY" = meanY, "offsetX" = "", "offsetY" = "", "label" = tissueColMap$PlotClass)
#write.table(clusterLabels, "Umap/lables.txt", sep = "\t", quote = F, row.names = F)


str(clusterLabels)


#sampleAnnotation$bg[sampleAnnotation$PlotClass == "osteoblastic cell"] <- adjustcolor("red", alpha.f = 0.5)


sampleAnnotation$plotOrder <- 1
sampleAnnotation$plotOrder[sampleAnnotation$col != defaultCol] <- 2


sampleAnnotation <- sampleAnnotation[ order(sampleAnnotation$plotOrder),]
sampleUmap <- sampleUmap[order(sampleAnnotation$plotOrder),]



clusterLabels <- read.delim("Umap/lables.txt",stringsAsFactors = F, comment.char = "#")


for(i in 1:nrow(tissueColMap)){
  sampleAnnotation$col[sampleAnnotation$PlotClass == tissueColMap$PlotClass[i]] <- adjustcolor(tissueColMap$col[i], alpha.f = 0.1)
  sampleAnnotation$bg[sampleAnnotation$PlotClass == tissueColMap$PlotClass[i]] <- adjustcolor(tissueColMap$col[i], alpha.f = 0.3)
}

sampleAnnotation$col[sampleAnnotation$PlotClass == "spleen"] <- adjustcolor("black", alpha.f = 0.5)


#X11()
#rpng( width = 1200, height = 1000)
#pdf("Umap/sampleUmap.pdf", width = 10, height = 10, useDingbats = F)
layout(matrix(2:1, nrow =1),widths = c(5,1))
par(mar = c(0,0,0,0), xpd = NA)
plot.new()
plot.window(0:1,0:1)
legend(x = "center", legend = tissueColMap$PlotClass, fill = tissueColMap$col)
par(mar = c(3,3,1,0))
plot.new()
plot.window(xlim = c(-340,390), ylim = c(-380,390))
axis(side = 1, at = c(-200,0,200), col = "gray30", lwd = 1.5, col.axis = "gray30")
axis(side = 2, at = c(-200,0,200), col = "gray30", lwd = 1.5, col.axis = "gray30")
points(sampleUmap, col = sampleAnnotation$col, pch = 16, cex = 0.6)

text(clusterLabels$centerX + clusterLabels$offsetX, clusterLabels$centerY + clusterLabels$offsetY, labels = clusterLabels$label, col = "gray30")


#dev.off()

fhs(sampleUmap)

library(grid)
grid.locator(unit = "npc")
  
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
str(ced)

cedExp <- ced$Expression
cedExp <- cedExp[match(row.names(sampleUmap), cedExp$Sample),]

all(cedExp$Sample == row.names(sampleUmap))

str(cedExp)
#colCed <- brewer.pal(9, "GnBu")[as.numeric(cut(cedExp$Enrichment.Z.score[match(row.names(sampleUmap), cedExp$Sample)],breaks = 9))]

colCed <- rep("grey90", nrow(sampleUmap)) #lightblue1
#colCed[cedExp$Enrichment.Z.score > 0 & cedExp$FDR.5..significant] <- "orange"
#colCed[cedExp$Enrichment.Z.score > 0 & cedExp$Bonferroni.significant] <- "red4"


fdrZthreshold <- min(abs(cedExp$Enrichment.Z.score[cedExp$FDR.5..significant]))
maxZ <- max(cedExp$Enrichment.Z.score)
breaks <- seq(fdrZthreshold, maxZ, length.out = 21)

greyBreakCount <- round((fdrZthreshold - -fdrZthreshold ) / (breaks[2] - breaks[1]))

colfunc<-colorRampPalette(c("gold2","orange","orangered","red3"))
colMapEnrich <- colfunc(20)
colCed[cedExp$Enrichment.Z.score > 0 & cedExp$FDR.5..significant] <- colMapEnrich[as.numeric(cut(cedExp$Enrichment.Z.score[cedExp$Enrichment.Z.score > 0 & cedExp$FDR.5..significant],breaks))]


colfunc<-colorRampPalette(c("lightblue1","dodgerblue4"))
colMapDepl <- colfunc(20)
colCed[cedExp$Enrichment.Z.score < 0 & cedExp$FDR.5..significant] <- colMapDepl[as.numeric(cut(abs(cedExp$Enrichment.Z.score[cedExp$Enrichment.Z.score < 0 & cedExp$FDR.5..significant]),breaks))]

fullColGradient <- c(rev(colMapDepl), rep("grey90", greyBreakCount), colMapEnrich)

library( plotfunctions)

pdf("MS-UMAP.pdf", width = 10, height = 10, useDingbats = F)
plot(sampleUmap[order(abs(cedExp$Enrichment.Z.score),decreasing = F),], xlim = c(-380,390), ylim = c(-410,390), col = adjustcolor(colCed[order(abs(cedExp$Enrichment.Z.score),decreasing = F)], alpha.f = 0.5), pch = 16, cex = 0.4, xlab = "UMAP-1", ylab = "UMAP-2", bty = "n")#, bg = adjustcolor(colCed[order(abs(cedExp$Enrichment.Z.score),decreasing = F)], alpha.f = 0.05)

gradientLegend(
  c(-maxZ,maxZ),
  color = fullColGradient,
  pos = 0.6,
  side = 3,
  dec = 1,
  length = 0.4,
  depth = 0.025,
  inside = T,
  coords = FALSE,
  pos.num = NULL,
  n.seg = c(-maxZ,-bonfZscore,-fdrZthreshold,0,fdrZthreshold,bonfZscore,maxZ),
  border.col = "black",
  tick.col = NULL,
  fit.margin = TRUE,
)



dev.off()


