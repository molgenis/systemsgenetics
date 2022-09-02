#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)


remoter::client("localhost", port = 55501, password = "laberkak")



library(uwot)


setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

tissueCol <- read.delim("umap/col.txt", row.names = 1, na.strings = "")

#save(pcs, combinedMeta, tissueCol, file = "umap/dataForUmap.RData")
load(file = "DataForPredictions.RData")

#load(file = "combinedMeta_2022_08_30.RData", verbose = T)
#str(combinedMeta)
#updatedAnnotations <- combinedMeta[,c("Tissue",	"Tissue2",	"Cellline",	"CelllineName",	"Cancer",	"Cohort", "Fetal")]

#all(rownames(pcsAndMeta) %in% rownames(updatedAnnotations))
#updatedAnnotations <- updatedAnnotations[rownames(pcsAndMeta),]
#all(rownames(pcsAndMeta) == rownames(updatedAnnotations))

#pcsAndMeta[,colnames(updatedAnnotations)] <- updatedAnnotations

#pcsAndMeta$selectedSamples <- !pcsAndMeta$excludeBasedOnPredictionCellline2 & !pcsAndMeta$excludeBasedOnPredictionCancer & !(!is.na(pcsAndMeta$Cancer) & pcsAndMeta$Cancer) & !(!is.na(pcsAndMeta$Cellline) & pcsAndMeta$Cellline)

table(pcsAndMeta$selectedSamples, useNA = "a")


clusterAnnotations <- read.delim("umap/annotationsBasedOnOldUmap.txt", row.names = 1)
pcsAndMeta <- merge(pcsAndMeta, clusterAnnotations, by = 0, all.x = T)
table(pcsAndMeta$ClusterAnnotation)




#pcsAndMeta[!is.na(pcsAndMeta$study) & (pcsAndMeta$study== "ERP104864") & (grepl("synovium", pcsAndMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= ""


tissueSamples <- pcsAndMeta[pcsAndMeta$selectedSamples,]


tissueSamples$class <- tissueSamples$Tissue

hasT2 <- tissueSamples$Tissue2 != ""
tissueSamples$class[hasT2] <- paste0(tissueSamples$class[hasT2], "-", tissueSamples$Tissue2[hasT2])

isFetal <- !is.na(tissueSamples$Fetal) & tissueSamples$Fetal
tissueSamples$class[isFetal] <- paste0(tissueSamples$class[isFetal], "-Fetal")

noTbutCluster <- tissueSamples$class == "" & !is.na(tissueSamples$ClusterAnnotation)
table(noTbutCluster, useNA = "a")
tissueSamples$class[noTbutCluster] <- tissueSamples$ClusterAnnotation[noTbutCluster]

table(tissueSamples$class)
write.table(table(tissueSamples$class, useNA = "always"), file = "umap/tissues.txt", sep = "\t", quote = F, row.names = F)

str(tissueSamples)



mapping <- read.delim("umap/tissuesMapping.txt")
str(mapping)

all(tissueSamples$class %in% mapping$Class)


tissueSamples$umapFactor <- as.factor(mapping$ClassificationClass[match(tissueSamples$class, mapping$Class)])

table(tissueSamples$umapFactor, useNA = "always")


defaultCol <- adjustcolor("grey", alpha.f = 0.6)
tissueCol <- read.delim("umap/col.txt", row.names = 1)


tissueSamples$TissueCol <- defaultCol
sum(unique(tissueSamples$umapFactor) %in% rownames(tissueCol))
tissueSamples$TissueCol[tissueSamples$umapFactor %in% rownames(tissueCol)] <- adjustcolor(tissueCol[as.character(tissueSamples$umapFactor[tissueSamples$umapFactor %in% rownames(tissueCol)]),1], alpha.f = 0.5)
#tissueSamples$TissueCol[tissueSamples$umapFactor %in% rownames(tissueCol)] <- tissueCol[as.character(tissueSamples$umapFactor[tissueSamples$umapFactor %in% rownames(tissueCol)]),1]
table(tissueSamples$TissueCol)

plotOrderTissues <- order((tissueSamples$TissueCol != defaultCol) + 1)

#, n_threads = 22

compsToUseForUmap <- 100
init <- as.matrix(tissueSamples[,paste0("PC_",1:2)])
umapInput <- as.matrix(tissueSamples[,paste0("PC_",1:compsToUseForUmap)])


sampleUmap <- umap(
  umapInput, 
  n_epochs = 500, 
  init = init, 
  n_neighbors = 500, 
  min_dist = 2, init_sdev = 1e-4, learning_rate = 1, spread = 40 ,scale = "scale", nn_method = "fnn",
  local_connectivity = 10)

rownames(sampleUmap) <- rownames(tissueSamples)
colnames(sampleUmap) <- c("UMAP1", "UMAP2")
umapAndMeta <- merge(sampleUmap, tissueSamples, by = 0)
dim(umapAndMeta)


rpng()

par(mar = c(3,5,0.1,0.1), xpd = NA)
plot(umapAndMeta[plotOrderTissues,"UMAP1"], umapAndMeta[plotOrderTissues,"UMAP2"], col = umapAndMeta$TissueCol[plotOrderTissues], cex = 0.2, pch = 16)

dev.off()


locator(n =2, type = "l")




write.table(umapAndMeta,file = "umaptest.txt", sep = "\t", quote = F, col.names = NA)

#save.image( file="umap_tmp.RData")
load("umap_tmp.RData") 

rpng()

par(mar = c(3,5,0.1,0.1), xpd = NA)
plot(umapAndMeta[plotOrder,"UMAP1"], umapAndMeta[plotOrder,"UMAP2"], col = umapAndMeta$TissueCol[plotOrder], cex = 0.8, pch = 16, xlim = c(-25,25), ylim = c(-25,25))

dev.off()



#png(file = "umaptest.png", width = 1600, height = 800)

pdf(file = "umaptest.pdf", width = 16, height = 8)
#rpng()

layout(matrix(1:2,ncol = 2))

par(mar = c(3,5,0.1,0.1), xpd = NA)
plot(umapAndMeta[plotOrder,"UMAP1"], umapAndMeta[plotOrder,"UMAP2"], col = umapAndMeta$TissueCol[plotOrder], cex = 0.8, pch = 16)

par(mar = c(0,0,0,0), xpd = NA)
plot.new()
plot.window(xlim = 0:1, ylim = 0:1)
legend("center", fill = tissueCol[,1], legend = row.names(tissueCol), bty = "n", ncol = 2)


dev.off()





