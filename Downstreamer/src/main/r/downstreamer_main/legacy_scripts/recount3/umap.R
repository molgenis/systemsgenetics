#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)


remoter::client("localhost", port = 55501, password = "laberkak")



library(uwot)


setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

tissueCol <- read.delim("umap/col.txt", row.names = 1, na.strings = "")

#save(pcs, combinedMeta, tissueCol, file = "umap/dataForUmap.RData")
load(file = "DataForLasso.RData")

load(file = "combinedMeta_2022_08_30.RData", verbose = T)
str(combinedMeta)
updatedAnnotations <- combinedMeta[,c("Tissue",	"Tissue2",	"Cellline",	"CelllineName",	"Cancer",	"Cohort", "Fetal")]

all(rownames(pcsAndMeta) %in% rownames(updatedAnnotations))
updatedAnnotations <- updatedAnnotations[rownames(pcsAndMeta),]
all(rownames(pcsAndMeta) == rownames(updatedAnnotations))

pcsAndMeta[,colnames(updatedAnnotations)] <- updatedAnnotations

pcsAndMeta$selectedSamples <- !pcsAndMeta$excludeBasedOnPredictionCellline2 & !pcsAndMeta$excludeBasedOnPredictionCancer & !(!is.na(pcsAndMeta$Cancer) & pcsAndMeta$Cancer) & !(!is.na(pcsAndMeta$Cellline) & pcsAndMeta$Cellline)

table(pcsAndMeta$selectedSamples, useNA = "a")



  
  
tissueSamples <- pcsAndMeta[pcsAndMeta$selectedSamples,]
str(tissueSamples)

init <- as.matrix(tissueSamples[,paste0("PC_",1:2)])

mapping <- read.delim("umap/tissuesMapping.txt")
str(mapping)

tissueSamples$TissueCombined <- paste0(tissueSamples$Tissue, "-", tissueSamples$Tissue2)
all(tissueSamples$TissueCombined %in% mapping$Tissue)

tissueSamples$umapFactor <- as.factor(mapping$ClassificationClass[match(tissueSamples$TissueCombined, mapping$Tissue)])

table(tissueSamples$umapFactor, useNA = "always")


#write.table(table(tissueSamples$umapFactor, useNA = "always"), file = "tissues.txt", sep = "\t", quote = F)

#, n_threads = 22

compsToUseForUmap <- 100

umapInput <- as.matrix(tissueSamples[,paste0("PC_",1:compsToUseForUmap)])





rpng()

par(mar = c(3,5,0.1,0.1), xpd = NA)
plot(umapAndMeta[plotOrder,"UMAP1"], umapAndMeta[plotOrder,"UMAP2"], col = umapAndMeta$TissueCol[plotOrder], cex = 0.5, pch = 16)

dev.off()

umapAndMeta$TissueCol[umapAndMeta$umapFactor == "Muscle"]

write.table(umapAndMeta,file = "umaptest.txt", sep = "\t", quote = F, col.names = NA)

save.image( file="umap_tmp.RData")

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





