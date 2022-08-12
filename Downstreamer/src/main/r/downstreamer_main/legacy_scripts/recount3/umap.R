library(uwot)


#tissueCol <- read.delim("Recount3_QC_2ndRun/SRA_Studies_Annotations_Patrick/Annotations_color2.txt", row.names = 1)

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")
#save(pcs, combinedMeta, tissueCol, file = "umap/dataForUmap.RData")
load(file = "umap/dataForUmap.RData")

init <- pcs[,c(1,2)]

table(combinedMeta$CelllineName)

table(combinedMeta$Tissue)

combinedMeta$umapFactor <- paste0(combinedMeta$Tissue, "-", combinedMeta$Tissue2)

table(combinedMeta$umapFactor, useNA = "always")



dim(pcs)
sampleUmap <- umap(pcs[,1:200], n_threads = 22, n_epochs = 100000, init = init,  n_neighbors = 500, min_dist = 2, init_sdev = 1e-4, learning_rate = 1, spread = 40 ,scale = "none", nn_method = "fnn")


rownames(sampleUmap) <- rownames(pcs)
colnames(sampleUmap) <- c("UMAP1", "UMAP2")
str(sampleUmap)
umapAndMeta <- merge(sampleUmap, combinedMeta, by = 0, all.x = T)
dim(umapAndMeta)



defaultCol <- adjustcolor("grey", alpha.f = 0.6)
umapAndMeta$col <- defaultCol

tissueAndCol <- umapAndMeta[,"Tissue"] %in% tissueCol$PlotClass
umapAndMeta$col[tissueAndCol] <-  adjustcolor(tissueCol$col[match(umapAndMeta[tissueAndCol,"Tissue"], tissueCol$PlotClass)], alpha.f = 0.6)
tissue2AndCol <- umapAndMeta[,"Tissue2"] %in% tissueCol$PlotClass
umapAndMeta$col[tissue2AndCol] <-  adjustcolor(tissueCol$col[match(umapAndMeta[tissue2AndCol,"Tissue2"], tissueCol$PlotClass)], alpha.f = 0.6)


plotOrder <- order((umapAndMeta$col != defaultCol) + 1)

png(width = 800, height = 800)
plot(umapAndMeta[plotOrder,"UMAP1"], umapAndMeta[plotOrder,"UMAP2"], col = umapAndMeta$col[plotOrder], cex = 0.3, pch = 16)
dev.off()
