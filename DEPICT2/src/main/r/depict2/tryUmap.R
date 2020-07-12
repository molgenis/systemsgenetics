
projectDir <- getwd()
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

source(paste0(projectDir, "/depict2_functions.r"))

platelets <- read.depict2("excelsForEnrichment/30080_raw.gwas.imputed_v3.both_sexes_enrichtments_exHla_v69.xlsx")
str(platelets$Coregulation)


geneP <- read.de



library(readr)
table_tmp <- read_delim("reference_datasets/human_b37/pathway_databases/gene_network/hpo_prepared_predictions.txt", delim = "\t", quote = "")
hpoPredictions <- as.matrix(table_tmp[,-1])
rownames(hpoPredictions) <- table_tmp[,1][[1]]
rm(table_tmp)


library(RColorBrewer)
library(umap)

str(hpoPredictions)

custom.settings = umap.defaults
custom.settings$n_neighbors = 100
custom.settings
hpoPredictionsUmap = umap(t(hpoPredictions[row.names(hpoPredictions) %in% genes$Ensembl.Gene.ID,]), config = custom.settings)


hpoClass <- read.delim("hpoClass", stringsAsFactors = T)
hpoClass <- hpoClass[match(row.names(hpoPredictionsUmap$layout), hpoClass$term),]
hpoClass$class <- factor(hpoClass$class, levels=c(levels(hpoClass$class), "Other"))
hpoClass$class[hpoClass$class == ""] <- "Other"
hpoClass$class[grep("Neoplasm", hpoClass$class)] <- "Neoplasm"
hpoClass$class[hpoClass$class %in% names(table(hpoClass$class)[table(hpoClass$class) <= 10])] <- "other"


hpoClass$class <- droplevels(hpoClass$class)
str(hpoClass)

View(table(hpoClass$class))

hpoClass$col <- "black"

hpoClass$col[grep("", hpoClass$class)] <- "red3"





str(hpoPredictionsUmap)

plot(hpoPredictionsUmap$layout, col = hpoClass$class)

genes <- read.delim("ensgR75_protein_coding.txt", stringsAsFactors = F)

#This adds a column of color values
# based on the y values
x <- brewer.pal(9, "GnBu")[as.numeric(cut(platelets$HPO$Enrichment.Z.score[match(rownames(hpoPredictionsUmap$layout), platelets$HPO$Gene.set)],breaks = 9))]

plot(hpoPredictionsUmap$layout, col = x)

row.names(hpoPredictionsUmap$layout)[hpoPredictionsUmap$layout[,2]> 20]

plot(hpoPredictionsUmap$layout[], bg = adjustcolor("dodgerblue2", alpha.f = 0.1), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.1))

