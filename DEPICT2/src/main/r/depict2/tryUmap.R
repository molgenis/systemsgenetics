
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
custom.settings$n_neighbors = 200
custom.settings$min_dist = 0.01
custom.settings$n_epochs = 1000
custom.settings$knn_repeats = 3
hpoPredictionsUmap = umap(t(hpoPredictions[row.names(hpoPredictions) %in% genes$Ensembl.Gene.ID,]), config = custom.settings)
plot(hpoPredictionsUmap$layout, col = hpoClass$col)

hpoClass <- read.delim("hpoClass", stringsAsFactors = T)
hpoClass <- hpoClass[match(row.names(hpoPredictionsUmap$layout), hpoClass$term),]
hpoClass$class <- factor(hpoClass$class, levels=c(levels(hpoClass$class), "Other"))
hpoClass$class[hpoClass$class == ""] <- "Other"
hpoClass$class[grep("Neoplasm", hpoClass$class)] <- "Neoplasm"
hpoClass$class[grep("Abnormality of the skeletal system", hpoClass$class)] <- "Abnormality of the skeletal system"
hpoClass$class[hpoClass$class %in% names(table(hpoClass$class)[table(hpoClass$class) <= 10])] <- "other"


hpoClass$class <- droplevels(hpoClass$class)
str(hpoClass)

View(table(hpoClass$class))

hpoClass$col <- "black"

hpoClass$col[c(grep("immune system", hpoClass$class), grep("blood", hpoClass$class))] <- "red3"
hpoClass$col[c(grep("immune system", hpoClass$class), grep("blood", hpoClass$class))] <- "red3"
hpoClass$col[ hpoClass$class == "Abnormality of the eye" ] <- "dodgerblue2"
hpoClass$col[ hpoClass$class == "Abnormality of head or neck" ] <- "yellow2"
hpoClass$col[ hpoClass$class == "Abnormality of the cardiovascular system" ] <- "orange1"
hpoClass$col[ hpoClass$class == "Abnormality of the nervous system" ] <- "purple"
hpoClass$col[ hpoClass$class == "Abnormality of the skeletal system" ] <- "ivory3"
hpoClass$col[ hpoClass$class == "Abnormality of the genitourinary system" ] <- "pink1"
hpoClass$col[ hpoClass$class == "Abnormality of the digestive system" ] <- "coral4"
hpoClass$col[ hpoClass$class == "Abnormality of head or neck" ] <- "ivory3"
hpoClass$col[ hpoClass$class == "Neoplasm" ] <- "green"
hpoClass$col[ hpoClass$class == "Abnormality of the musculature" ] <- "deeppink"
hpoClass$col[ hpoClass$class == "Abnormality of metabolism/homeostasis" ] <- "skyblue"


plot(hpoPredictionsUmap$layout, col = hpoClass$col)

genes <- read.delim("ensgR75_protein_coding.txt", stringsAsFactors = F)

#This adds a column of color values
# based on the y values
x <- brewer.pal(9, "GnBu")[as.numeric(cut(platelets$HPO$Enrichment.Z.score[match(rownames(hpoPredictionsUmap$layout), platelets$HPO$Gene.set)],breaks = 9))]

plot(hpoPredictionsUmap$layout, col = x)

row.names(hpoPredictionsUmap$layout)[hpoPredictionsUmap$layout[,1] < -2.5 & hpoPredictionsUmap$layout[,2] < -2.5]


plot(hpoPredictionsUmap$layout[], bg = adjustcolor("dodgerblue2", alpha.f = 0.1), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.1))

