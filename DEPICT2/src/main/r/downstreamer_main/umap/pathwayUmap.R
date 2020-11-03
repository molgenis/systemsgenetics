#remoter::server(verbose = T, port = 55557, password = "laberkak", sync = T)
remoter::client("localhost", port = 55557, password = "laberkak")

setwd("/groups/umcg-wijmenga/tmp04/projects/depict2/")

library(readr)
table_tmp <- read_delim("depict2_bundle/reference_datasets/human_b37/pathway_databases/gene_network/hpo_prepared_predictions.txt", delim = "\t", quote = "")
hpoPredictions <- as.matrix(table_tmp[,-1])
rownames(hpoPredictions) <- table_tmp[,1][[1]]
rm(table_tmp)

str(hpoPredictions)

hpoPredictionsPca <- prcomp(hpoPredictions, scale. = TRUE)
hpoPredictionsPca$sdevPercentage <- (hpoPredictionsPca$sdev * 100) / sum(hpoPredictionsPca$sdev)
compsToUse <- min(which(cumsum(hpoPredictionsPca$sdevPercentage) >= 50))
hpoUmap <- umap(hpoPredictionsPca$rotation[,1:compsToUse], n_epochs = 5000, init = hpoPredictionsPca$rotation[,1:2])

plot(hpoUmap)
dev.off()



hpoClass <- read.delim("Umap/hpoClass", stringsAsFactors = T)
hpoClass <- hpoClass[match(row.names(hpoPredictionsPca$rotation), hpoClass$term),]
hpoClass$class <- factor(hpoClass$class, levels=c(levels(hpoClass$class), "Other"))
hpoClass$class[hpoClass$class == ""] <- "Other"
hpoClass$class[grep("Neoplasm", hpoClass$class)] <- "Neoplasm"
hpoClass$class[grep("Abnormality of the skeletal system", hpoClass$class)] <- "Abnormality of the skeletal system"
hpoClass$class[hpoClass$class %in% names(table(hpoClass$class)[table(hpoClass$class) <= 10])] <- "other"


table(droplevels(hpoClass$class[grep("Abnormality of metabolism/homeostasis", hpoClass$class)]))

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


plot(hpoUmap, col = hpoClass$col)
dev.off()

cat(levels(hpoClass$class), sep = "\n")
