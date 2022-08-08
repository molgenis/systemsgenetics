#srun --cpus-per-task=1 --mem=200gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)



remoter::client("localhost", port = 55508, password = "laberkak")

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")


#save.image(file='celllineCancer.RData')
load('celllineCancer.RData')

library(pheatmap)

library(readr)


colTypes <- cols(
  .default = col_double(),
  `02/08/2022` = col_character()
)



table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/PCA/pc-scores.txt.gzip", delim = "\t", quote = "", col_types = colTypes)
pcs <- as.matrix(table_tmp[,-1])
rownames(pcs) <- table_tmp[,1][[1]]
rm(table_tmp)

str(pcs)

rpng()
heatmap3::heatmap3(cor(pcs[,1:10]), scale ="none", balanceColor=T, Rowv = NA, Colv = NA)
dev.off()

table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Metadata/metadataSRA.txt", delim = "\t", quote = "")
customMeta <- as.data.frame(table_tmp[,-1])
rownames(customMeta) <- table_tmp[,3][[1]]
rm(table_tmp)

table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Metadata/metadata_sra1.txt", delim = "\t", quote = "", guess_max = 20000)
sraMeta1 <- as.data.frame(table_tmp[,-1])
rm(table_tmp)

table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Metadata/metadata_sra2.txt", delim = "\t", quote = "", guess_max = 20000)
sraMeta2 <- as.data.frame(table_tmp[,-1])
rm(table_tmp)

sraSharedCol <- intersect(colnames(sraMeta2), colnames(sraMeta1))

  sraMeta <- rbind(sraMeta1[,sraSharedCol], sraMeta2[,sraSharedCol])
length(unique(sraMeta$external_id))
write.table("")

#For some reason some runs are duplicated in the meta data file. 
#Quick inspection showed that they have the same values
#Solution exclude duplicate row
sraUniqueIds <- unique(sraMeta$external_id)
str(sraUniqueIds)
sraMeta <- sraMeta[ match(sraUniqueIds, sraMeta$external_id), ]
rownames(sraMeta) <- sraMeta$external_id


sum(rownames(metaAndPcs) %in% sraMeta$external_id)

str(sraMeta)

sum(rownames(pcs) %in% rownames(customMeta))


metaAndPcs <- merge(customMeta, pcs, all.y = T, by = 0)
rownames(metaAndPcs) <- metaAndPcs[,"Row.names"]

sum(rownames(metaAndPcs) %in% rownames(sraMeta))

allMetaAndPcs <- merge(sraMeta, metaAndPcs, all.y = T, by = 0)
str(allMetaAndPcs)


table(allMetaAndPcs$Tissue)



allMetaAndPcs$colVecCellline <- adjustcolor("grey", alpha.f = 0.2)
allMetaAndPcs$colVecCellline[!is.na(allMetaAndPcs$CellLine)] <- adjustcolor("red", alpha.f = 0.5)


allMetaAndPcs$col <- adjustcolor("grey", alpha.f = 0.2)
allMetaAndPcs$col[!is.na(allMetaAndPcs$sra.library_layout) & allMetaAndPcs$sra.library_layout == "single"] <- adjustcolor("red", alpha.f = 0.5)
allMetaAndPcs$col[!is.na(allMetaAndPcs$sra.library_layout) & allMetaAndPcs$sra.library_layout == "paired"] <- adjustcolor("blue", alpha.f = 0.5)

rpng()
boxplot(allMetaAndPcs[,"PC_1"] ~ allMetaAndPcs$sra.library_layout)
t.test(allMetaAndPcs[,"PC_1"] ~ allMetaAndPcs$sra.library_layout)
dev.off()

table(allMetaAndPcs$sra.library_layout)
rpng()

pairs(allMetaAndPcs[,paste0("PC_",1:5)], upper.panel = NULL, col = allMetaAndPcs$col, cex = 0.2 )

dev.off()
#plot(pcs[,1], pcs[,2], col = metaAndPcs$colVecCellline )
dev.off()

allMetaAndPcs[allMetaAndPcs$PC_4 > 200,1:3]


pc1Sd <- sd(pcs[,1])
pc1Mean <- mean(pcs[,1])

rpng(width = 800 , height = 800)
plot(pcs[,1], pcs[,2], col = adjustcolor("grey", alpha.f = 0.8))
abline(v=pc1Mean - (4*pc1Sd), col = "red")
abline(v=pc1Mean + (4*pc1Sd), col = "red")
dev.off()

table(allMetaAndPcs$sra.run_published)
rpng()
plot(pcs[,1], allMetaAndPcs$sra.run_published, col = adjustcolor("grey", alpha.f = 0.8))
dev.off()

numColumns <- unlist(lapply(allMetaAndPcs, is.numeric))


allMetaAndPcsNumMatrix <- as.matrix(allMetaAndPcs[,numColumns])
str(allMetaAndPcsNumMatrix)

pc1Cor <- cor(allMetaAndPcsNumMatrix[,"PC_1"], allMetaAndPcsNumMatrix, use = "pairwise.complete.obs")
tail(sort(abs(pc1Cor[1,])), n = 30)
pc1Cor[1,"recount_qc.star.average_input_read_length"]


pc4Cor <- cor(allMetaAndPcsNumMatrix[,"PC_4"], allMetaAndPcsNumMatrix, use = "pairwise.complete.obs")
tail(sort(abs(pc4Cor[1,])), n = 10)



pcCor <- cor(allMetaAndPcsNumMatrix[,paste0("PC_", 1:50)], allMetaAndPcsNumMatrix, use = "pairwise.complete.obs")
str(pcCor)

rpng()
pdf("tmp.pdf", height = 12, width =20)
pheatmap::pheatmap(pcCor[,!colnames(pcCor) %in% c("sra.spot_length", paste0("PC_", 1:1000))], scale = "none", cluster_cols= FALSE, cluster_rows = FALSE)
dev.off()
tail(sort(abs(pcCor[1,])), n = 10)

rpng()
pairs(allMetaAndPcsNumMatrix[,c("recount_seq_qc.%a","recount_seq_qc.%c","recount_seq_qc.%t","recount_seq_qc.%g", "recount_seq_qc.%n", "recount_qc.intron_sum_%", "recount_qc.bc_auc.unique_%")])
dev.off()

rpng()
plot(allMetaAndPcs$`PC_1`, allMetaAndPcs$`recount_seq_qc.avg_len`, col = adjustcolor("dodgerblue2", alpha.f = 0.2))
cor.test(allMetaAndPcs$`PC_1`, allMetaAndPcs$`recount_seq_qc.avg_len`)
dev.off()


rpng()
plot(allMetaAndPcs$`recount_qc.star.mapping_speed,_million_of_reads_per_hour`, allMetaAndPcs$`recount_seq_qc.avg_len`, col = adjustcolor("dodgerblue2", alpha.f = 0.2))
cor.test(allMetaAndPcs$`recount_qc.star.mapping_speed,_million_of_reads_per_hour`, allMetaAndPcs$`recount_seq_qc.avg_len`)
dev.off()

summary(lm(allMetaAndPcs$`0` ~ allMetaAndPcs$`recount_qc.star.uniquely_mapped_reads_%` + allMetaAndPcs$`recount_seq_qc.avg_len` + allMetaAndPcs$`recount_qc.star.mapping_speed,_million_of_reads_per_hour`))


sum(allMetaAndPcs$`0` >= 0.0017)
sum(allMetaAndPcs$`0` < 0.0017)

rpng()
plot(allMetaAndPcs$`0`, allMetaAndPcs$`sra.sample_spots`*allMetaAndPcs$`recount_qc.star.uniquely_mapped_reads_%`, col = adjustcolor("dodgerblue2", alpha.f = 0.2))
cor.test(allMetaAndPcs$`0`, allMetaAndPcs$`sra.sample_spots`*allMetaAndPcs$`recount_qc.star.uniquely_mapped_reads_%`)
dev.off()



rpng()
plot(allMetaAndPcs$`recount_qc.star.average_input_read_length`, allMetaAndPcs$`recount_seq_qc.avg_len`, col = allMetaAndPcs$col)
dev.off()


rpng()
plot(allMetaAndPcs$`PC_1`, allMetaAndPcs$`recount_qc.star.%_mapped_reads_both`)
dev.off()



rpng()
plot(allMetaAndPcs$`PC_1`, allMetaAndPcs$`recount_qc.star.uniquely_mapped_reads_%`)
dev.off()


rpng()
plot(allMetaAndPcs$`recount_qc.star.%_mapped_reads_both`, allMetaAndPcs$`recount_qc.star.uniquely_mapped_reads_%`)
cor.test(allMetaAndPcs$`recount_qc.star.%_mapped_reads_both`, allMetaAndPcs$`recount_qc.star.uniquely_mapped_reads_%`)
dev.off()



sum(is.na(allMetaAndPcs$`recount_qc.star.%_of_reads_mapped_to_too_many_loci`))
sum(is.na(allMetaAndPcs$`recount_qc.star.%_reads_mapped_to_too_many_loci`))

  rpng()
hist( log10(allMetaAndPcs$`sra.sample_spots`), breaks = 100)
abline(v = log10(2e8), col = "red", lwd  =3)
abline(v = 6, col = "red", lwd  =3)
dev.off()

sum(allMetaAndPcs$`sra.sample_spots`>2e8 | allMetaAndPcs$`sra.sample_spots`<1e6, na.rm=T)

min(allMetaAndPcs$`sra.sample_spots`, na.rm=T)

rpng()
hist( allMetaAndPcs$`recount_qc.star.deletion_average_length`, breaks = 100) 
abline(v = 3, col = "red", lwd  =3)
dev.off()

dim(allMetaAndPcs)

sum(allMetaAndPcs$`recount_seq_qc.%n`>2, na.rm=T)


rpng()
hist( allMetaAndPcs$`recount_qc.star.uniquely_mapped_reads_%`, breaks = 100)
abline(v = 60, col = "red", lwd  =3)
dev.off()

rpng()
plot(allMetaAndPcs$PC_1, allMetaAndPcs$`recount_qc.intron_sum_%`)
cor.test(allMetaAndPcs$PC_1, allMetaAndPcs$`recount_qc.intron_sum_%`)
dev.off()


rpng()
plot(allMetaAndPcs$PC_1, allMetaAndPcs$`recount_qc.aligned_reads%.chrm`)
dev.off()


rpng()
plot(allMetaAndPcs$`recount_qc.aligned_reads%.chrx`, allMetaAndPcs$`recount_qc.aligned_reads%.chrm`)
dev.off()

rpng()
hist(allMetaAndPcs$`recount_qc.aligned_reads%.chrm`, breaks =100)
abline(v=20, col = "red", lwd  =3)
dev.off()

rpng()
hist(allMetaAndPcs$`recount_qc.aligned_reads%.chrx`, breaks =100)
abline(v=6, col = "red", lwd  =3)
dev.off()

rpng()
hist(allMetaAndPcs$`recount_qc.aligned_reads%.chry`, breaks =100)
abline(v=0.5, col = "red", lwd  =3)
dev.off()


cor.test(allMetaAndPcs$PC_1, allMetaAndPcs$`recount_qc.star.uniquely_mapped_reads_%`)

sum(allMetaAndPcs$`recount_qc.star.uniquely_mapped_reads_%` < 60, na.rm = TRUE)


rpng()
hist( allMetaAndPcs$`recount_qc.star.%_reads_mapped_to_multiple_loci_both`, breaks = 100)

dev.off()

tail(allMetaAndPcs$`sra.read_info`)
rpng()
x <- hist( allMetaAndPcs$`sra.spot_length`, breaks = c(seq(0,1000, by = 10),max( allMetaAndPcs$`sra.paired_nominal_length`, na.rm= T)), xlim = c(0,600))
dev.off()

rpng()
plot(allMetaAndPcs$`sra.paired_nominal_length`, allMetaAndPcs$`recount_qc.star.average_input_read_length`, xlim = c(0,600))
dev.off()

tail(sort(allMetaAndPcs$`sra.paired_nominal_length`), n = 1000)


meanStatistics <- apply(allMetaAndPcsNumMatrix, 2, mean, na.rm = T)
medianStatistics <- apply(allMetaAndPcsNumMatrix, 2, median, na.rm = T)

strangeSamples <- allMetaAndPcs[allMetaAndPcs$PC_4 > 200,numColumns]

write.table(rbind(strangeSamples, meanStatistics, medianStatistics, allMetaAndPcs[10000:10010,numColumns]), file = "tmp.txt", sep = "\t", quote = FALSE, col.names = NA)

str(meanStatistics)
  
  rpng()
  plot(log10(allMetaAndPcs[,"sra.sample_spots"]), allMetaAndPcs[,"PC_3"])
  dev.off()

cor.test(-log10(allMetaAndPcs[,"sra.sample_spots"]), allMetaAndPcs[,"PC_3"])



oldAnnotations <- read.delim("celllinesAndCancer/oldAnnotations/sampleAnnotations.txt")
tissueCol <- read.delim("celllinesAndCancer/oldAnnotations/tissueCol5.txt")

sum(oldAnnotations$Sample %in% allMetaAndPcs$SampleID)
sum(oldAnnotations$Sample %in% row.names(allMetaAndPcs))

allMetaAndPcs[10000:10010,c("SampleID","external_id")]


sum(allMetaAndPcs$external_id %in% oldAnnotations$Sample)
str(oldAnnotations)
pcsAndOldAnnotations <- merge( pcs, oldAnnotations, by.x = 0, by.y = "Sample", all.x = T, suffixes = c("", "_old"))

pcsAndOldAnnotations$col <- tissueCol$col[match(pcsAndOldAnnotations[,"PlotClass"], tissueCol$PlotClass, nomatch = nrow(tissueCol))]

plotOrder <- order(pcsAndOldAnnotations[,"PlotClass"] %in% tissueCol$PlotClass + 1)

rpng(width = 800, height = 800)
plot(pcsAndOldAnnotations[plotOrder,"PC_1"], pcsAndOldAnnotations[plotOrder,"PC_2"], col = pcsAndOldAnnotations$col[plotOrder], cex = 0.4)
dev.off()

 
rpng(width = 800, height = 800)
plot(pcsAndOldAnnotations[plotOrder,"PC_4"], pcsAndOldAnnotations[plotOrder,"PC_5"], col = pcsAndOldAnnotations$col[plotOrder], cex = 0.4)
dev.off()

rpng(width = 800, height = 800)
plot(pcsAndOldAnnotations[plotOrder,"PC_2"], pcsAndOldAnnotations[plotOrder,"PC_5"], col = pcsAndOldAnnotations$col[plotOrder], cex = 0.4)
dev.off()


rpng(width = 800, height = 800)
plot(pcsAndOldAnnotations$`recount_qc.star.%_mapped_reads_both`[plotOrder], pcsAndOldAnnotations[plotOrder,"PC_1"], col = pcsAndOldAnnotations$col[plotOrder], cex = 0.4)
dev.off()

cor.test(pcsAndOldAnnotations$`recount_qc.star.%_mapped_reads_both`[plotOrder], pcsAndOldAnnotations[plotOrder,"PC_1"])

rpng(width = 800, height = 800)
plot(pcsAndOldAnnotations[plotOrder,paste0("PC_",1:3)], upper.panel = NULL, col = pcsAndOldAnnotations$col[plotOrder], cex = 0.4)
dev.off()

pcsAndOldAnnotations$CellLine2 <- pcsAndOldAnnotations$CellLine_old != ""
table(pcsAndOldAnnotations$CellLine2)


pcsAndOldAnnotations$CellLineCol <-  adjustcolor("grey", alpha.f = 0.2)
pcsAndOldAnnotations$CellLineCol[pcsAndOldAnnotations$CellLine2] <-  adjustcolor("magenta", alpha.f = 0.5)

plotOrder <- order((adjustcolor("grey", alpha.f = 0.2) == pcsAndOldAnnotations$CellLineCol) * -1)


rpng(width = 800, height = 800)
pairs(pcsAndOldAnnotations[plotOrder,paste0("PC_",1:5)], col = pcsAndOldAnnotations$CellLineCol[plotOrder], upper.panel = NULL, cex = 0.3)
dev.off()

dim(pcsAndOldAnnotations)
str(as.factor(pcsAndOldAnnotations$CellLine2))
str(pcsAndOldAnnotations$PC_1)
x <- wilcox.test(pcsAndOldAnnotations$PC_1~as.factor(pcsAndOldAnnotations$CellLine2)) 

library(pROC)
test <- roc(as.factor(pcsAndOldAnnotations$CellLine2), pcsAndOldAnnotations$PC_1)
str(auc(test))


library(vioplot)


tmp <- pcsAndOldAnnotations[pcsAndOldAnnotations[,"PC_1"] > 200 | pcsAndOldAnnotations[,"PC_1"] < -200,]
str(tmp)
write.table(tmp, file = "tmp.txt", sep = "\t", quote = F)



columnsToCheckForSu4 <- c("sra.library_construction_protocol", 
                          "sra.study_abstract", 
                          "sra.experiment_title", 
                          "sra.design_description",
                          "sra.sample_description",
                          "sra.library_construction_protocol",
                          "sra.sample_attributes",
                          "sra.sample_title")


match4su <- Reduce("|", lapply(columnsToCheckForSu4, function(col, pattern, ...){grepl(pattern,pcsAndOldAnnotations[,col],...)}, pattern = "4su|thiouridine", ignore.case = T))
sum(match4su)

#sum(grepl("4su|thiouridine", pcsAndOldAnnotations$`sra.run_alias`, ignore.case = T))

pcsAndOldAnnotations$su4Col <-  adjustcolor("grey", alpha.f = 0.2)
pcsAndOldAnnotations$su4Col[match4su] <-  adjustcolor("magenta", alpha.f = 0.5)

plotOrder <- order((adjustcolor("grey", alpha.f = 0.2) == pcsAndOldAnnotations$su4Col) * -1)
str(pcsAndOldAnnotations)

write.table(pcsAndOldAnnotations[match4su,"external_id",drop =F], file = "4su_samples.txt", sep = "\t", quote  = FALSE, row.names = F)


rpng(width = 800, height = 800)
plot(pcsAndOldAnnotations[plotOrder,"PC_3"], pcsAndOldAnnotations[plotOrder,"PC_7"], col = pcsAndOldAnnotations$su4Col[plotOrder], cex = 0.4)
dev.off()

rpng(width = 800, height = 800)
pairs(pcsAndOldAnnotations[plotOrder,paste0("PC_",c(31,17,3,54,18))], col = pcsAndOldAnnotations$su4Col[plotOrder], cex = 0.4)
dev.off()

library(pROC)

x <- 
as.numeric(x)


su4TestP <- apply(pcsAndOldAnnotations[,paste0("PC_",1:100)], 2, function(x){t.test(x ~match4su)$p.value})
sort(su4TestP)

su4TestP <- apply(pcsAndOldAnnotations[,paste0("PC_",1:100)], 2, function(x){as.numeric(auc(response = match4su, predictor = x))})
sort(su4TestP)

