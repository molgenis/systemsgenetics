setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

library(eulerr)


set.seed(1)

humanSkeletal <- read.delim("heightVenn/humanSkeletal.txt")[,1]
mouseLethal <- read.delim("heightVenn/mouseLethal.txt")[,1]
mouseGrowth <- read.delim("heightVenn/mouseGrowth.txt")[,1]
ehlersDanlos <- read.delim("heightVenn/EhlersDanlos.txt")[,1]
oi <- read.delim("heightVenn/oi.txt")[,1]
heightKey <- read.delim("heightVenn/heightSignificant.txt")[,1]
heightTested <- read.delim("heightVenn/heightTested.txt")[,1]
collagenFibrilOrganization <- read.delim("heightVenn/collagenFibrilOrganization.txt")[,1]
extracellularMatrixOrganization <- read.delim("heightVenn/extracellularMatrixOrganization.txt")[,1]
embryonicDigitMorphogenesis <- read.delim("heightVenn/embryonicDigitMorphogenesis.txt")[,1]
osteoblastDifferentiation <- read.delim("heightVenn/osteoblastDifferentiation.txt")[,1]
endodermalCellDifferentiation <- read.delim("heightVenn/endodermalCellDifferentiation.txt")[,1]
cellMigration <- read.delim("heightVenn/cellMigration.txt")[,1]
cellMatrixAdhesion <- read.delim("heightVenn/cellMatrixAdhesion.txt")[,1]
basementMembraneOrganization <- read.delim("heightVenn/basementMembraneOrganization.txt")[,1]
positiveRegulationOfTranscriptionByRNAPolymeraseII <- read.delim("heightVenn/positiveRegulationOfTranscriptionByRNAPolymeraseII.txt")[,1]
proteinPeptidyl_prolylIsomerization <- read.delim("heightVenn/proteinPeptidyl_prolylIsomerization.txt")[,1]
skeletalSystemDevelopment <- read.delim("heightVenn/skeletalSystemDevelopment.txt")[,1]
abnormalityOfConnectiveTissue <- read.delim("heightVenn/AbnormalityOfConnectiveTissue.txt")[,1]
growthAbnormality <- read.delim("heightVenn/GrowthAbnormality.txt")[,1]
heightPascalSig <- read.delim("heightVenn/HeightPascalSig.txt")[,1]




library(readr)

#change X1 in case of specified header for row names. Use ...1 if first colun has no ID
colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\GO\\goa_human_2020_06_01.gaf_P_2020_06_01_matrix.txt", delim = "\t", quote = "", col_types = colTypes)
goP <- as.matrix(table_tmp[,-1])
rownames(goP) <- table_tmp[,1][[1]]
rm(table_tmp)


#change X1 in case of specified header for row names. Use ...1 if first colun has no ID
colTypes <- cols( .default = col_double(),  `-` = col_character())
table_tmp <- read_delim("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\2023_06_17\\HPO_2023_06_17.txt.gz", delim = "\t", quote = "", col_types = colTypes)
hpo <- as.matrix(table_tmp[,-1])
rownames(hpo) <- table_tmp[,1][[1]]
rm(table_tmp)




goP <- goP[rownames(goP) %in% heightTested,]

str(humanSkeletal)

str(mouseLethal)

str(mouseGrowth)
str( 	)
geneLists <- list(heightKey = heightKey, humanSkeletal= humanSkeletal, mouseGrowth = mouseGrowth, mouseLethal= mouseLethal, ehlersDanlos = ehlersDanlos, collagenFibrilOrganization= collagenFibrilOrganization)

str(geneLists)

plot(euler(geneLists, shape = "ellipse"))


geneLists <- list(heightDownstreamer = heightKey, ehlersDanlos = ehlersDanlos, oi = unique(oi[oi %in% heightTested]))

geneLists <- list(heightDownstreamer = heightKey, ehlersDanlos = ehlersDanlos[ehlersDanlos %in% heightKey],oi = unique(oi[oi %in% heightKey]), extracellularMatrixOrganization=  rownames(goP)[rownames(goP) %in% heightKey & goP[,"GO:0030198"]>0], ossification= rownames(goP)[rownames(goP) %in% heightKey & goP[,"GO:0001503"]>0])



geneLists <- list(ehlersDanlos = ehlersDanlos,oi = unique(oi), extracellularMatrixOrganization=  rownames(goP)[goP[,"GO:0030198"]>0])
geneLists <- list(ehlersDanlos = ehlersDanlos,oi = unique(oi), ossification= rownames(goP)[goP[,"GO:0001503"]>0])



rownames(hpo)[rownames(hpo) %in% heightKey & hpo[,"HP:0000002"]>0]

str(heightKey)


mgiGenes <- readRDS(file = "../GeneNetwork/MGI/2020_10_20/MGI_2020_10_20_PhenoGenoMP_matrix.rds")




geneLists <- list(heightDownstreamer = heightKey, AbnormalityOfBodyHeight  = rownames(hpo)[rownames(hpo) %in% heightKey & hpo[,"HP:0000002"]>0], MouseOrganogenesis  = rownames(mgiGenes)[rownames(mgiGenes) %in% heightKey & mgiGenes[,"MP:0011098"]>0], MouseSlowEmbryonicDevelopment  = rownames(mgiGenes)[rownames(mgiGenes) %in% heightKey & mgiGenes[,"MP:0003984"]>0])
plot(euler(geneLists, shape = "ellipse"))






geneLists <- list(heightDownstreamer = heightKey, ehlersDanlos = ehlersDanlos[ehlersDanlos %in% heightKey],oi = unique(oi[oi %in% heightKey]), extracellularMatrixOrganization=  rownames(goP)[rownames(goP) %in% heightKey & goP[,"GO:0030198"]>0],  AbnormalityOfBodyHeight  = rownames(hpo)[rownames(hpo) %in% heightKey & hpo[,"HP:0000002"]>0] )
plot(euler(geneLists, shape = "ellipse"))





#
pdf("heightVenn/vennOtherPathways.pdf")
geneLists2 <- list(heightKey = heightKey, humanSkeletal= humanSkeletal, extracellularMatrixOrganization = extracellularMatrixOrganization, embryonicDigitMorphogenesis= embryonicDigitMorphogenesis, osteoblastDifferentiation= osteoblastDifferentiation, cellMigration= cellMigration)
plot(euler(geneLists2), quantities = TRUE)
#plot(euler(geneLists2, shape = "ellipse"), quantities = TRUE)
dev.off()

pdf("heightVenn/vennMouse.pdf")
geneLists3 <- list(heightKey = heightKey, humanSkeletal= humanSkeletal, mouseGrowth = mouseGrowth, mouseLethal= mouseLethal)
plot(euler(geneLists3), quantities = TRUE)
#plot(euler(geneLists3, shape = "ellipse"), quantities = TRUE)
dev.off()


geneLists4 <- list(heightKey = heightKey, 
                   embryonicDigitMorphogenesis= embryonicDigitMorphogenesis, 
                   cellMigration= cellMigration, 
                   collagenFibrilOrganization = collagenFibrilOrganization, 
                   extracellularMatrixOrganization = extracellularMatrixOrganization, 
                   endodermalCellDifferentiation = endodermalCellDifferentiation, 
                   basementMembraneOrganization = basementMembraneOrganization, 
                   proteinPeptidyl_prolylIsomerization = proteinPeptidyl_prolylIsomerization, 
                   positiveRegulationOfTranscriptionByRNAPolymeraseII = positiveRegulationOfTranscriptionByRNAPolymeraseII)
pdf("heightVenn/vennPathways.pdf")
plot(euler(geneLists4, shape = "ellipse"), quantities = TRUE)
dev.off()

hpoMatrix <- readRDS("phenotype_to_genes_V1268_OMIMandORPHA.txt_matrix.rds")

#hpoMatrix <- as.matrix(read.delim("heightVenn/hpoSelection.txt", row.names = 1))
str(hpoMatrix)
hpoList <- sapply(c("HP.0001507","HP.0000924"), function(term){
  row.names(hpoMatrix)[hpoMatrix[,term]>0]
})
str(hpoList)
#

geneLists <- list(heightKey = heightKey,heightPascal= heightPascalSig, heightRare = unique(do.call(c, hpoList)))
#geneLists <- c(geneLists, hpoList)
str(geneLists)
plot(euler(geneLists, shape = "ellipse"), quantities = TRUE)


cat(intersect(intersect(geneLists[[1]],geneLists[[2]]), geneLists[[3]]), sep = "\n")










ibdKey <- read.delim("heightVenn/ibdKey.txt")[,1]
ibdGwas <- read.delim("heightVenn/ibdGwas.txt")[,1]

hpoList <- sapply(c("HP:0002715"), function(term){
  row.names(hpoMatrix)[hpoMatrix[,term]>0]
})
str(hpoList)
#

geneLists <- list(ibdKey = ibdKey,ibdGwas= ibdGwas, rare = unique(do.call(c, hpoList)))
#geneLists <- c(geneLists, hpoList)
str(geneLists)
plot(euler(geneLists, shape = "ellipse"), quantities = TRUE)


cat(intersect(intersect(geneLists[[1]],geneLists[[2]]), geneLists[[3]]), sep = "\n")





gwasGenes <- read.delim("keyCisTransRare/all_gws_genes.txt")[,1]
hpoGenes <- read.delim("keyCisTransRare/all_hpo_genes_gwas_matched.txt")[,1]
keyGenes <- read.delim("keyCisTransRare/all_key_genes.txt")[,1]

geneLists <- list(gwasGenes = gwasGenes, hpoGenes = hpoGenes, keyGenes = keyGenes)
str(geneLists)
pdf("keyCisTransRare/eulerKeyCisRare.pdf")
plot(euler(geneLists, shape = "ellipse"), quantities = TRUE)
dev.off()

sum(!(keyGenes %in% gwasGenes)) / length(keyGenes) 


1-(sum(!(keyGenes %in% hpoGenes)) / length(keyGenes) )
1-(sum(!(gwasGenes %in% hpoGenes)) / length(gwasGenes) )

