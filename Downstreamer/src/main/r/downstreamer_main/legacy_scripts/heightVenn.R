setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

library(eulerr)


set.seed(1)

humanSkeletal <- read.delim("heightVenn/humanSkeletal.txt")[,1]
mouseLethal <- read.delim("heightVenn/mouseLethal.txt")[,1]
mouseGrowth <- read.delim("heightVenn/mouseGrowth.txt")[,1]
ehlersDanlos <- read.delim("heightVenn/EhlersDanlos.txt")[,1]
heightKey <- read.delim("heightVenn/heightKey.txt")[,1]
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


str(humanSkeletal)

str(mouseLethal)

str(mouseGrowth)


geneLists <- list(heightKey = heightKey, humanSkeletal= humanSkeletal, mouseGrowth = mouseGrowth, mouseLethal= mouseLethal, ehlersDanlos = ehlersDanlos, collagenFibrilOrganization= collagenFibrilOrganization)

str(geneLists)

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

