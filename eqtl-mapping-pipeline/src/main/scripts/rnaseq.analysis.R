#
# rnaseq.analysis.R
#
# Copyright (c) 2013-2013 GBIC: Danny Arends, Harm-Jan Westra, Lude Franke and Ritsert C. Jansen
# last modified Jul, 2013
# first written Jun, 2013
#
setwd("~/Github/systemsgenetics/eqtl-mapping-pipeline/src/main/scripts")
source("helper.functions.R")

setwd("~/Github/Juha/")
# Load the RNA seq data
RNASeq      <- read.csv("expression_table.genes.exonic_v69.0.3.rawCounts_Named.txt", sep='\t', row.names=1)

# Affy neutrophils
Neutr       <- read.csv("GPL570_Neutrophil_Good.txt.txt", sep='\t', row.names=1)
translation <- read.affy.translation()

Neutr <- annotate.affy.neutrohils(Neutr, translation)
NeutrRNASeq <- match.annotated.affy.rnaseq(Neutr, RNASeq)

# LOG2 transform the RNA Seq data
RNASeqLog  <- log2(NeutrRNASeq[,7])
tokeep      <- which(is.finite(RNASeqLog)) #Keep only the finite ones
NeutrRNASeq <- NeutrRNASeq[keep, ]

RNASeqLog <- log2(NeutrRNASeq[,7])   # LOG 2 transform the data

# Illumina Celltype data
Illu <- read.illumina.celltypes()
Illu <- annotate.illumina.celltypes(Illu)
Illu <- add.illumina.probes.information(Illu)

inIllu <- which(Illu[,1] %in% rownames(NeutrRNASeq))
inAffy <- which(rownames(NeutrRNASeq) %in% Illu[,1])
Illu <- Illu[inIllu, ]
NeutrRNASeq <- NeutrRNASeq[inAffy, ]

# Add the illumina data to the affy & RNA seq data
RnaAffyIllu <- NULL
for(gene in rownames(NeutrRNASeq)){
  IlluMean <- mean(as.numeric(Illu[which(Illu[,1] == gene), 2]))
  AffyRowID <- which(rownames(NeutrRNASeq)==gene)
  RnaAffyIllu <- rbind(RnaAffyIllu, c(gene, IlluMean, NeutrRNASeq[AffyRowID,-8]))
}

# Write the full table
write.table(RnaAffyIllu,"NeutrophilIllumina_AllRNAseq_AllNeutroAffy_ByHUGO.txt", sep='\t', quote=FALSE)

AffyMean  <- apply(RnaAffyIllu[,9:ncol(Neutr)], 1, mean)
IlluMean  <- apply(RnaAffyIllu[,2], 1, mean)

CorAffyRNASeqMean <- round(cor(AffyMean, RNASeqLog, method="spearman"), d = 2)
plot(AffyMean, RNASeqLog, xlab = "Affy", ylab = "RNAseq", main=paste0("Mean Cor: ",CorAffyRNASeqMean), cex=0.7)

