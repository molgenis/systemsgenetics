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
Neutr       <- read.csv("GPL570_Neutrophil_Good.txt", sep='\t', row.names=1)
translation <- read.affy.translation()

Neutr       <- annotate.affy.by.rownames(Neutr, translation)
NeutrRNASeq <- match.annotated.affy.rnaseq(Neutr, RNASeq)

# LOG2 transform the RNA Seq data
RNASeqLog  <- log2(NeutrRNASeq[,"granulocytes"])
tokeep      <- which(is.finite(RNASeqLog)) #Keep only the finite ones
NeutrRNASeq <- NeutrRNASeq[tokeep, ]

NeutrRNASeq[,"granulocytes"] <- log2(NeutrRNASeq[,"granulocytes"])   # LOG 2 transform the data

# Load the Illumina celltype data
setwd("/home/danny/Github/LudeNew")
Illu <- read.illumina.celltypes()
Illu <- annotate.illumina.celltypes(Illu)
Illu <- Illu[, which(colnames(Illu) == "Granulocyte")]    # Take only the Granulocytes
Illu <- add.illumina.probes.information(Illu)

inIllu <- which(Illu[,1] %in% rownames(NeutrRNASeq))
inAffy <- which(rownames(NeutrRNASeq) %in% Illu[,1])
Illu <- Illu[inIllu, ]
NeutrRNASeq <- NeutrRNASeq[inAffy, ]

setwd("/home/danny/Github/Juha")
RnaAffyIllu <- NULL
if(!file.exists("NeutrophilIllumina_AllRNAseq_AllNeutroAffy_ByHUGO.txt")){
  # Add the illumina data to the affy & RNA seq data
  cnt <- 0; n <- length(rownames(NeutrRNASeq))
  for(gene in rownames(NeutrRNASeq)){
    IlluMean <- mean(as.numeric(Illu[which(Illu[,1] == gene), 2]))
    AffyRowID <- which(rownames(NeutrRNASeq)==gene)
    RnaAffyIllu <- rbind(RnaAffyIllu, c(gene, IlluMean, NeutrRNASeq[AffyRowID,-8]))
    if(cnt %% 500 == 0)cat("Done",cnt,"out of",n,"\n")
    cnt <- cnt+1
  }
  colnames(RnaAffyIllu)[1:2] <- c("HUGO","Illumina(G)")
  # Write the full table
  write.table(RnaAffyIllu,"NeutrophilIllumina_AllRNAseq_AllNeutroAffy_ByHUGO.txt", sep = '\t', quote = FALSE)
}else{
  RnaAffyIllu <- read.csv("NeutrophilIllumina_AllRNAseq_AllNeutroAffy_ByHUGO.txt", sep = '\t', row.names = 1)
  colnames(RnaAffyIllu)[1:2] <- c("HUGO","Illumina(G)")
}

plot.AffyIllu(RnaAffyIllu)

geneAnnotations <- read.csv("EnsemblGeneAnnotation.txt",sep='\t')
RnaAffyIllu <- RnaAffyIllu[which(RnaAffyIllu[,1] %in% geneAnnotations[,2]),]
inAnnot <- match(RnaAffyIllu[,1],geneAnnotations[,2])
RnaAffyIlluAnnotated <- cbind(geneAnnotations[inAnnot,],RnaAffyIllu)

XYgenes <- which(!(RnaAffyIlluAnnotated[,"Chromosome.Name"]== "X"| RnaAffyIlluAnnotated[,"Chromosome.Name"]== "Y"))

plot.AffyIllu(RnaAffyIllu, XYgenes)

AffyMean  <- get.affy.mean(RnaAffyIllu, XYgenes)
IlluMean  <- get.illu.mean(RnaAffyIllu, XYgenes)
RNASeqLog <- get.rnaseq.mean(RnaAffyIllu, XYgenes)

aa <- is.outlier(AffyMean, RNASeqLog)
cor(AffyMean, RNASeqLog,use='spearman')
is.outlier(IlluMean, RNASeqLog)

