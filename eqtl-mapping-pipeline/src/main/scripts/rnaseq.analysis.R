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
fn <- "NeutrophilIllumina_AllRNAseq_AllNeutroAffy_ByHUGO.txt"
if(!file.exists(fn)){
  # Add the illumina data to the affy & RNA seq data
  cat("HUGO","Illumina(G)", unlist(colnames(NeutrRNASeq)[-8]),"\n",sep='\t', file = fn)
  cnt <- 0; n <- length(rownames(NeutrRNASeq))
  for(gene in rownames(NeutrRNASeq)){
    IlluMean <- mean(as.numeric(Illu[which(Illu[,1] == gene), -1]))
    AffyRowID <- which(rownames(NeutrRNASeq)==gene)
    cat(gene, IlluMean, unlist(NeutrRNASeq[AffyRowID,-8]),"\n", file=fn, append=TRUE, sep='\t')
    if(cnt %% 500 == 0)cat("Done",cnt,"out of",n,"\n")
    cnt <- cnt+1
  }
}
RnaAffyIllu <- read.csv("NeutrophilIllumina_AllRNAseq_AllNeutroAffy_ByHUGO.txt", sep = '\t')

plot.AffyIllu(RnaAffyIllu)

geneAnnotations <- read.csv("EnsemblGeneAnnotation.txt",sep='\t')
RnaAffyIllu <- RnaAffyIllu[which(RnaAffyIllu[,1] %in% geneAnnotations[,2]),]
inAnnot <- match(RnaAffyIllu[,1],geneAnnotations[,2])
RnaAffyIlluAnnotated <- cbind(geneAnnotations[inAnnot,],RnaAffyIllu)

XYgenes <- which(!(RnaAffyIlluAnnotated[,"Chromosome.Name"]== "X"| RnaAffyIlluAnnotated[,"Chromosome.Name"]== "Y"))

plot.AffyIllu(RnaAffyIllu, XYgenes)

highRNAseq <- which(RnaAffyIllu[,"granulocytes"] > 3)
plot.AffyIllu(RnaAffyIllu, highRNAseq)

AffyMean  <- get.affy.mean(RnaAffyIllu)
IlluMean  <- get.illu.mean(RnaAffyIllu)
RNASeqLog <- get.rnaseq.mean(RnaAffyIllu)

aa <- is.outlier(AffyMean, RNASeqLog, 0.2)
plot.single(AffyMean, RNASeqLog, col=aa)

aa <- is.outlier(AffyMean+2, IlluMean, 0.2)
plot.single(IlluMean, AffyMean, N1="Illu", N2= "Affy", col=aa)

cor(AffyMean, RNASeqLog,use='spearman')
is.outlier(IlluMean, RNASeqLog)

