#
# helper.functions.R
#
# Copyright (c) 2013-2013 GBIC: Danny Arends, Harm-Jan Westra, Lude Franke and Ritsert C. Jansen
# last modified Jul, 2013
# first written Jun, 2013
#

read.illumina.celltypes <- function(){
  CellTypeDATA <- read.table("E-TABM-633/HaemAtlasMKEBNormalizedIntensities.csv",sep='\t',row.names=1,header=TRUE)
  CellTypeDATA <- CellTypeDATA[-1,]
  tmp <- apply(CellTypeDATA, 2, as.numeric)
  rownames(tmp) <- rownames(CellTypeDATA)
  CellTypeDATA  <- tmp
  return(CellTypeDATA)
}

add.illumina.probes.information(CellTypeDATA){
  ProbeAnnotation <- read.csv("GPL6102-11574.txt",sep="\t",skip=27) # Probe annotation
  ProbeAnnotation <- ProbeAnnotation[which(ProbeAnnotation[,"Symbol"] != ""),]
  inAnnot <- which(rownames(CellTypeDATA) %in% ProbeAnnotation[,1])
  CellTypeDATA <- CellTypeDATA[inAnnot,]
  sortAnnot <- match(names(CellTypeDATA), ProbeAnnotation[,1]) # Align
  CellTypeDATA <- cbind(as.character(ProbeAnnotation[sortAnnot,"Symbol"]), CellTypeDATA)
  return(CellTypeDATA)
}

annotate.illumina.celltypes <- function(CellTypeDATA){
  CellTypeAnnotation <- read.csv("dataDescr/E-TABM-633.txt",sep="\t", header=TRUE)
  cellTypeNames <- as.character(CellTypeAnnotation[, "Hybridization.Name"])
  cellTypeTypes <- as.character(CellTypeAnnotation[, "Extract.Name"])
  cellTypes     <- unlist(lapply(strsplit(cellTypeTypes,"-"),"[",3))

  cellTypeAnnot <- cbind(cellTypeNames, cellTypes) # FInally the annotation we want Hyb ref -> Cell-type
  CellTypeDATA  <- CellTypeDATA[, match(cellTypeAnnot[,1], gsub("X","",colnames(CellTypeDATA)))] # ARRANGE

  colnames(CellTypeDATA) <- as.character(cellTypeAnnot[,2])
  return(CellTypeDATA)
}

match.annotated.affy.rnaseq <- function(Neutr, RNASeq){
  inSeq   <- which(as.character(Neutr[,1]) %in% rownames(RNASeq))
  Neutr   <- Neutr[inSeq,]
  sortSeq <- match(as.character(Neutr[,1]), rownames(RNASeq)) # Align
  Neutr <- cbind(RNASeq[sortSeq,], Neutr)
  return(Neutr)
}

annotate.affy.neutrohils <- function(Neutr, translation){
  inTrans <- which(rownames(Neutr) %in% translation[,1])
  Neutr <- Neutr[inTrans,]
  sortTrans <- match(rownames(Neutr), translation[,1]) # Align
  Neutr <- cbind(translation[sortTrans,9], Neutr)
  return(Neutr)
}

read.affy.translation <- function(){
  translation <- read.csv("GPL570ProbeENSGInfo+HGNC.txt",sep='\t',row.names=1)
  translation <- translation[which(translation[,9] != "-"),]
  return(translation)
}

annotate.RNASeq <- function(){
  RNASeq           <- read.csv("expression_table.genes.exonic_v69.0.3.rawCounts_Named.txt", sep='\t', row.names=1)
  sampleNames      <- read.csv("SampleDetails.txt", sep='\t', row.names=1)
  colnames(RNASeq) <- sampleNames[,1]
  write.table(RNASeq,file="expression_table.genes.exonic_v69.0.3.rawCounts_Named.txt",sep='\t',quote=FALSE)
}

remove.bad.Neutrophils <- function(){
  Neutr      <- read.csv("GPL570_Neutrophil.txt", sep='\t', row.names=1)
  badneutros <- c("GSM141250", "GSM141251", "GSM141252", "GSM141253", "GSM141254", "GSM141255", "GSM141256",
                "GSM141257", "GSM549581", "GSM549582", "GSM549583", "GSM549584")
  Neutr  <- Neutr[,-which(colnames(Neutr) %in% badneutros)]  # Remove remaining bad neutrophil samples
  write.table(Neutr,file="GPL570_Neutrophil_Good.txt",sep='\t',quote=FALSE)
}

