#
# helper.functions.R
#
# Copyright (c) 2013-2013 GBIC: Danny Arends, Harm-Jan Westra, Lude Franke and Ritsert C. Jansen
# last modified Jul, 2013
# first written Jun, 2013
#

read.harmjan.zscores <- function(){
  if(file.exists("metaRes.Rdata")){
    load("metaRes.Rdata")
  }else{
    metaRes <- read.csv("MetaAnalysisZScoreMatrix-Ensembl.txt",sep='\t',row.names=1)
    save(metaRes, file="metaRes.Rdata")
  }
  return(metaRes)
}

read.illumina.dir <- function(folder, annotation){
  data <- NULL
  for(x in paste0(folder, "/", annotation[,1], ".txt")){
    cat(x,"\n")
    m <- read.csv(x, header=TRUE, row.names=1,sep="\t")
    data <- cbind(data, apply(m,2,as.numeric)[,2])
  }
  rownames(data) <- rownames(m)
  data <- cbind(m[,1], data)
  colnames(data) <- c("Array_Address_Id", annotation[,1])
  invisible(data)
}

read.illumina.wb <- function(){
  GSE19790DATA    <- read.csv("GSE19790/GSE19790_non-normalized.txt",sep="\t",row.names=1)
  GSE19790DATA    <- GSE19790DATA[, grep("AVG_Signal", colnames(GSE19790DATA))]
  colnames(GSE19790DATA) <- gsub(".AVG_Signal","", colnames(GSE19790DATA))
  colnames(GSE19790DATA) <- gsub("X","Ind", colnames(GSE19790DATA))

  GSE24757DATA <- read.csv("GSE24757/GSE24757_non-normalized_data.txt",sep="\t",row.names=1, skip=4,header=TRUE)
  sortGSE24757 <- match(rownames(GSE19790DATA), rownames(GSE24757DATA))
  ALL <- cbind(GSE24757DATA[sortGSE24757,], GSE19790DATA)

  GSE13255A <- read.csv("dataDescr/GSE13255.txt", sep="\t", header=FALSE, colClasses=("character")) # Annotation
  GSE13255DATA <- read.illumina.dir("GSE13255", GSE13255A)

  GSE13255DATA <- GSE13255DATA[which(rownames(GSE13255DATA) %in% rownames(ALL)),]
  ALL <- ALL[which(rownames(ALL) %in% rownames(GSE13255DATA)),]

  sortGSE13255 <- match(rownames(ALL), rownames(GSE13255DATA))

  ALL <- cbind(GSE13255DATA[sortGSE13255,-1], ALL)
  return(ALL)
}

read.illumina.celltypes <- function(){
  CellTypeDATA <- read.table("HaemAtlasMKEBNormalizedIntensities.csv",sep='\t',row.names=1,header=TRUE)
  CellTypeDATA <- CellTypeDATA[-1,]
  tmp <- apply(CellTypeDATA, 2, as.numeric)
  rownames(tmp) <- rownames(CellTypeDATA)
  CellTypeDATA  <- tmp
  return(CellTypeDATA)
}

read.illumina.probes.information <- function(){
  ProbeAnnotation <- read.csv("GPL6102-11574.txt",sep="\t",skip=27) # Probe annotation
  ProbeAnnotation <- ProbeAnnotation[which(ProbeAnnotation[,"Symbol"] != ""),]
  return(ProbeAnnotation)
}

add.illumina.probes.information <- function(CellTypeDATA){
  ProbeAnnotation <- read.illumina.probes.information()
  inAnnot <- which(rownames(CellTypeDATA) %in% ProbeAnnotation[,1])
  CellTypeDATA <- CellTypeDATA[inAnnot,]
  sortAnnot <- match(rownames(CellTypeDATA), ProbeAnnotation[,1]) # Align
  CellTypeDATA <- cbind(as.character(ProbeAnnotation[sortAnnot,"Symbol"]), CellTypeDATA)
  colnames(CellTypeDATA)[1] <- "Symbol"
  return(CellTypeDATA)
}

annotate.illumina.celltypes <- function(CellTypeDATA){
  CellTypeAnnotation <- read.csv("E-TABM-633.txt",sep="\t", header=TRUE)
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

add.metares.to.tmatrix <- function(metaRes, tscores_Annotated){
  # Bind the tscores_Annotated columns Neutrophil, Bcell and Tcell to the MetaAnalysis data
  metaRes <- metaRes[which(rownames(metaRes) %in% tscores_Annotated[,1]),]     # match MetaRes
  tscores_Annotated <- tscores_Annotated[which(tscores_Annotated[,1] %in% rownames(metaRes)),] # match ResVector
  sortRes <- match(rownames(metaRes), tscores_Annotated[,1]) # Align
  metaRes <- cbind(tscores_Annotated[sortRes, 3:7], metaRes)
  return(metaRes)
}

annotate.affy.by.rownames <- function(Neutr, translation){
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

get.affy.mean <- function(RnaAffyIllu, selection = 1:nrow(RnaAffyIllu)){
  res <- apply(RnaAffyIllu[selection,10:ncol(Neutr)], 1, function(x){mean(as.numeric(x))})
  names(res) <- RnaAffyIllu[selection,"HUGO"]
  res
}

get.illu.mean <- function(RnaAffyIllu, selection = 1:nrow(RnaAffyIllu)){
  res <- as.numeric(RnaAffyIllu[selection,2])
  names(res) <- RnaAffyIllu[selection,"HUGO"]
  res
}

get.rnaseq.mean <- function(RnaAffyIllu, selection = 1:nrow(RnaAffyIllu)){
  res <- as.numeric(RnaAffyIllu[selection,"granulocytes"])
  names(res) <- RnaAffyIllu[selection,"HUGO"]
  res
}

plot.single <- function(M1, M2, N1 = "Affy", N2 = "RNAseq", col = rep(1, length(M1))){
  corrr <- round(cor(M1, as.numeric(M2), method="spearman"), d = 2)
  plot(M1, M2, xlab = N1, ylab = N2, main=paste0("Mean Cor: ",corrr), cex=0.7, col=col)
}

plot.AffyIllu <- function(RnaAffyIllu, selection = 1:nrow(RnaAffyIllu)){
  AffyMean  <- get.affy.mean(RnaAffyIllu, selection)
  IlluMean  <- get.illu.mean(RnaAffyIllu, selection)
  RNASeqLog <- get.rnaseq.mean(RnaAffyIllu, selection)
  #op <- par(mfrow=c(1,3))
  CorAffyRNASeqMean <- round(cor(AffyMean, as.numeric(RNASeqLog), method="spearman"), d = 2)
  plot(AffyMean, RNASeqLog, xlab = "Affy", ylab = "RNAseq", main=paste0("Mean Cor: ",CorAffyRNASeqMean), cex=0.7)

  CorIlluRNASeqMean <- round(cor(IlluMean, as.numeric(RNASeqLog), method="spearman"), d = 2)
  plot(IlluMean, RNASeqLog, xlab = "Illu", ylab = "RNAseq", main=paste0("Mean Cor: ",CorIlluRNASeqMean), cex=0.7)

  CorIlluAffyMean <- round(cor(IlluMean, AffyMean, method="spearman"), d = 2)
  plot(IlluMean, AffyMean, xlab = "Illu", ylab = "Affy", main=paste0("Mean Cor: ",CorIlluAffyMean), cex=0.7)
}

annotate.RNASeq <- function(){
  RNASeq           <- read.csv("expression_table.genes.exonic_v69.0.3.rawCounts_Named.txt", sep='\t', row.names=1)
  sampleNames      <- read.csv("SampleDetails.txt", sep='\t', row.names=1)
  colnames(RNASeq) <- sampleNames[,1]
  write.table(RNASeq,file="expression_table.genes.exonic_v69.0.3.rawCounts_Named.txt",sep='\t',quote=FALSE)
}

is.outlier <- function(d1, d2, cutoff = 0.3){
  up <- as.numeric((d1/max(d1)) - (d2/max(d2)) > cutoff)
  down <- as.numeric((d1/max(d1)) - (d2/max(d2)) < -cutoff)
  res <- up+(2*down) +1
  names(res) <- names(d1)
  res
}

remove.bad.Neutrophils <- function(){
  Neutr      <- read.csv("GPL570_Neutrophil.txt", sep='\t', row.names=1)
  badneutros <- c("GSM141250", "GSM141251", "GSM141252", "GSM141253", "GSM141254", "GSM141255", "GSM141256",
                "GSM141257", "GSM549581", "GSM549582", "GSM549583", "GSM549584")
  Neutr  <- Neutr[,-which(colnames(Neutr) %in% badneutros)]  # Remove remaining bad neutrophil samples
  write.table(Neutr,file="GPL570_Neutrophil_Good.txt",sep='\t',quote=FALSE)
}

