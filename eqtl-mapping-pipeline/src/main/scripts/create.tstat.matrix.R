#
# create.tstat.matrix.R
#
# Copyright (c) 2013-2013 GBIC: Danny Arends, Harm-Jan Westra, Lude Franke and Ritsert C. Jansen
# last modified Jul, 2013
# first written Jun, 2013
#
setwd("~/Github/Juha/")

# Load the RNA seq data
RNASeq           <- read.csv("expression_table.genes.exonic_v69.0.3.rawCounts.txt", sep='\t', row.names=1)
sampleNames      <- read.csv("SampleDetails.txt", sep='\t', row.names=1)
colnames(RNASeq) <- sampleNames[,1]

# Load the Z scores from Meta Analysis
#metaRes <- read.csv("MetaAnalysisZScoreMatrix-Ensembl.txt",sep='\t',row.names=1)
#save(metaRes, file="metaRes.Rdata")
load("metaRes.Rdata")

# Load the cell type expression measured on GPL570
wholeblood <- read.csv("GPL570_WholeBlood.txt", sep='\t', row.names=1)
Neutr      <- read.csv("GPL570_Neutrophil.txt", sep='\t', row.names=1)
Bcell      <- read.csv("GPL570_Bcell.txt", sep='\t', row.names=1)
Tcell      <- read.csv("GPL570_Tcell.txt", sep='\t', row.names=1)
NKcell     <- read.csv("GPL570_NKcell.txt", sep='\t', row.names=1)
RBC        <- read.csv("GPL570_RBC.txt", sep='\t', row.names=1)

# Load the probe annotation 'GPL570', remove probes not targetting a genes
translation <- read.csv("GPL570ProbeENSGInfo+HGNC.txt",sep='\t',row.names=1)
translation <- translation[which(translation[,9] != "-"),]

# Load the probe annotation for Illumina, remove probes not targetting a genes
ivectorAnn <- read.table("IlluminaProbeTranslationTable.txt",sep='\t',header=TRUE)
ivectorAnn <- ivectorAnn[which(ivectorAnn[,5] != "-"),]

# Load the vector from Lude and the interaction vector (from Harm-Jan)
LudeVec <- read.csv("NeutrophilVectorLude.txt",sep="\t")
IntrVec <- read.csv("2013-06-21-EGCUT-Vector-rs12057769-2000128.txt",sep='\t',row.names=NULL)

# Preprocessing remove some 'Bad' samples from neutrophils. Use correlation for the other cell-types
badneutros <- c("GSM141250", "GSM141251", "GSM141252", "GSM141253", "GSM141254", "GSM141255", "GSM141256",
                "GSM141257", "GSM549581", "GSM549582", "GSM549583", "GSM549584")

Neutr  <- Neutr[,-which(colnames(Neutr) %in% badneutros)]  # Remove remaining bad neutrophil samples
Bcell  <- Bcell[,-which(apply(cor(Bcell), 2, median) < 0.9)] # Use correlation for Bcell
Tcell  <- Tcell[,-which(apply(cor(Tcell), 2, median) < 0.9)] # Use correlation for Tcell

#Only take the interaction vector elements that match the Illumina annotation
IntrVec <- IntrVec[which(as.character(IntrVec[,1]) %in% as.character(ivectorAnn[,6])),]

row <- 0
signLVL <- 0.05/nrow(wholeblood)
tscores <- t(apply(wholeblood, 1, function(x){ # Using basis T statistics actually
  row <<- row + 1
  meanWB <- mean(x)
  sdWB <- sd(x)

  tNeutr  <- t.test(as.numeric(Neutr[row,]), x)
  sNeutr  <- tNeutr$statistic
  tBcell  <- t.test(as.numeric(Bcell[row,]), x) 
  sBcell  <- tBcell$statistic
  tTcell  <- t.test(as.numeric(Tcell[row,]), x)
  sTcell  <- tTcell$statistic
  tNKcell <- t.test(as.numeric(NKcell[row,]),x)
  sNKcell <- tNKcell$statistic
  tRBC    <- t.test(as.numeric(RBC[row,])   ,x)
  sRBC    <- tRBC$statistic

  if(tNeutr$p.value > signLVL)    sNeutr  <- NA
  if(tBcell$p.value > signLVL)    sBcell  <- NA
  if(tTcell$p.value > signLVL)    sTcell  <- NA
  if(tNKcell$p.value > signLVL)   sNKcell <- NA
  if(tRBC$p.value > signLVL)      sRBC    <- NA
  if(row %% 500 == 0)cat("Done: ", 5 * row," t-tests\n")
  return( cbind(sNeutr,sBcell,sTcell,sNKcell,sRBC) )
}))
colnames(tscores) <- c("Neutrophil", "Bcell", "Tcell", "NKcell", "RBC")
write.table(tscores, file="TScores5CellTypes_SignOnly.txt", quote = FALSE, sep='\t')

# Add annotation to the cell type vector
tscores <- tscores[which(rownames(tscores) %in% translation[,1]),]
sortZ <- match(rownames(tscores), translation[,1])
scores <- cbind(translation[sortZ, c(1,9)], tscores) # Add annotation

# Add annotation to the EGCUT vector
sortA <- match(IntrVec[,1], ivectorAnn[,6])
IntrVec <- cbind(ivectorAnn[sortA, 5:6], IntrVec)

# Merge the two vectors
IntrVec <- IntrVec[which(IntrVec[,1] %in% scores[,2]),]
sortT <- match(IntrVec[,1], scores[,2])

tscores_Annotated <- cbind(scores[sortT,], IntrVec[,c(1,4)])

# Bind the tscores_Annotated columns Neutrophil, Bcell and Tcell to the MetaAnalysis data
metaRes <- metaRes[which(rownames(metaRes) %in% tscores_Annotated[,1]),]     # match MetaRes
tscores_Annotated <- tscores_Annotated[which(tscores_Annotated[,1] %in% rownames(metaRes)),] # match ResVector
sortRes <- match(rownames(metaRes), tscores_Annotated[,1]) # Align

metaRes <- cbind(tscores_Annotated[sortRes, 3:7], metaRes)

corrs <- NULL  # Correlation analysis between 1:5 Neutro, B, Tcell, and all other columns
cnt <- 1
snpnames <- colnames(metaRes[, 6:length(metaRes)])
r <- apply(metaRes[, 6:length(metaRes)], 2, function(qtl){
  corrs <<- rbind(corrs, cor(qtl, metaRes[,1:5],use="pair", method="spearman"))
  rownames(corrs) <<- snpnames[1:cnt]
  cnt <<- cnt + 1
  if(cnt %% 100 == 0)cat("Cor of", cnt,"(QTL:cType):SNP pairs\n")
})

write.table(corrs, file="COR_MetaAnalysisZScore_CellTypeTScore_SignOnly.txt", quote = FALSE, sep='\t')

summary(corrs[,1]) # Neutrophil -0.6 to 0.6
summary(corrs[,2]) # Bcell      -0.4 to 0.4
summary(corrs[,3]) # Tcell      -0.6 to 0.6
summary(corrs[,4]) # NKcell     -0.5 to 0.5
summary(corrs[,5]) # RBC        -0.2 to 0.2


# Create some plots of the top 100 cell type specific expressions
top100_Neutro <- rownames(corrs)[sort(abs(corrs[,1]),index.ret=TRUE, decreasing=TRUE)$ix[1:100]]
top100_Bcell <- rownames(corrs)[sort(abs(corrs[,2]),index.ret=TRUE, decreasing=TRUE)$ix[1:100]]
top100_Tcell <- rownames(corrs)[sort(abs(corrs[,3]),index.ret=TRUE, decreasing=TRUE)$ix[1:100]]
top100_NKcell <- rownames(corrs)[sort(abs(corrs[,4]),index.ret=TRUE, decreasing=TRUE)$ix[1:100]]
top100_RBC <- rownames(corrs)[sort(abs(corrs[,5]),index.ret=TRUE, decreasing=TRUE)$ix[1:100]]

heatmap(corrs[c(top100_Neutro,top100_Bcell,top100_Tcell,top100_NKcell,top100_RBC), ], col=c("red", "white", "green"))

top <- which(apply(corrs,1,function(x){max(abs(x)) > 0.4}))
tree <- heatmap(corrs[top,], keep.dendro=TRUE)
heatmap.2(t(corrs[top[tree$rowInd],]), trace="none", col=c("red","pink","white","lightgreen","green"), main="Correlation CellTypes vs QTL Z-scores", key=FALSE, margins=c(5,10), scale="none", dendrogram="row", labCol="",xlab="Zscore Snp:Probe")

# Show the correlations between the annotated T vectors
cor(tscores_Annotated[,c(3, 4, 5, 7)], use="pair",method="spearman")

# Save the annotated Tvector
write.csv(tscores_Annotated[,-6], file="TScores5CellTypes_SignOnly_Annotated.txt", quote = FALSE)

inTrans <- which(rownames(metaRes) %in% rownames(translation))
metaRes <- metaRes[inTrans,]
sortTrans <- match(rownames(metaRes), rownames(translation)) # Align
mmRes <- cbind(translation[sortTrans,9], metaRes)

inSeq   <- which(as.character(mmRes[,1]) %in% rownames(RNASeq))
mmRes   <- mmRes[inSeq,]
sortSeq <- match(as.character(mmRes[,1]), rownames(RNASeq)) # Align

mmRes <- cbind(RNASeq[sortSeq,7], mmRes)

#NEUTRO GPL570 to SEBO RNA seq
ccrr <- function(Neutr,RNASeq,translation, col=7, type="Neutrophil"){
  inTrans <- which(rownames(Neutr) %in% translation[,1])
  Neutr <- Neutr[inTrans,]
  sortTrans <- match(rownames(Neutr), translation[,1]) # Align
  Neutr <- cbind(translation[sortTrans,9], Neutr)

  inRNASeq <- which(Neutr[,1] %in% rownames(RNASeq))
  Neutr <- Neutr[inRNASeq,]
  sortRNASeq <- match(Neutr[,1], rownames(RNASeq)) # Align
  Neutr <- cbind(RNASeq[sortRNASeq, col], Neutr)

  means <- apply(Neutr[,3:ncol(Neutr)],1,function(x){mean(as.numeric(x))})
  cat("Correlation:", cor(means, Neutr[,1],method="spearman"),"\n")
  name <- colnames(RNASeq)[col]
  png(paste0("plot",type,".png"),width=1024, height=800)
  plot(means, log2(Neutr[,1]),pch=19,cex=0.4,xlab=paste0("Affy GPL570 (",type,")"), ylab=paste0("RNA Seq (",name,")"))
  dev.off()
  means
}

nmean <- ccrr(Neutr,RNASeq,translation,7,"Neutrophil")
nmean <- ccrr(Neutr,RNASeq,translation,6,"Neutrophil_6")
nmean <- ccrr(Neutr,RNASeq,translation,5,"Neutrophil_5")
nmean <- ccrr(Neutr,RNASeq,translation,4,"Neutrophil_4")
nmean <- ccrr(Neutr,RNASeq,translation,3,"Neutrophil_3")
nmean <- ccrr(Neutr,RNASeq,translation,2,"Neutrophil_2")
nmean <- ccrr(Neutr,RNASeq,translation,1,"Neutrophil_1")


nmeans <- apply(Neutr[,3:ncol(Neutr)],1,function(x){mean(as.numeric(x))})
bmeans <- apply(Bcell[,3:ncol(Bcell)],1,function(x){mean(as.numeric(x))})
tmeans <- apply(Tcell[,3:ncol(Tcell)],1,function(x){mean(as.numeric(x))})
nkmeans <- apply(NKcell[,3:ncol(NKcell)],1,function(x){mean(as.numeric(x))})
rbcmeans <- apply(RBC[,3:ncol(RBC)],1,function(x){mean(as.numeric(x))})

MeanMatrix <- cbind(nmeans,bmeans,tmeans,nkmeans,rbcmeans)

colnames(MeanMatrix) <- c("Neutrophil(A)","Bcell(A)","Tcell(A)","nkCell(A)","RBC(A)")

inTrans <- which(rownames(MeanMatrix) %in% translation[,1])
MeanMatrix <- MeanMatrix[inTrans,]
sortTrans <- match(rownames(MeanMatrix), translation[,1]) # Align
MeanMatrix <- cbind(as.character(translation[sortTrans,9]), MeanMatrix)

inRNASeq <- which(MeanMatrix[,1] %in% rownames(RNASeq))
MeanMatrix <- MeanMatrix[inRNASeq,]
sortRNASeq <- match(MeanMatrix[,1], rownames(RNASeq)) # Align
MeanMatrix <- cbind(RNASeq[sortRNASeq, 1:7], MeanMatrix)

AffyRNASeqcorrs <- cor(apply(MeanMatrix[,-8],2,as.numeric))

write.csv(AffyRNASeqcorrs, file="COR_Affy_RNASeq.txt", quote = FALSE)

#plot(means, Neutr[,1],pch=19,cex=0.4,xlab="Affy GPL570 (Neutrophil)", ylab="RNA Seq (Granulocyte)")

