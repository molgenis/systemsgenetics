#
# create.tstat.matrix.R
#
# Copyright (c) 2013-2013 GBIC: Danny Arends, Harm-Jan Westra, Lude Franke and Ritsert C. Jansen
# last modified Jul, 2013
# first written Jun, 2013
#
setwd("~/Github/systemsgenetics/eqtl-mapping-pipeline/src/main/scripts")
source("helper.functions.R")

setwd("~/Github/Juha/")
# Load the Z scores from Meta Analysis
metaRes    <- read.harmjan.zscores()

# Load the cell type expression measured on GPL570
wholeblood <- read.csv("GPL570_WholeBlood.txt", sep='\t', row.names=1)
Neutr      <- read.csv("GPL570_Neutrophil_Good.txt", sep='\t', row.names=1)
Bcell      <- read.csv("GPL570_Bcell.txt", sep='\t', row.names=1)
Tcell      <- read.csv("GPL570_Tcell.txt", sep='\t', row.names=1)
NKcell     <- read.csv("GPL570_NKcell.txt", sep='\t', row.names=1)
RBC        <- read.csv("GPL570_RBC.txt", sep='\t', row.names=1)

# Load the probe annotation 'GPL570', remove probes not targetting a genes
translation <- read.affy.translation()

# Load the probe annotation for Illumina, remove probes not targetting a genes
setwd("/home/danny/Github/LudeNew")
ivectorAnn <- read.illumina.probes.information()
setwd("~/Github/Juha/")

# Load the vector from Lude and the interaction vector (from Harm-Jan)
LudeVec <- read.csv("NeutrophilVectorLude.txt",sep="\t")
IntrVec <- read.csv("2013-06-21-EGCUT-Vector-rs12057769-2000128.txt",sep='\t',row.names=NULL)

# Preprocessing remove uncorrelated arrays for the other cell-types
Bcell  <- Bcell[,-which(apply(cor(Bcell), 2, median) < 0.9)] # Use correlation for Bcell
Tcell  <- Tcell[,-which(apply(cor(Tcell), 2, median) < 0.9)] # Use correlation for Tcell

#Only take the interaction vector elements that match the Illumina annotation
IntrVec <- IntrVec[which(as.character(IntrVec[,1]) %in% as.character(ivectorAnn[,6])),]

if(!file.exists("tstat.matrix.affymetrix.txt")){
  row <- 0
  signLVL <- 0.05/nrow(wholeblood)

  tscores <- apply(wholeblood, 1, function(x){      # Using basis T statistics
    row <<- row + 1

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
    return( c(sNeutr,sBcell,sTcell,sNKcell,sRBC) )
  })
  colnames(tscores) <- c("Neutrophil", "Bcell", "Tcell", "NKcell", "RBC")
  write.table(tscores, file="tstat.matrix.affymetrix.txt", quote = FALSE, sep='\t')
}else{
  tscores <- read.csv("tstat.matrix.affymetrix.txt", sep='\t', row.names = 1)
}

scores <- annotate.affy.by.rownames(tscores, translation)

# Add annotation to the EGCUT vector
sortA <- match(IntrVec[,1], ivectorAnn[,6])
IntrVec <- cbind(ivectorAnn[sortA, 5:6], IntrVec)

# Merge the two vectors
IntrVec <- IntrVec[which(IntrVec[,1] %in% scores[,2]),]
sortT <- match(IntrVec[,1], scores[,2])

TscoresAnnotated <- cbind(scores[sortT,], IntrVec[,c(1,4)])

# Show the correlations between the annotated T vectors
cor(TscoresAnnotated[,c(3, 4, 5, 7)], use="pair",method="spearman")
# Save the annotated Tvector
write.csv(TscoresAnnotated[,-6], file="TScores5CellTypes_SignOnly_Annotated.txt", quote = FALSE)

# Bind the TscoresAnnotated columns Neutrophil, Bcell and Tcell to the MetaAnalysis data
metaRes <- add.metares.to.tmatrix(metaRes, TscoresAnnotated)

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

