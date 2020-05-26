#
# create.tstat.matrix.illumina.R
#
# Copyright (c) 2013-2013 GBIC: Danny Arends, Harm-Jan Westra, Lude Franke and Ritsert C. Jansen
# last modified Jul, 2013
# first written Jun, 2013
#
setwd("~/Github/systemsgenetics/eqtl-mapping-pipeline/src/main/scripts")
source("helper.functions.R")

setwd("~/Github/Juha/")
# Load the Illumina celltype data
Illu <- read.illumina.celltypes()
Illu <- annotate.illumina.celltypes(Illu)

WB <- read.illumina.wb()
WBAnnot <- add.illumina.probes.information(WB)

WBAnnot <- WBAnnot[which(rownames(WBAnnot) %in% rownames(Illu)),]
Illu <- Illu[which(rownames(Illu) %in% rownames(WBAnnot)),]

IlluSort <- match(rownames(WBAnnot), rownames(Illu))
Illu <- Illu[IlluSort,]

cellTypes <- unique(colnames(Illu))

if(!file.exists("tstat.matrix.illumina.txt")){
  signLVL <- 0.05/nrow(Illu)
  full <- NULL
  for(p in 1:nrow(Illu)){
    data <- NULL  
    for(celltype in cellTypes){
      cols <-  which(colnames(Illu)==celltype)
      tC  <- t.test(unlist(Illu[p,cols]), log2(unlist(WBAnnot[p,-1])))
      sC  <- tC$statistic
#      if(tC$p.value > signLVL) sC  <- NA
      data <- cbind(data, sC)
    }
    if(p %% 100 == 0) cat("Done",p,"Out of",nrow(Illu),"\n")
    full <- rbind(full,data)
    colnames(full) <- cellTypes
    rownames(full) <- rownames(Illu)[1:p]
  }
  write.table(full, file="tstat.matrix.illumina.txt", sep = '\t',quote=FALSE)
}else{
  full <- read.csv("tstat.matrix.illumina.txt", sep='\t', row.names = 1)
}

onlyonce <- which(apply(full,1,function(x){ sum(!is.na(x)) })==1)
onceCell <- full[onlyonce,]

