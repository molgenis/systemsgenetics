#
# tstat.postprocessing.R
#
# Copyright (c) 2013-2013 GBIC: Danny Arends, Harm-Jan Westra, Lude Franke and Ritsert C. Jansen
# last modified Jul, 2013
# first written Jun, 2013
#
setwd("~/Github/systemsgenetics/eqtl-mapping-pipeline/src/main/scripts")
source("helper.functions.R")

setwd("~/Github/Juha/")

affyTrans <- read.affy.translation()

affymetrix <- read.csv("tstat.matrix.affymetrix.txt", sep='\t', row.names = 1)
illumina   <- read.csv("tstat.matrix.illumina.txt", sep='\t', row.names = 1)

getUnique <- function(matrix, nTypes = 1){
  onlyonce <- which(apply(matrix,1,function(x){ sum(!is.na(x)) }) == nTypes)
  return(matrix[onlyonce,])
}

affymetrix <- annotate.affy.by.rownames(affymetrix, affyTrans)
illumina <- add.illumina.probes.information(illumina)


affyUni <- getUnique(affymetrix)
illuUni <- getUnique(illumina)


IlluGenes <- apply(illuUni[,-1], 2, function(x){ illuUni[which(!is.na(x)),1] })
AffyGenes <- apply(affyUni[,-1], 2, function(x){ affyUni[which(!is.na(x)),1] })

