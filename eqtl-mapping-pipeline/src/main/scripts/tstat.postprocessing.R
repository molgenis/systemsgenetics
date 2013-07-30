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

affyUni <- getUnique(affymetrix)
illuUni <- getUnique(illumina)

affyUni <- annotate.affy.by.rownames(affyUni, affyTrans)
illuUni <- add.illumina.probes.information(illuUni)
