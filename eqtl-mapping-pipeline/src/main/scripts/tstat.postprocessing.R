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
ProbeAnnotation <- read.csv("GPL6102-11574.txt",sep="\t",skip=27) # Probe annotation
IlluminaMapping <- read.csv("2013-07-18-ProbeAnnotationFile.txt",sep='\t')
#ID = IlluminaArrayAdress HT12V3
CCM <- read.csv("CellCountCorrelationMatrix.txtEndophenotypeVsExpressionCorrelationMatrix.txt")

affyTrans <- read.affy.translation()

affymetrix <- read.csv("tstat.matrix.affymetrix.txt", sep='\t', row.names = 1)
illumina   <- read.csv("tstat.matrix.illumina.txt", sep='\t', row.names = 1)

arrayIDs   <- illuProbeToArrayID(rownames(illumina), ProbeAnnotation)
HUGOsIllu  <- ArrayIdToHugo(arrayIDs, IlluminaMapping, "HumanWG.6v3.txt")
HUGOsCCM   <- ArrayIdToHugo(rownames(CCM), IlluminaMapping)

#getUnique <- function(matrix, nTypes = 1){
#  onlyonce <- which(apply(matrix,1,function(x){ sum(!is.na(x)) }) == nTypes)
#  return(matrix[onlyonce,])
#}

affymetrix <- annotate.affy.by.rownames(affymetrix, affyTrans)
illumina <- add.illumina.probes.information(illumina)

#affyUni <- getUnique(affymetrix)
#illuUni <- getUnique(illumina)


IlluGenes <- apply(illuUni[,-1], 2, function(x){ illuUni[which(!is.na(x)),1] })
AffyGenes <- apply(affyUni[,-1], 2, function(x){ affyUni[which(!is.na(x)),1] })

