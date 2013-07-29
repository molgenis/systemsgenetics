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
WBMean  <- apply(AllAnnot[,-1], 1, mean)

