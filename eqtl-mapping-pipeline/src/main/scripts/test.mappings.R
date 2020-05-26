#
# test.mappings.R
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

WBillu <- read.illumina.wb()
WBillu <- add.illumina.probes.information(WBillu)
        
translation <- read.affy.translation()
WBaffy <- read.csv("GPL570_WholeBlood.txt", sep='\t', row.names=1)
WBaffy <- annotate.affy.by.rownames(WBaffy, translation)