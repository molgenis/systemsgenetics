library(readxl)
library(ggplot2)
source("~/Documents/projects/pr_integration/wd_integration/Code/UtillityScripts/PlottingFunctions.r")
library(gridExtra)


ced <- data.frame(read_excel("/home/work/Documents/tmp/MS_enrichtments_exHla_patrick.xlsx", sheet="Coregulation"), row.names=1)
ced.ichip <- data.frame(read_excel("/home/work/Documents/tmp/MS_enrichtments_exHla.xlsx", sheet="Coregulation_eQTLGen"), row.names = 1)
p <- theme.nature(plot.lm(ced.ichip[rownames(ced),]$Enrichment.Z.score,
                          ced$Enrichment.Z.score,
                          xlab=paste0("MS Zscore ", "Coregulation"),
                          ylab=paste0("MS Zscore ", "Coregulation eQTLGen"), horiz=F, fix.axes=T))
p


ced <- data.frame(read_excel("/home/work/Documents/tmp/SLE_enrichtments_exHla_4_patrick.xlsx", sheet="Coregulation"), row.names=1)
ced.ichip <- data.frame(read_excel("/home/work/Documents/tmp/SLE_enrichtments_exHla.xlsx", sheet="Coregulation_eQTLGen"), row.names = 1)
p <- theme.nature(plot.lm(ced.ichip[rownames(ced),]$Enrichment.Z.score,
                          ced$Enrichment.Z.score,
                          xlab=paste0("SLE Zscore ", "Coregulation"),
                          ylab=paste0("SLE Zscore ", "Coregulation eQTLGen"), horiz=F, fix.axes=T))
p
