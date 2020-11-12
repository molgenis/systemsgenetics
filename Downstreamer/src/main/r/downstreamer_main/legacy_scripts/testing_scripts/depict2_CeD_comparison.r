library(readxl)
library(ggplot2)
source("~/Documents/projects/pr_integration/wd_integration/Code/UtillityScripts/PlottingFunctions.r")
library(gridExtra)

make.plot <- function(trait) {
  ced <- data.frame(read_excel("/home/work/Documents/tmp/CeD_Z-scores_enrichtments_exHla.xlsx", sheet=trait), row.names=1)
  ced.ichip <- data.frame(read_excel("/home/work/Documents/tmp/CeD_ichip_Z-scores_enrichtments_exHla.xlsx", sheet=trait), row.names = 1)
  p <- theme.nature(plot.lm(ced.ichip[rownames(ced),]$Enrichment.Z.score,
                       ced$Enrichment.Z.score,
                       xlab=paste0("CeD Zscore ", trait),
                       ylab=paste0("CeD ichip Zscore ", trait), horiz=F, fix.axes=T))
  return(p)
}

plots <- list()
plots[[1]] <- make.plot("Coregulation")
plots[[2]] <- make.plot("Reactome")
plots[[3]] <- make.plot("GO_P")
plots[[4]] <- make.plot("GO_C")
plots[[5]] <- make.plot("GO_F")
plots[[6]] <- make.plot("KEGG")
plots[[7]] <- make.plot("HPO")
plots[[8]] <- make.plot("gtex")
plots[[9]] <- make.plot("eigen")
plots[[10]] <- make.plot("eigen_eQTLGen")

grid.arrange(grobs=plots)



make.plot <- function(trait) {
  ced <- data.frame(read_excel("/home/work/Documents/tmp/CeD_Z-scores_enrichtments_exHla.xlsx", sheet=trait), row.names=1)
  ced.ichip <- data.frame(read_excel("~/Downloads/CeD_enrichtments_exHla_4.xlsx", sheet=trait), row.names = 1)
  p <- theme.nature(plot.lm(ced.ichip[rownames(ced),]$Enrichment.Z.score,
                            ced$Enrichment.Z.score,
                            xlab=paste0("CeD Zscore ", trait),
                            ylab=paste0("CeD  Patrick Zscore ", trait), horiz=F, fix.axes=T))
  return(p)
}

plots <- list()
plots[[1]] <- make.plot("Coregulation")
plots[[2]] <- make.plot("Reactome")
plots[[3]] <- make.plot("GO_P")
plots[[4]] <- make.plot("GO_C")
plots[[5]] <- make.plot("GO_F")
plots[[6]] <- make.plot("KEGG")
plots[[7]] <- make.plot("HPO")
plots[[8]] <- make.plot("gtex")
plots[[9]] <- make.plot("eigen")
plots[[10]] <- make.plot("eigen_eQTLGen")

grid.arrange(grobs=plots)


make.plot <- function(trait) {
  ced <- data.frame(read_excel("/home/work/Documents/tmp/CeD_Z-scores_enrichtments_exHla.xlsx", sheet="Coregulation"), row.names=1)
  ced.ichip <- data.frame(read_excel("/home/work/Documents/tmp/CeD_Z-scores_enrichtments_exHla_2.xlsx", sheet="Coregulation_eQTLGen"), row.names = 1)
  p <- theme.nature(plot.lm(ced.ichip[rownames(ced),]$Enrichment.Z.score,
                            ced$Enrichment.Z.score,
                            xlab=paste0("CeD Zscore ", "Coregulation"),
                            ylab=paste0("CeD eqtlgen coreg Zscore ", "Coregulation eQTLGen"), horiz=F, fix.axes=T))
  return(p)
}

grid.arrange(grobs=plots)


