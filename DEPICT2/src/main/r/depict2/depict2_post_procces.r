library(ggplot2)
library(gridExtra)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

source("depict2_functions.r")

path <- "/home/work/Desktop/depict2/maf_filtered/metabolites_2016_27005778_hg19_48/"

path <- "/home/work/Desktop/depict2/output/height_paper_v2/v56/"
files <- list.files(path, pattern="*.xlsx")

datasets <- list()
for (file in files) {
  name <- gsub("\\_hg19\\_enrichtments\\_exHla\\.xlsx", "", file)
  name <- gsub("\\_hg19\\_enrichtments\\_exHla\\_1\\.xlsx", "", name)
  name <- gsub("\\_enrichtments\\_exHla\\.xlsx", "", name)
  name <- gsub("\\_hg19\\.txt\\_exHla\\.xlsx", "", name)
  
  datasets[[name]] <- read.depict2(paste0(path, file))
}

# Manual 
dataset     <- "run_41_enrichtments_exHla.xlsx"
cur.dataset <- datasets[[dataset]]
p1          <- make.tsne.plot(cur.dataset$expression, dataset, x="Annotation1", y="Annotation2") +
  xlab("t-SNE component 1") +
  ylab("t-SNE component 2")

p1

pdf(width=10, height=10, file="/home/work/Desktop/depict2/plots/v56_tsne_plots_blood.pdf")
for (dataset in 1:length(datasets)) {
  cur.dataset <- datasets[[dataset]]
  p1 <- make.tsne.plot(cur.dataset$expression, paste0(names(datasets)[dataset], " expression"), x="Annotation1", y="Annotation2", limits=c(-6, 6))
  #p2 <- make.tsne.plot(cur.dataset$expression_scP3, paste0(names(datasets)[dataset], " expression_scP3"), x="Annotation2", y="Annotation3")
  #p3 <- make.tsne.plot(cur.dataset$expression_brain,paste0(names(datasets)[dataset], " expression_brain"), x="Annotation6", y="Annotation7")
  plot(p1)
  #grid.arrange(grobs=list(p1,p2,p3), ncol=3)
}
dev.off()


pdf(width=30, height=30, file="~/Desktop/depict2/plots/v56_correlation_heatmaps.pdf")
trait <- "Coregulation"
bla <- make.zscore.matrix(datasets, trait=trait)
make.correlation.heatmap(bla, trait)

#trait <- "Coregulation_brain"
#bla <- make.zscore.matrix(datasets, trait=trait)
#make.correlation.heatmap(bla, trait)

#trait <- "Coregulation_eQTLGen"
#bla <- make.zscore.matrix(datasets, trait=trait)
#make.correlation.heatmap(bla, trait)

trait <- "expression"
bla <- make.zscore.matrix(datasets, trait=trait)
make.correlation.heatmap(bla, trait)

#trait <- "expression_brain"
#bla <- make.zscore.matrix(datasets, trait=trait)
#make.correlation.heatmap(bla, trait)

#trait <- "expression_scP3"
#bla <- make.zscore.matrix(datasets, trait=trait)
#make.correlation.heatmap(bla, trait)

trait <- "KEGG"
bla <- make.zscore.matrix(datasets, trait=trait)
make.correlation.heatmap(bla, trait)

trait <- "HPO"
bla <- make.zscore.matrix(datasets, trait=trait)
make.correlation.heatmap(bla, trait)

trait <- "Reactome"
bla <- make.zscore.matrix(datasets, trait=trait)
make.correlation.heatmap(bla, trait)

dev.off()



############# Scratchpad

pbc <- datasets$primary_biliary_cirrhosis_2012_22961000
alz <- datasets$alzheimers_disease_2018_29777097

p1 <- make.tsne.plot(pbc$expression," expression", x="Annotation1", y="Annotation2")
p2 <- make.tsne.plot(alz$expression, " expression", x="Annotation1", y="Annotation2")
p3 <- theme.nature(plot.lm(pbc$expression$Enrichment.Z.score, alz$expression$Enrichment.Z.score, ylab="PBC zscore", xlab="ALZ zscore"))

grid.arrange(grobs=list(p1,p2,p3), ncol=3)


library(data.table)

x <- fread("gunzip -c ~/Desktop/depict2/Alzh_4_UK_Biobank_IGAP_17May2018_formatted.txt.gz", data.table=F)
rownames(x) <- make.unique(x$SNP)
y <- fread("gunzip -c ~/Desktop/depict2/primary_biliary_cirrhosis_2012_22961000_hg19.txt.gz", data.table=F)
rownames(y) <- make.unique(y[,1])
y <- y[-1,]

ol <- intersect(rownames(x), rownames(y))

c <- cor(-log10(x[ol,2]), -log10(y[ol,2]), use="complete.obs")

