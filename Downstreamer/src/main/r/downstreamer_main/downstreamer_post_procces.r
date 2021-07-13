source("downstreamer_functions.r")

# Read excel files
path  <- "/home/work/Desktop/depict2/output/final_paper/excels/"
files <- list.files(path, pattern=".*_enrichtments.*.xlsx")

datasets <- list()
for (file in files) {
  name <- gsub("\\_enrichtments\\_exHla\\.xlsx", "", file)
  name <- gsub("\\_enrichtments\\_exHla\\_1\\.xlsx", "", name)
  name <- gsub("\\_enrichtments\\_exHla\\.xlsx", "", name)
  name <- gsub("\\.txt\\_exHla\\.xlsx", "", name)
  
  datasets[[name]] <- read.depict2(paste0(path, file))
}


path  <- "/home/work/Desktop/depict2/output/simulated_gwas//"
files <- list.files(path, pattern=".*_enrichtments.*_1.xlsx")
files <- c(files, list.files(path, pattern=".*_enrichtments.*.xlsx"))

for (file in files) {
  name <- gsub("\\_enrichtments\\_exHla\\.xlsx", "", file)
  name <- gsub("\\_enrichtments\\_exHla\\_1\\.xlsx", "", name)
  name <- gsub("\\_enrichtments\\_exHla\\.xlsx", "", name)
  name <- gsub("\\.txt\\_exHla\\.xlsx", "", name)
  
  datasets[[name]] <- read.depict2(paste0(path, file))
}

save(datasets, file= "depict2/output/downstreamer_v75_results_cache_2020-10-29.RData")
load("depict2/output/downstreamer_v75_results_cache_2020-10-29.RData")

#----------------------------------------------------------------------
# T-SNE plots

pdf(width=10, height=10, file="/home/work/Desktop/depict2/plots/v75_tsne_plots_blood.pdf")
for (dataset in 1:length(datasets)) {
  cur.dataset <- datasets[[dataset]]
  p1 <- make.tsne.plot(cur.dataset$expression, paste0(names(datasets)[dataset], " expression"), x="Annotation1", y="Annotation2", limits=c(-6, 6))
  #p2 <- make.tsne.plot(cur.dataset$expression_scP3, paste0(names(datasets)[dataset], " expression_scP3"), x="Annotation2", y="Annotation3")
  #p3 <- make.tsne.plot(cur.dataset$expression_brain,paste0(names(datasets)[dataset], " expression_brain"), x="Annotation6", y="Annotation7")
  plot(p1)
  #grid.arrange(grobs=list(p1,p2,p3), ncol=3)
}
dev.off()

#----------------------------------------------------------------------
# Correlation heatmaps

pdf(width=15, height=15, file="~/Desktop/depict2/plots/v75_correlation_heatmaps.pdf")
trait <- "Coregulation"
bla <- make.zscore.matrix(datasets, trait=trait)
make.correlation.heatmap(bla, trait)

trait <- "eigenvectors_1588"
bla <- make.zscore.matrix(datasets, trait=trait)
make.correlation.heatmap(bla, trait)

trait <- "expression"
bla <- make.zscore.matrix(datasets, trait=trait)
make.correlation.heatmap(bla, trait)

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

#----------------------------------------------------------------------
# qq plots

plots.coregulation <- lapply(names(datasets), function(dataset) {
    cur.dataset <- datasets[[dataset]]$Coregulation$Enrichment.P.value
  
  return(theme.nature(fancy.qq.plot(cur.dataset, main=dataset)))
})

png(width=2500, height=2500, file="~/Desktop/depict2/plots/v75_qqplots_coregulation.png")
grid.arrange(grobs=plots.coregulation)
dev.off()


plots.expression <- lapply(names(datasets), function(dataset) {
  cur.dataset <- datasets[[dataset]]$expression$Enrichment.P.value
  
  return(theme.nature(fancy.qq.plot(cur.dataset, main=dataset)))
})

png(width=2500, height=2500, file="~/Desktop/depict2/plots/v75_qqplots_expression.png")
grid.arrange(grobs=plots.expression)
dev.off()



  for(dataset in names(datasets)) {
  cur.dataset <- datasets[[dataset]]$Coregulation$Enrichment.P.value
  hist(cur.dataset)
}

