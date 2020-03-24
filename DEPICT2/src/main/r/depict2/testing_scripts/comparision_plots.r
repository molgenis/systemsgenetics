library(readxl)

# Loading data
ensembl <- unique(read.table("~/Documents/data/reference/ensembl/ensembl_gene_position_export.bed", stringsAsFactors = F)[,c(4,5)])
rownames(ensembl) <- make.names(ensembl[,2], unique=T)

# Function definitions
simpleQQPlot = function (observedPValues) {
  plot(-log10(1:length(observedPValues)/length(observedPValues)), 
       -log10(sort(observedPValues)))
  abline(0, 1, col = "red")
}

# ------------------------------------------------------ 
read.depict2 <- function(path) {
  potential_traits <- c("Coregulation")
  output <- list()
  for (sheet in potential_traits) {
    tmp <- tryCatch({data.frame(read_excel(path, sheet=sheet, col_types ="guess", trim_ws = T), stringsAsFactors=F)},
                    error=function(a){return(NA)},
                    warn=function(a){return(NA)})
    if(!is.na(tmp)) {
      for (i in 1:ncol(tmp)) {
        if (class(tmp[,i]) == "character"){
          tmp[,i] <- type.convert(tmp[,i], as.is=T)
        }
      }
      rownames(tmp) <- tmp[,1]
      output[[sheet]] <- tmp
    }
  }
  return(output)
}


# -----------------------------------------------------------------------

# Height HPO AUC
x <- read.table("~/Desktop/depict2/hpo_predictions/height_2018_30124842_hg19_coreGene_hpoAUC_hpo.txt", header=T, row.names = 1)
y <- read.table("~/Desktop/depict2/output/height_paper/height_2018_30124842_hg19_coreGene_hpoAUC_hpo.txt", header=T, row.names=1)
z <- read.table("~/Desktop/depict2/output/height_paper_no_gc/height_2018_30124842_hg19_coreGene_hpoAUC_hpo.txt", header=T, row.names=1)

# height coreg
x.coreg <- read.depict2("~/Desktop/depict2/output/maf_filtered/v51/height_2018_30124842_hg19_enrichtments_exHla.xlsx")$Coregulation
y.coreg <- read.depict2("~/Desktop/depict2/output/height_paper/height_2018_30124842_hg19_enrichtments_exHla.xlsx")$Coregulation
z.coreg <- read.depict2("~/Desktop/depict2/output/height_paper_no_gc/height_2018_30124842_hg19_enrichtments_exHla.xlsx")$Coregulation

# Overlapping HPO terms with old and new results
ol <- intersect(rownames(x), rownames(y))

# Remove zeroes
x[x[,1]==0,] <- NA
y[y[,1]==0,] <- NA
z[z[,1]==0,] <- NA

x <- x[!is.na(x[,1]),,drop=F]
y <- y[!is.na(y[,1]),,drop=F]
z <- z[!is.na(z[,1]),,drop=F]

# -----------------------------------------------------------------------
# Boxplots
par(mfrow=c(1,2))
ol <- intersect(rownames(x), rownames(y))
boxplot(x[ol,1], y[ol,1], z[ol, 1],
        main="HPO AUC predictions for coregulation",
        names = c("Skewness 0.25", "Protein coding with GC", "Protein coding no GC"),
        ylab="HPO AUC")

ol <- intersect(rownames(x.coreg), rownames(y.coreg))
boxplot(x.coreg[ol,"Enrichment.Z.score"],
        y.coreg[ol,"Enrichment.Z.score"],
        z.coreg[ol,"Enrichment.Z.score"],
        main="Core gene score comparison",
        names = c("Skewness 0.25", "Protein coding with GC", "Protein coding no GC"),
        ylab="Core gene Zscore")

# -----------------------------------------------------------------------
# Scatterplots
par(mfrow=c(2,2))
ol <- intersect(rownames(x), rownames(y))
plot(x[ol,1],
     y[ol,1],
     xlab="HPO AUC Skewness 0.25 1000M rescue",
     ylab="HPO AUC Protein coding 10M rescue with GC",
     ylim=c(0,1), xlim=c(0,1),
     main=paste0("Pearson R: ", format(cor(x[ol,1], y[ol,1], use="complete.obs"), digits=2)))
abline(a=0, b=1, col="red")

plot(y[ol,1],
     z[ol,1],
     xlab="HPO AUC Protein coding 10M rescue with GC",
     ylab="HPO AUC Protein coding 10M rescue no GC",
     ylim=c(0,1), xlim=c(0,1),
     main=paste0("Pearson R: ", format(cor(y[ol,1], z[ol,1], use="complete.obs"), digits=2)))
abline(a=0, b=1, col="red")

ol <- intersect(rownames(x.coreg), rownames(y.coreg))
plot(x.coreg[ol,"Enrichment.Z.score"],
     y.coreg[ol,"Enrichment.Z.score"],
     xlab="Core gene score Skewness 0.25 1000M rescue",
     ylab="Core gene score Protein coding 10M rescue with GC",
     main=paste0("Pearson R: ", format(cor(x.coreg[ol,"Enrichment.Z.score"], y.coreg[ol,"Enrichment.Z.score"], use="complete.obs"), digits=2)))
abline(a=0, b=1, col="red")

plot(y.coreg[ol,"Enrichment.Z.score"],
     z.coreg[ol,"Enrichment.Z.score"],
     xlab="Core gene score Protein coding 10M rescue with GC",
     ylab="Core gene score Protein coding 10M rescue no GC",
     main=paste0("Pearson R: ", format(cor(y.coreg[ol,"Enrichment.Z.score"], z.coreg[ol,"Enrichment.Z.score"], use="complete.obs"), digits=2)))
abline(a=0, b=1, col="red")

# -----------------------------------------------------------------------
# Kidney vs height vs blood

par(mfrow=c(1,3))
x.coreg <- read.depict2("~/Desktop/depict2/output/height_paper/glomerular_filtration_rate_EUR_2019_31152163_hg19_enrichtments_exHla.xlsx")$Coregulation
y.coreg <- read.depict2("~/Desktop/depict2/output/height_paper/height_2018_30124842_hg19_enrichtments_exHla.xlsx")$Coregulation
x.coreg <- read.depict2("~/Desktop/depict2/output/height_paper_no_gc/vucovik_cellcounts_enrichtments_lymph_p_exHla.xlsx")$Coregulation

ol <- intersect(rownames(x.coreg), rownames(y.coreg))
plot(x.coreg[ol,"Enrichment.Z.score"],
     y.coreg[ol,"Enrichment.Z.score"],
     xlab="Lymp p core gene score",
     ylab="Height core gene score",
     main=paste0("Core gene Pearson R: ", format(cor(x.coreg[ol,"Enrichment.Z.score"], y.coreg[ol,"Enrichment.Z.score"], use="complete.obs"), digits=2)))
abline(a=0, b=1, col="red")

x <- read.table("~/Desktop/depict2/output/height_paper/glomerular_filtration_rate_EUR_2019_31152163_hg19_coreGene_hpoAUC_hpo.txt", header=T, row.names = 1)
y <- read.table("~/Desktop/depict2/output/height_paper/height_2018_30124842_hg19_coreGene_hpoAUC_hpo.txt", header=T, row.names=1)

ol <- intersect(rownames(x), rownames(y))

# Remove zeroes
x[x[,1]==0,] <- NA
y[y[,1]==0,] <- NA

x <- x[!is.na(x[,1]),,drop=F]
y <- y[!is.na(y[,1]),,drop=F]

plot(x[ol,1],
     y[ol,1],
     xlab="GFR HPO AUC",
     ylab="Height HPO AUC",
     main=paste0("HPO AUC Pearson R: ", format(cor(x[ol,1], y[ol,1], use="complete.obs"), digits=2)))
abline(a=0, b=1, col="red")

x <- read.table("~/Desktop/depict2/output/height_paper/glomerular_filtration_rate_2019_31152163_hg19_genePvalues.txt", header=T, row.names = 1)
y <- read.table("~/Desktop/depict2/output/height_paper/height_2018_30124842_hg19_genePvalues.txt", header=T, row.names=1)

ol <- intersect(rownames(x), rownames(y))

# Remove zeroes
x <- x[!is.na(x[,1]),,drop=F]
y <- y[!is.na(y[,1]),,drop=F]

x[x[,1]==0,] <- NA
y[y[,1]==0,] <- NA

x <- x[!is.na(x[,1]),,drop=F]
y <- y[!is.na(y[,1]),,drop=F]

plot(-log10(x[ol,1]),
     -log10(y[ol,1]),
     xlab="-log10(GFR gene p)",
     ylab="-log10(Height gene p )",
     main=paste0("Gene P Pearson R: ", format(cor(x[ol,1], y[ol,1], use="complete.obs"), digits=2)))
abline(a=0, b=1, col="red")



# ---------------------------------------------
# Gene Pvalue plots
x <- read.table("~/Desktop/depict2/output/height_paper/height_2018_30124842_hg19_genePvalues.txt", header=T, row.names=1)
y <- read.table("~/Desktop/depict2/output/height_paper_no_gc/height_2018_30124842_hg19_genePvalues.txt", header=T, row.names=1)
#z <- read.table("~/Desktop/depict2/output/", header=T, stringsAsFactors = F)
#rownames(z) <- make.names(ensembl[z$gene_symbol,1], unique=T)
ol <- intersect(rownames(x), rownames(y))

par(mfrow=c(1,1))
plot(-log10(x[ol, 1]),
     -log10(y[ol, 1]),
     xlab="-log10(Pascal gene Pvalue)",
     ylab="-log10(Depict gene Pvalue)",
     main="Pascal vs Depict gene P comparison height 2018")
abline(a=0, b=1, col="red",)
