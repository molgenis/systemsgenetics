
projectDir <- "C:/Users/patri/Documents/GitHub/systemsgenetics/DEPICT2/src/main/r/depict2"
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")


library(ggsignif)
source(paste0(projectDir,"/depict2_functions.r"))

# Read reference datasets
ensembl <- read.table("ensembl_gene_position_b37_biotype.txt",sep="\t", header=T, row.names = 1, stringsAsFactors = F)
ensembl$gene.length = ensembl$Gene.end..bp. - ensembl$Gene.start..bp.
rownames(ensembl) <- make.names(ensembl$Gene.name, unique=T)

# Gnomad Pli
gnomad <- read.table("gnomad.v2.1.1.lof_metrics.by_gene.txt",sep="\t", header=T, stringsAsFactors = F)
gnomad <- gnomad[!duplicated(gnomad$gene),]
rownames(gnomad) <- make.names(gsub("\\.\\d+", "", ensembl[gnomad$gene, ]$Gene.stable.ID.version), unique=T)

# HPO genes
hpo.genes <- read.table("genes_to_phenotype.txt", stringsAsFactors = F, sep="\t")
hpo.genes <- gsub("\\.\\d+", "", ensembl[unique(hpo.genes$V2), "Gene.stable.ID.version"])

path <- "excelsForEnrichment"

# Coregulation
files    <- list.files(path, pattern="*.xlsx", full.names = T)
datasets <- list()
for (file in files) {
  name <- gsub("\\_hg19\\_enrichtments\\_exHla\\.xlsx", "", basename(file))
  name <- gsub("\\_hg19\\_enrichtments\\_exHla\\_1\\.xlsx", "", name)
  name <- gsub("\\_enrichtments\\_exHla\\.xlsx", "", name)
  name <- gsub("\\_hg19\\.txt\\_exHla\\.xlsx", "", name)
  
  if (length(grep("v55", file)) > 0) {
    name <- paste0(name, "_v55")
  }
  datasets[[name]] <- read.depict2(file)
}

# Genepvalues
files <- list.files(path, pattern="*_genePvalues.txt", full.names = T)
genep <- read.enrichments(files)
genep[is.na(genep)] <- 1


# HPO genes and Pli

hpo.genes             <- intersect(hpo.genes, rownames(gnomad))
bonf.sig.gwas.genes   <- rownames(genep)[rowSums(genep < (0.05 / (nrow(genep) * ncol(genep)))) >=1]

bonf.sig.coreg.genes <- unique(unlist(lapply(datasets[grep("run", names(datasets), invert=T, value=T)], function(dataset){
  coreg <- dataset$Coregulation
  return(coreg$Gene.set[coreg$Bonferroni.significant])
})))

others                <- rownames(gnomad)[!rownames(gnomad) %in% unique(c(hpo.genes, bonf.sig.coreg.genes))]

df.plot <- data.frame(pli=c(gnomad[others,"pLI"],
                            gnomad[bonf.sig.coreg.genes,"pLI"],
                            gnomad[bonf.sig.gwas.genes, "pLI"],
                            gnomad[hpo.genes, "pLI"]),
                      annot=factor(c(rep(paste0("Not core|HPO gene N=", length(others)), length(others)), 
                              rep(paste0("Core genes N=", length(bonf.sig.coreg.genes)), length(bonf.sig.coreg.genes)),
                              rep(paste0("GWAS genes N=", length(bonf.sig.gwas.genes)), length(bonf.sig.gwas.genes)),
                              rep(paste0("HPO genes N=", length(hpo.genes)), length(hpo.genes))),
                              levels=c(paste0("Core genes N=", length(bonf.sig.coreg.genes)),
                                       paste0("HPO genes N=", length(hpo.genes)),
                                       paste0("GWAS genes N=", length(bonf.sig.gwas.genes)),
                                       rep(paste0("Not core|HPO gene N=", length(others))))))

df.plot2 <- data.frame(pli=c(gnomad[others,"mis_z"],
                            gnomad[bonf.sig.coreg.genes,"mis_z"],
                            gnomad[bonf.sig.gwas.genes, "mis_z"],
                            gnomad[hpo.genes, "mis_z"]),
                      annot=factor(c(rep(paste0("Not core|HPO gene N=", length(others)), length(others)), 
                                     rep(paste0("Core genes N=", length(bonf.sig.coreg.genes)), length(bonf.sig.coreg.genes)),
                                     rep(paste0("GWAS genes N=", length(bonf.sig.gwas.genes)), length(bonf.sig.gwas.genes)),
                                     rep(paste0("HPO genes N=", length(hpo.genes)), length(hpo.genes))),
                                   levels=c(paste0("Core genes N=", length(bonf.sig.coreg.genes)),
                                            paste0("HPO genes N=", length(hpo.genes)),
                                            paste0("GWAS genes N=", length(bonf.sig.gwas.genes)),
                                            rep(paste0("Not core|HPO gene N=", length(others))))))



pdf(width=8, height=3, file="pli_comparrison_all.pdf")
p <- ggplot(df.plot, aes(y=pli, x=annot, fill=annot)) +
geom_violin(color="white", scale="width") +
  geom_boxplot(width=0.05, color="black") +
  ylab("Gnomad Pli score (LoF intolerance)") +
  xlab("")
  t.test(gnomad[others,"pLI"], gnomad[bonf.sig.coreg.genes,"pLI"]) 
  t.test(gnomad[hpo.genes,"pLI"], gnomad[bonf.sig.coreg.genes,"pLI"]) 
theme.nature(p) + scale_fill_manual(values=c("dodgerblue3", "#3BB273", "goldenrod2", "#8576B6"))

p <- ggplot(df.plot2, aes(y=pli, x=annot, fill=annot)) +
  geom_violin(color="white", scale="width") +
  geom_boxplot(width=0.05, color="black") +
  ylab("Gnomad misZ score") +
  xlab("")
t.test(gnomad[others,"mis_z"], gnomad[bonf.sig.coreg.genes,"mis_z"]) 
t.test(gnomad[hpo.genes,"mis_z"], gnomad[bonf.sig.coreg.genes,"mis_z"]) 
theme.nature(p) + scale_fill_manual(values=c("dodgerblue3", "#3BB273", "goldenrod2", "#8576B6"))


dev.off()

sum(gnomad[hpo.genes,"pLI"] >= 0.9, na.rm = T)
sum(gnomad[bonf.sig.coreg.genes,"pLI"] >= 0.9, na.rm = T )

sum(gnomad[hpo.genes,"pLI"] < 0.9, na.rm = T)
sum(gnomad[bonf.sig.coreg.genes,"pLI"] < 0.9, na.rm = T )

(x <- matrix(c(490, 1773, 102, 241), nrow =2, byrow = T))
(x <- matrix(c(102, 241, 490, 1773), nrow =2))


mis_z             

fisher.test(x)
