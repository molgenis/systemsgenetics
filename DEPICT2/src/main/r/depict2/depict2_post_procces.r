library(ggplot2)
library(gridExtra)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

# Function definitions
# ------------------------------------------------------ 
# NPG friendly ggplot theme
theme.nature <- function(p, base_size = 11, base_family = "ArialMT") {
  p <- p + theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black", size=0.5),
          axis.ticks = element_line(size=0.5),
          axis.text = element_text(size=base_size, family="ArialMT", face="plain"),
          strip.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size=base_size, family = "ArialMT", face="plain"),
          complete = TRUE,
          plot.title = element_text(hjust=0.5))
  return(p)
}
# ------------------------------------------------------ 
# Plot depict2 scatterplot by enrichment score
make.tsne.plot <- function(data, trait, x="Annotation1", y="Annotation2", colour="Enrichment.Z.score") {
  data <- data[order(abs(data[,colour])),]
  # low  <- colorRampPalette(c('firebrick3', '#F6DCDC'))
  # mid  <- colorRampPalette(c("#F6DCDC", "#DDE6EE"))
  # high <- colorRampPalette(c('#DDE6EE', 'dodgerblue4'))
  # bins= 100
  # breaks = seq(-15, 15, length.out = bins)
  # colRamp <- c(low((bins-10)/2), mid(10), high((bins-10)/2))
  # cols <- colRamp[as.numeric(cut(data[,colour], breaks = breaks))]
  # 
  # pdf(width=7.5, height=7.5, file="~/Desktop/")
  # 
  # 
  # layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
  # plot(data$Annotation1, data$Annotation2, pch = 20, col=adjustcolor(cols, alpha.f = 0.8), bty='n', xlab="t-SNE component 1", ylab="t-SNE compenent 2")
  # 
  # legend_image <- as.raster(matrix(colRamp, ncol=1))
  # plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-score')
  # text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
  # rasterImage(legend_image, 0, 0, 1,1)
  # 
  #  legend("topright",legend=c("15", "0", "-15"),fill=c("dodgerblue4", "grey", "firebrick3"), title="Z-score", border=F, bty="n")
  # #legend.scale(c(-15, 15), col=cols)
  # 
  # p <- ggplot(aes(x=as.numeric(data[,x]),
  #                 y=as.numeric(data[,y]),
  #                 fill=data[,colour]),
  #             data=data) +
  #   geom_point(color=cols) +
  #   labs(title=trait) +
  #   xlab(x) +
  #   ylab(y) +
  #   scale_fill_continuous(name="Z-score", limits=c(-15, 15), low="firebrick3", mid="grey", high="dodgerblue4")
  
 p <- ggplot(aes(x=as.numeric(data[,x]), y=as.numeric(data[,y])), data=data) +
    labs(title=trait) + xlab(x) + ylab(y)

 p <- ggplot(aes(x=as.numeric(data[,x]), y=as.numeric(data[,y]), colour=data[,colour]), data=data) +
    geom_point() +
    labs(title=trait) + xlab(x) + ylab(y) +
   scale_colour_gradient(low=adjustcolor("red", alpha.f = 0.5),
                         high=adjustcolor("blue", alpha.f = 0.5),
                         limits=c(-15, 15),
                         name="Z-score")
  
  p <- theme.nature(p)
  return(p)
}


# ------------------------------------------------------ 
read.depict2 <- function(path) {
  potential_traits <- c("expression","expression_scP3","expression_brain","expression_gs_tcell","Coregulation","Coregulation_eQTLGen","Coregulation_brain","Reactome","GO_P","GO_C","GO_F","KEGG","HPO","GO_P_brain","GO_C_brain","GO_F_brain","KEGG_brain","HPO_brain","Reactome_brain","gtexV8","gtex")
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

# ------------------------------------------------------ 
make.zscore.matrix <- function(datasets, trait="Coregulation") {
  out <- matrix()
  i=0
  for (dataset in names(datasets)) {
    tmp <- datasets[[dataset]][[trait]]
    if (i==0) {
      out <- matrix(tmp$Enrichment.Z.score)
      rownames(out) <- tmp[,1]
    } else {
      out <- cbind(out, tmp[rownames(out),]$Enrichment.Z.score)
    }
    i <- i+1
  }
  colnames(out) <- names(datasets)
  
  return(out)
}

# ------------------------------------------------------ 
cor.test.p <- function(x, method){
  FUN <- function(x, y) cor.test(x, y, method=method, use="complete.obs")[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

# ------------------------------------------------------ 
make.correlation.heatmap <- function(matrix, trait) {
  c <- cor(matrix, method="spearman", use="pairwise.complete.obs")
  #p <- cor.test.p(matrix, method="spearman")
  
  #c[p > 0.05] <- 0
  
  hm(c, main=trait)
}

# ------------------------------------------------------ 
hm <- function(data, cellwidth=12, cellheight=12, limit=NULL, ...) {
  break.list <- seq(-max(abs(data)), max(abs(data)), by=max(abs(data))/100)
  pheatmap(data,
           breaks=break.list,
           col=colorRampPalette(rev(brewer.pal(n=7, name ="RdBu")))(length(break.list)),
           cellwidth=cellwidth,
           cellheight=cellheight,
           ...)
}

# ------------------------------------------------------ 

path <- "/home/work/Desktop/depict2/maf_filtered/metabolites_2016_27005778_hg19_48/"

path <- "/home/work/Desktop/depict2/output/maf_filtered/v51/"
files <- list.files(path)

datasets <- list()
for (file in files) {
  name <- gsub("\\_hg19\\_enrichtments\\_exHla\\.xlsx", "", file)
  datasets[[name]] <- read.depict2(paste0(path, file))
}

# Manual 
dataset <- "educational_attainment_2018_30038396"
p1 <- make.tsne.plot(cur.dataset$expression, "Educational attainment", x="Annotation1", y="Annotation2") +
  xlab("t-SNE component 1") +
  ylab("t-SNE component 2")

p1

pdf(width=10, height=10, file="/home/work/Desktop/depict2/plots/v51_tsne_plots_blood.pdf")
for (dataset in 1:length(datasets)) {
  cur.dataset <- datasets[[dataset]]
  p1 <- make.tsne.plot(cur.dataset$expression, paste0(names(datasets)[dataset], " expression"), x="Annotation1", y="Annotation2")
  #p2 <- make.tsne.plot(cur.dataset$expression_scP3, paste0(names(datasets)[dataset], " expression_scP3"), x="Annotation2", y="Annotation3")
  #p3 <- make.tsne.plot(cur.dataset$expression_brain,paste0(names(datasets)[dataset], " expression_brain"), x="Annotation6", y="Annotation7")
  plot(p1)
  #grid.arrange(grobs=list(p1,p2,p3), ncol=3)
}
dev.off()


pdf(width=30, height=30, file="~/Desktop/depict2/plots/v51_correlation_heatmaps.pdf")
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

