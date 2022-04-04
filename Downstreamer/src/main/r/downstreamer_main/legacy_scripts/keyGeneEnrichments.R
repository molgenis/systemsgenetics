
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

source(paste0("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\Downstreamer\\src\\main\\r\\downstreamer_main/downstreamer_functions.r"))


traits <- read.delim("traits.txt", stringsAsFactors = F)

pdf("Enrichment_comparisons.pdf", width = 7, height = 7, useDingbats = F)

par(mar= c(5,5,4,4), xpd = NA)

for(i in 1:nrow(traits)){
  
  fileNameGwasEnrich <- ""
  fileNamePrioEnrich <- ""
  if(traits[i,"pheno"]==""){
    fileNameGwasEnrich = paste0(traits[i,"machine_friendly_id"],"_GenePvalues_Enrichment.xlsx")
    fileNamePrioEnrich = paste0(traits[i,"machine_friendly_id"],"_GenePrioritization_Enrichment.xlsx")
  } else {
    fileNameGwasEnrich = paste0(traits[i,"machine_friendly_id"],"_GenePvalues_Enrichment_", traits[i,"pheno"] , ".xlsx")
    fileNamePrioEnrich = paste0(traits[i,"machine_friendly_id"],"_GenePrioritization_Enrichment_", traits[i,"pheno"] , ".xlsx")
  }
  
  gwasEnrichFile = paste0("D:\\UMCG\\FrankeSwertz - Documents\\Projects\\downstreamer\\Results\\enrichments\\",fileNameGwasEnrich)
  prioEnrichFile = paste0("D:\\UMCG\\FrankeSwertz - Documents\\Projects\\downstreamer\\Results\\enrichments\\",fileNamePrioEnrich)
  name = traits[i,"Name"]
  
  gwasEnrich <- read.depict2(gwasEnrichFile)
  prioEnrich <- read.depict2(prioEnrichFile)
  
  names(gwasEnrich)
  
  gwasEnrichHpo <- gwasEnrich$HPO
  prioEnrichHpo <- prioEnrich$HPO
  
  row.names(gwasEnrichHpo) <- gwasEnrichHpo$Human.phenotype
  row.names(prioEnrichHpo) <- prioEnrichHpo$Human.phenotype

  
  axisRange <- range(-log10(gwasEnrichHpo$P.bonf), -log10(prioEnrichHpo$P.bonf))
  #plot.new()
  #plot.window(xlim = c(-150,100), ylim = c(-200,-50))
  #plot.window(xlim = axisRange, ylim = axisRange, asp = 1)
  #axis(side = 1, col = "gray30", lwd = 1, col.axis = "gray30")
  #axis(side = 2, col = "gray30", lwd = 1, col.axis = "gray30")
  #lines(c(0,axisRange[2]),c(0,axisRange[2]), lwd = 2, col = "gray30")
  #title(main = paste0("Comparison of HPO-term enrichment of\nGWAS genes vs key-genes of ", name), xlab = "-log10(Enrichment GWAS genes)", ylab = "-log10(Enrichment key-genes genes)", col = "gray30")
  #points(-log10(gwasEnrichHpo$P.bonf), -log10(prioEnrichHpo[gwasEnrichHpo$Human.phenotype,"P.bonf"]), bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
  
  grid.arrange(grobs=list(p1, p2, p3), ncol=3, nrow=2)
  
  plotData <- data.frame(term = gwasEnrichHpo[,1], description = gwasEnrichHpo[,2], gwasEnrichLog10P = -log10(gwasEnrichHpo$P.bonf), keygeneEnrichLog10P = -log10(prioEnrichHpo[gwasEnrichHpo[,1],"P.bonf"]))
  
  
  labelRows <- unique(c(tail(order(plotData$gwasEnrichLog10P), n = 10), tail(order(plotData$keygeneEnrichLog10P), n = 10)))
  
  dim(plotData[, ])
  plotData[-labelRows, "description"] <- NA
  
  p1 <- ggplot(data=plotData, mapping=aes(x=gwasEnrichLog10P,
                                          y=keygeneEnrichLog10P,
                                          label=description)) +
    geom_point(alpha=0.2) +
    geom_text_repel(min.segment.length = 0, point.padding = 0.2, force = 10) +
    xlim(axisRange) +
    ylim(axisRange) +
    geom_abline() +
    xlab("piet") +
    ggtitle("henkie")
    
  
  theme.nature(p1) + theme(plot.title = element_text(size = 20, face = "bold"))
  
  d#
  
  
  
  
  
  
}

dev.off()
x <- merge(gwasEnrichHpo, prioEnrichHpo[,c("OR.bonf", "P.bonf")], by= 0)
View(x)
write.table(x, file = "tmp.txt", sep ="\t", quote  =F)
