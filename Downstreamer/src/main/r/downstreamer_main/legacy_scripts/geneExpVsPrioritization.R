setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")
library(plotfunctions)
sampleUmap <- read.delim("Umap/umap2.txt")


ensg <- read.delim("ensgR75_protein_coding.txt", stringsAsFactors = F)
row.names(ensg) <- ensg$Ensembl.Gene.ID
str(ensg)

exp <- readRDS("31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.rds")
exp <- t(as.matrix(exp))
exp <- exp[match(row.names(sampleUmap), row.names(exp)),,drop=F]
str(exp)



tissueColMap <- read.delim("Umap/tissueCol5.txt", stringsAsFactors = F)
sampleAnnotation <- read.delim("Umap/sampleAnnotations.txt", stringsAsFactors = F, row.names = 1 )
str(sampleAnnotation)

all(row.names(sampleAnnotation)  %in% row.names(sampleUmap))

sampleAnnotation <- sampleAnnotation[match(row.names(sampleUmap), rownames(sampleAnnotation) ),]

all(row.names(sampleAnnotation) == row.names(sampleUmap))


#sampleAnnotation$pch = 1
#sampleAnnotation$pch[sampleAnnotation$CellLine != ""] = 5
#sampleAnnotation$pch[sampleAnnotation$TissueType != ""] = 21

defaultCol <- adjustcolor("gray90", alpha.f = 0.1)
sampleAnnotation$col = defaultCol
#sampleAnnotation$col = NA
#sampleAnnotation$bg = NA



#meanX <- sapply(tissueColMap$PlotClass, function(plotClass){
#mean(sampleUmap[sampleAnnotation$PlotClass == plotClass,1])
#})

#meanY <- sapply(tissueColMap$PlotClass, function(plotClass){
#mean(sampleUmap[sampleAnnotation$PlotClass == plotClass,2])
#})

#clusterLabels <- data.frame("PlotClass" = tissueColMap$PlotClass, "centerX" = meanX, "centerY" = meanY, "offsetX" = "", "offsetY" = "", "label" = tissueColMap$PlotClass)
#write.table(clusterLabels, "Umap/lables.txt", sep = "\t", quote = F, row.names = F)



#sampleAnnotation$bg[sampleAnnotation$PlotClass == "osteoblastic cell"] <- adjustcolor("red", alpha.f = 0.5)


clusterLabels <- read.delim("Umap/lables.txt",stringsAsFactors = F, comment.char = "#")
clusterLabelsZoom <- read.delim("Umap/lablesZoom.txt",stringsAsFactors = F, comment.char = "#")

zoomLines <- read.delim("Umap/zoomLines.txt",stringsAsFactors = F, comment.char = "#")

lines <- read.delim("Umap/lines.txt",stringsAsFactors = F, comment.char = "#")


for(i in 1:nrow(tissueColMap)){
  sampleAnnotation$col[sampleAnnotation$PlotClass == tissueColMap$PlotClass[i]] <- adjustcolor(tissueColMap$col[i], alpha.f = 0.2)
  sampleAnnotation$bg[sampleAnnotation$PlotClass == tissueColMap$PlotClass[i]] <- adjustcolor(tissueColMap$col[i], alpha.f = 0.3)
}

#sampleAnnotation$col[sampleAnnotation$PlotClass == "skin"] <- adjustcolor("black", alpha.f = 0.5)

sampleAnnotation$plotOrder <- 1
sampleAnnotation$plotOrder[sampleAnnotation$colPredictedCell != defaultCol] <- 2

plotOrder <- order(sampleAnnotation$plotOrder)

sampleUmapPlot <- sampleUmap[plotOrder,]
sampleCol <- sampleAnnotation$col[plotOrder]



#X11()
#rpng( width = 1200, height = 1000)









genesToPlot <- c("ENSG00000167034", "ENSG00000142515", "ENSG00000165929", "ENSG00000066813", "ENSG00000137731", "ENSG00000198398", "ENSG00000117091", "ENSG00000145287", "ENSG00000090104", "ENSG00000166589", "ENSG00000156413")
#ENSG00000099985
#ENSG00000090104
#ENSG00000156413

for(gene in genesToPlot){

    gene <- "ENSG00000109182" 
  geneSymbol <- ensg[gene, "Associated.Gene.Name"]
  
      
  #pdf(paste0("expressionPlots/", gene, "_", geneSymbol , "_umap.pdf"), width = 14, height = 7.5, useDingbats = F, title = paste0("GeneNetwork UMAP, ", name, " enrichments"))

  geneExp <- exp[,gene]
  
  #plot(geneExp, cedExp$Enrichment.Z.score, xlab = paste0(geneSymbol, "expression"), ylab = "Sample enrichment Z-score for prostate cancer", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.2))
  
  cor.test(geneExp, cedExp$Enrichment.Z.score)
  
  
  breaks <- quantile(geneExp, probs = seq(0,1,0.01))
  
  all(names(geneExp) == cedExp$Sample)
  
  colfunc<-colorRampPalette(c("gold2","orange","orangered","red3"))
  colMapEnrich <- adjustcolor(colfunc(20), alpha.f = 0.5)
  colfunc<-colorRampPalette(c("lightblue1","dodgerblue4"))
  colMapDepl <- adjustcolor(colfunc(20), alpha.f = 0.5)
  
  fullColGradient <- c(rev(colMapDepl), rep(defaultCol, 20), colMapEnrich)
  
  sampleCol <- fullColGradient[cut(geneExp, 60)]
  
  sampleAnnotation$plotOrder <- 1
  sampleAnnotation$plotOrder[sampleAnnotation$colPredictedCell != defaultCol] <- 2
  
  plotOrder <- order(sampleAnnotation$plotOrder)
  
  sampleUmapPlot <- sampleUmap[plotOrder,]
  sampleCol <- sampleCol[plotOrder]
  
  
  createUmap(sampleUmapPlot, sampleCol)
  
  
  posLegend = 0.89
  gradientLegend(
    c(-min(geneExp),max(geneExp)),
    color = fullColGradient,
    pos = c(0.25,posLegend,0.75,posLegend + 0.02),
    side = 3,
    dec = 1,
    length = 0.5,
    depth = 0.025,
    inside = T,
    coords = FALSE,
    pos.num = NULL,
    n.seg = 1,
    border.col = NA,
    tick.col = "gray30",
    fit.margin = TRUE,
  )
  
  mtext(paste0("Expression of: ", geneSymbol ), side = 3, line = -2.9, col = "gray30")

  dev.off()
}

createUmap <- function(sampleUmapPlot, sampleCol){
  
  layout(matrix(1:2, nrow =1), widths = c(1.2,1))
  
  par(mar = c(3,3,0,0), xpd = NA)
  
  plot.new()
  #plot.window(xlim = c(-150,100), ylim = c(-200,-50))
  plot.window(xlim = c(-358,380), ylim = c(-402,335), asp = 1)
  axis(side = 1, at = c(-200,0,200), col = "gray30", lwd = 1, col.axis = "gray30")
  axis(side = 2, at = c(-200,0,200), col = "gray30", lwd = 1, col.axis = "gray30")
  points(sampleUmapPlot, col = sampleCol, pch = 16, cex = 0.6)
  #title(xlab = "UMAP-1", ylab = "UMAP-2", line = 2)
  mtext("UMAP-1",side = 1,at = 0, line = 2, col = "gray30")
  mtext("UMAP-2",side = 2,at = 0, line = 2, col = "gray30")
  
  apply(lines, 1, function(coord){
    lines(x = coord[1:2], y = coord[3:4], col = "gray70")
  })
  
  
  text(clusterLabels$centerX + clusterLabels$offsetX, clusterLabels$centerY + clusterLabels$offsetY, labels = clusterLabels$label, col = "gray30")
  
  zoomX <- c(-120,54)
  zoomY <- c(-202,-66)
  
  rect(zoomX[1], zoomY[1], zoomX[2], zoomY[2], border = "gray70")
  
  zoomLineStartX <- grconvertX(zoomX[2], from = "user", to = "device")
  zoomLineStartY1 <- grconvertY(zoomY[1], from = "user", to = "device")
  zoomLineStartY2 <- grconvertY(zoomY[2], from = "user", to = "device")
  
  zoomPoints <- sampleUmapPlot$V1 >= zoomX[1] & sampleUmapPlot$V1 <= zoomX[2] & sampleUmapPlot$V2 >= zoomY[1] & sampleUmapPlot$V2 <= zoomY[2]
  
  par(mar = c(0,0,0,0), xpd = NA)
  plot.new()
  plot.window(xlim = c(-115,50), ylim = c(-202,-65), asp = 1)
  #axis(side = 1, at = c(-200,0,200), col = "gray30", lwd = 1.5, col.axis = "gray30")
  #axis(side = 2, at = c(-200,0,200), col = "gray30", lwd = 1.5, col.axis = "gray30")
  points(sampleUmapPlot[zoomPoints, ], col = sampleCol[zoomPoints], pch = 16, cex = 0.9)
  
  lines(c(grconvertX(zoomLineStartX, from = "device"), zoomX[1]), c(grconvertY(zoomLineStartY1, from = "device"),zoomY[1]), col = "gray70")
  lines(c(grconvertX(zoomLineStartX, from = "device"), zoomX[1]), c(grconvertY(zoomLineStartY2, from = "device"),zoomY[2]), col = "gray70")
  rect(zoomX[1], zoomY[1], zoomX[2], zoomY[2], border = "gray70")
  
  
  text(clusterLabelsZoom$centerX + clusterLabelsZoom$offsetX, clusterLabelsZoom$centerY + clusterLabelsZoom$offsetY, labels = clusterLabelsZoom$label, col = "gray30")
  
  apply(zoomLines, 1, function(coord){
    lines(x = coord[1:2], y = coord[3:4], col = "gray70")
  })
  
  
}





source(paste0("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\Downstreamer\\src\\main\\r\\downstreamer_main/downstreamer_functions.r"))

traits <- read.delim("traits.txt", stringsAsFactors = F)

i <- which(traits$trait == "prostate_cancer")

for(i in 1:nrow(traits)){
  
  fileName <- ""
  if(traits[i,"pheno"]==""){
    fileName = paste0(traits[i,"machine_friendly_id"],"_enrichtments.xlsx")
  } else {
    fileName = paste0(traits[i,"machine_friendly_id"],"_enrichtments_", traits[i,"pheno"] , ".xlsx")
  }
  
  file = paste0("D:\\UMCG\\FrankeSwertz - Documents\\Projects\\downstreamer\\Results\\",fileName)
  name = traits[i,"Name"]
  
  ced <- read.depict2(file)
  
  #ced <- read.depict2("multiple_sclerosis_patsopoulos_harm_jan_enrichtments_exHla_1.xlsx")
  str(ced)
  
  cedExp <- ced$SampleEnrichment
  cedExp <- cedExp[match(row.names(sampleUmap), cedExp$Sample),]
  
  all(cedExp$Sample == row.names(sampleUmap))
  
  str(cedExp)
  #colCed <- brewer.pal(9, "GnBu")[as.numeric(cut(cedExp$Enrichment.Z.score[match(row.names(sampleUmap), cedExp$Sample)],breaks = 9))]
  
  colCed <- rep("grey90", nrow(sampleUmap)) #lightblue1
  #colCed[cedExp$Enrichment.Z.score > 0 & cedExp$FDR.5..significant] <- "orange"
  #colCed[cedExp$Enrichment.Z.score > 0 & cedExp$Bonferroni.significant] <- "red4"
  
  
  fdrZthreshold <- min(abs(cedExp$Enrichment.Z.score[cedExp$FDR.5..significant]))
  maxZ <- max(cedExp$Enrichment.Z.score)
  if(!is.infinite(fdrZthreshold)){
    
    bonfZscore <- min(abs(cedExp$Enrichment.Z.score[cedExp$Bonferroni.significant]))
    
    breaks <- seq(fdrZthreshold, maxZ, length.out = 21)
    
    greyBreakCount <- round((fdrZthreshold + fdrZthreshold ) / (breaks[2] - breaks[1]))
    
    colfunc<-colorRampPalette(c("gold2","orange","orangered","red3"))
    colMapEnrich <- colfunc(20)
    colCed[cedExp$Enrichment.Z.score > 0 & cedExp$FDR.5..significant] <- colMapEnrich[as.numeric(cut(cedExp$Enrichment.Z.score[cedExp$Enrichment.Z.score > 0 & cedExp$FDR.5..significant],breaks))]
    
    
    colfunc<-colorRampPalette(c("lightblue1","dodgerblue4"))
    colMapDepl <- colfunc(20)
    colCed[cedExp$Enrichment.Z.score < 0 & cedExp$FDR.5..significant] <- colMapDepl[as.numeric(cut(abs(cedExp$Enrichment.Z.score[cedExp$Enrichment.Z.score < 0 & cedExp$FDR.5..significant]),breaks))]
    
    fullColGradient <- c(rev(colMapDepl), rep("grey90", greyBreakCount), colMapEnrich)
    
    library( plotfunctions)
    
    sampleCol <- adjustcolor(colCed[order(abs(cedExp$Enrichment.Z.score),decreasing = F)], alpha.f = 0.5)
    
    sampleUmaporder <- sampleUmap[order(abs(cedExp$Enrichment.Z.score),decreasing = F),]
    
  }  else {
    fullColGradient <- rep("grey90", 52)
    sampleCol <- adjustcolor(colCed, alpha.f = 0.5)
    sampleUmaporder <- sampleUmap
    fdrZthreshold <- maxZ
    bonfZscore <- maxZ
  }
  
  
  pdf(paste0("Umap/sampleUmap", name , ".pdf"), width = 14, height = 7.5, useDingbats = F, title = paste0("GeneNetwork UMAP, ", name, " enrichments"))
  
  createUmap(sampleUmaporder, sampleCol)
  
  posLegend = 0.89
  
  gradientLegend2(
    c(-maxZ,maxZ),
    color = fullColGradient,
    pos = c(0.25,posLegend,0.75,posLegend + 0.02),
    side = 3,
    dec = 1,
    length = 0.5,
    depth = 0.025,
    inside = T,
    coords = FALSE,
    pos.num = NULL,
    n.seg = c(-maxZ,-bonfZscore,-fdrZthreshold,0,fdrZthreshold,bonfZscore,maxZ),
    border.col = NA,
    tick.col = "gray30",
    fit.margin = TRUE,
  )
  
  mtext(paste0("Sample enrichment for: ", name ), side = 3, line = -2.9, col = "gray30")
  
  dev.off()
}


#addapted from plotfunctions
gradientLegend2 <- function (valRange, color = "terrain", nCol = 30, pos = 0.875, 
                             side = 4, dec = NULL, length = 0.25, depth = 0.05, inside = FALSE, 
                             coords = FALSE, pos.num = NULL, n.seg = 1, border.col = "black", 
                             tick.col = NULL, fit.margin = TRUE, ...) 
{
  loc <- c(0, 0, 0, 0)
  ticks <- c()
  labels <- c()
  if (side %in% c(1, 3)) {
    if (length(pos) == 1) {
      pos.other <- ifelse(side > 2, 1, 0)
      switch <- ifelse(inside, 0, 1)
      switch <- ifelse(side > 2, 1 - switch, switch)
      loc <- getCoords(c(pos - 0.5 * length, pos.other - 
                           switch * depth, pos + 0.5 * length, pos.other + 
                           (1 - switch) * depth), side = c(side, 2, side, 
                                                           2))
    }
    else if (length(pos) == 4) {
      if (coords) {
        loc <- pos
      }
      else {
        loc <- getCoords(pos, side = c(1, 2, 1, 2))
      }
    }
    if (length(n.seg) == 1) {
      ticks <- seq(loc[1], loc[3], length = n.seg + 2)
      labels <- seq(min(valRange), max(valRange), length = n.seg + 
                      2)
    }
    else {
      labels <- c(min(valRange), sort(n.seg[n.seg > min(valRange) & 
                                              n.seg < max(valRange)], decreasing = FALSE), 
                  max(valRange))
      b <- diff(loc[c(1, 3)])/diff(range(labels))
      ticks <- (labels - min(labels)) * b + loc[1]
    }
  }
  else if (side %in% c(2, 4)) {
    if (length(pos) == 1) {
      pos.other <- ifelse(side > 2, 1, 0)
      switch <- ifelse(inside, 0, 1)
      switch <- ifelse(side > 2, 1 - switch, switch)
      loc <- getCoords(c(pos.other - switch * depth, pos - 
                           0.5 * length, pos.other + (1 - switch) * depth, 
                         pos + 0.5 * length), side = c(1, side, 1, side))
    }
    else if (length(pos) == 4) {
      if (coords) {
        loc <- pos
      }
      else {
        loc <- getCoords(pos, side = c(1, 2, 1, 2))
      }
    }
    if (length(n.seg) == 1) {
      ticks <- seq(loc[2], loc[4], length = n.seg + 2)
      labels <- seq(min(valRange), max(valRange), length = n.seg + 
                      2)
    }
    else {
      labels <- c(min(valRange), sort(n.seg[n.seg > min(valRange) & 
                                              n.seg < max(valRange)], decreasing = FALSE), 
                  max(valRange))
      b <- diff(loc[c(2, 4)])/diff(range(labels))
      ticks <- (labels - min(labels)) * b + loc[2]
    }
  }
  if (is.null(pos.num)) {
    if (side %in% c(1, 3)) {
      if (inside) {
        pos.num <- ifelse(side == 1, 3, 1)
      }
      else {
        pos.num <- side
      }
    }
    else {
      if (inside) {
        pos.num <- ifelse(side == 2, 4, 2)
      }
      else {
        pos.num <- side
      }
    }
  }
  getcol <- get_palette(color, nCol = nCol)
  mycolors <- getcol[["color"]]
  if (is.null(tick.col)) {
    tick.col = border.col
  }
  vals <- seq(min(valRange), max(valRange), length = length(mycolors))
  im <- as.raster(mycolors[matrix(1:length(mycolors), ncol = 1)])
  n <- max(c(length(ticks) - 2, 0))
  if (side%%2 == 1) {
    im <- t(im)
    rasterImage(im, loc[1], loc[2], loc[3], loc[4], col = mycolors, 
                xpd = T)
    segments(x0 = ticks, x1 = ticks, y0 = rep(loc[2]-1.5, n), 
             y1 = rep(loc[2], n), col = tick.col, xpd = TRUE)
    if(is.na(border.col) & !is.na(tick.col)){
      segments(loc[1], loc[2], loc[3], loc[2], col = tick.col)
    } 
    rect(loc[1], loc[2], loc[3], loc[4], border = border.col, 
         xpd = T)
    
    
  }
  else {
    im <- rev(im)
    rasterImage(im, loc[1], loc[2], loc[3], loc[4], col = mycolors, 
                xpd = T)
    segments(x0 = rep(loc[1], n), x1 = rep(loc[3], n), y0 = ticks, 
             y1 = ticks, col = tick.col, xpd = TRUE)
    rect(loc[1], loc[2], loc[3], loc[4], border = border.col, 
         xpd = T)
  }
  lab.loc.x <- lab.loc.y <- c()
  if (side %in% c(1, 3)) {
    lab.loc.x <- ticks
    if (pos.num == 1) {
      lab.loc.y <- rep(loc[2], length(lab.loc.x))
    }
    else if (pos.num == 3) {
      lab.loc.y <- rep(loc[4], length(lab.loc.x))
    }
    else {
      lab.loc.y <- rep((loc[2] + loc[4])/2, length(lab.loc.x))
    }
  }
  else if (side %in% c(2, 4)) {
    lab.loc.y <- ticks
    if (pos.num == 2) {
      lab.loc.x <- rep(loc[1], length(lab.loc.y))
    }
    else if (pos.num == 4) {
      lab.loc.x <- rep(loc[3], length(lab.loc.y))
    }
    else {
      lab.loc.x <- rep((loc[1] + loc[3])/2, length(lab.loc.y))
    }
  }
  determineDec <- function(x) {
    out = max(unlist(lapply(strsplit(as.character(x), split = "\\."), 
                            function(y) {
                              return(ifelse(length(y) > 1, nchar(gsub("^([^0]*)([0]+)$", 
                                                                      "\\1", as.character(y[2]))), 0))
                            })))
    return(out)
  }
  if (is.null(dec)) {
    dec <- min(c(6, determineDec(labels)))
  }
  eval(parse(text = sprintf("labels = sprintf('%s', round(labels, dec) )", 
                            paste("%.", dec, "f", sep = ""))))
  labels <- gsub("^(\\-)(0)([\\.0]*)$", "\\2\\3", labels)
  if (fit.margin == TRUE & inside == FALSE) {
    lab.height <- max(strheight(labels))
    lab.width <- max(strwidth(labels))
    lab.cor <- strheight("0") * 0.5
    max.pos <- getFigCoords("f")
    cex.f <- NA
    change <- NA
    if (pos.num == 1) {
      max.height = lab.loc.y[1] - max.pos[3]
      cex.f = max.height/(lab.height + lab.cor)
      change <- ifelse(cex.f < 0.8, (lab.height + lab.cor) * 
                         0.8 - max.height, NA)
    }
    else if (pos.num == 2) {
      max.width = lab.loc.x[1] - max.pos[1]
      cex.f = max.width/(lab.width + lab.cor)
      change <- ifelse(cex.f < 0.8, (lab.width + lab.cor) * 
                         0.8 - max.width, NA)
    }
    else if (pos.num == 3) {
      max.height = max.pos[4] - lab.loc.y[1]
      cex.f = max.height/(lab.height + lab.cor)
      change <- ifelse(cex.f < 0.8, (lab.height + lab.cor) * 
                         0.8 - max.height, NA)
    }
    else if (pos.num == 4) {
      max.width = max.pos[2] - lab.loc.x[1]
      cex.f = max.width/(lab.width + lab.cor)
      change <- ifelse(cex.f < 0.8, (lab.width + lab.cor) * 
                         0.8 - max.width, NA)
    }
    if (cex.f < 0.8) {
      margin <- c("bottom", "left", "top", "right")
      warning(sprintf("Increase %s margin to fit labels or decrease the number of decimals, see help(gradientLegend).", 
                      margin[pos.num]))
    }
    par <- list(...)
    if ("cex" %in% names(par)) {
      text(x = lab.loc.x, y = lab.loc.y, labels = labels, 
           col = tick.col, pos = pos.num, xpd = T, ...)
    }
    else {
      text(x = lab.loc.x, y = lab.loc.y, labels = labels, 
           col = tick.col, pos = pos.num, cex = min(c(0.8, 
                                                      cex.f)), xpd = T)
    }
  }
  else {
    par <- list(...)
    if ("cex" %in% names(par)) {
      text(x = lab.loc.x, y = lab.loc.y, labels = labels, 
           col = tick.col, pos = pos.num, xpd = T, ...)
    }
    else {
      text(x = lab.loc.x, y = lab.loc.y, labels = labels, 
           col = tick.col, pos = pos.num, cex = 0.8, xpd = T)
    }
  }
  invisible(list(loc = loc, ticks = ticks, labels = labels, 
                 im = im))
}

