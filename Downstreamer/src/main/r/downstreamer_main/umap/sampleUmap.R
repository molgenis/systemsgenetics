#remoter::server(verbose = T, port = 55557, password = "laberkak", sync = T)

remoter::client("localhost", port = 55557, password = "laberkak")

setwd("/groups/umcg-wijmenga/tmp04/projects/depict2/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

library(readr)
library(gatepoints)

#if(FALSE){
#  #don't run by accident
#  table_tmp <- read_delim("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.txt.gz", delim = "\t", quote = "")
#  gnExp <- as.matrix(table_tmp[,-1])
#  rownames(gnExp) <- table_tmp[,1][[1]]
#  rm(table_tmp)
#  saveRDS(gnExp, "/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.rds")
#}
#gnExp <- readRDS("/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle/reference_datasets/human_b37/expression_databases/31.07.pc1.illumina.wantedgenes.DESeqNorm_log2.duplicatesRemoved.ProbesWithZeroVarianceRemoved.CovariatesRemoved.rds")

ensg75pc <- read.delim("ensgR75_protein_coding.txt", stringsAsFactors = F)

table_tmp <- read_delim("/groups/umcg-wijmenga/prm02/data_projects/Gado/GeneNetwork_V2_01-02-2018/Covariates/PCA/pc-scores250.txt.gz", delim = "\t", quote = "")
pcs <- as.matrix(table_tmp[,-1])
rownames(pcs) <- table_tmp[,1][[1]]
rm(table_tmp)
pcs <- pcs[,1:165]
str(pcs)

library(uwot)

init <- pcs[,c(2,1)]
init <- init * -1



 
#sampleUmap <- umap(pcs, n_threads = 24, n_epochs = 1000, init = init,  n_neighbors = 500, min_dist = 2, init_sdev = 1e-4, learning_rate = 1, spread = 40 ,scale = "none", nn_method = "fnn")
#sampleUmap <- umap(pcs, n_threads = 24, n_epochs = 10000, init = init,  n_neighbors = 500, min_dist = 2, init_sdev = 1e-4, learning_rate = 2, spread = 40 ,scale = "none", nn_method = "fnn")


#rownames(sampleUmap) <- rownames(pcs)
#colnames(sampleUmap) <- colnames(c("UMAP1", "UMAP2"))

#write.table(sampleUmap, file = "Umap/umap2.txt", sep = "\t", quote = F)


sampleUmap <- read.delim("Umap/umap2.txt")



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



meanX <- sapply(tissueColMap$PlotClass, function(plotClass){
  mean(sampleUmap[sampleAnnotation$PlotClass == plotClass,1])
})

meanY <- sapply(tissueColMap$PlotClass, function(plotClass){
  mean(sampleUmap[sampleAnnotation$PlotClass == plotClass,2])
})

#clusterLabels <- data.frame("PlotClass" = tissueColMap$PlotClass, "centerX" = meanX, "centerY" = meanY, "offsetX" = "", "offsetY" = "", "label" = tissueColMap$PlotClass)
#write.table(clusterLabels, "Umap/lables.txt", sep = "\t", quote = F, row.names = F)


str(clusterLabels)


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
sampleAnnotation$plotOrder[sampleAnnotation$col != defaultCol] <- 2

plotOrder <- order(sampleAnnotation$plotOrder)

sampleUmapPlot <- sampleUmap[plotOrder,]
sampleCol <- sampleAnnotation$col[plotOrder]

#all(row.names(sampleAnnotation) == row.names(sampleUmap))

str(sampleAnnotation)

#X11()
#rpng( width = 1200, height = 1000)

pdf("Umap/sampleUmap.pdf", width = 14, height = 7.5, useDingbats = F, title = "GeneNetwork UMAP")


createUmap(sampleUmapPlot, sampleCol)

dev.off()











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







#lines()  

#lines <- locator(20)

#lines2 <- data.frame(
#x1 = lines$x[seq(1, to = length(lines$x), by = 2)],
#x2 = lines$x[seq(2, to = length(lines$x), by = 2)] ,
#y1 = lines$y[seq(1, to = length(lines$x), by = 2)],
#y2 = lines$y[seq(2, to = length(lines$x), by = 2)]
#)

#zoomLines <- rbind(zoomLines, zoomLines2)

#axis(1)
  #axis(2)

write.table(lines2, "Umap/lines.txt", sep = "\t", quote = F, row.names = F)




fhs(sampleUmap)

library(grid)
grid.locator(unit = "npc")
  
abline(h=-70)
abline(h=-100)
abline(v=-40)
abline(v=-60)

library(RColorBrewer)
display.brewer.all()
cat(RColorBrewer::brewer.pal(12, "Set3"), sep = "\n")



colfunc <- colorRampPalette(c("khaki2", "orangered", "brown4"))


cols = c("orangered1", "orchid1", "magenta", "indianred1", "hotpink2", "darkorange1", "darkorchid3", "deeppink4", "brown3", "darkred", "khaki2", "goldenrod1", "orange", "mediumorchid4")
length(cols)
plot(rep(1,14),col=cols,pch=19,cex=3)


cols = c("", "", "", "", "", "", "", "", "", "", "", "", "")


source(paste0("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\Downstreamer\\src\\main\\r\\downstreamer_main/downstreamer_functions.r"))
ced <- read.depict2("parkinsons_nalls2019_harm_jan_enrichtments_exHla.xlsx")
name <- "Parkinsons"

#ced <- read.depict2("multiple_sclerosis_patsopoulos_harm_jan_enrichtments_exHla_1.xlsx")
str(ced)

cedExp <- ced$Expression
cedExp <- cedExp[match(row.names(sampleUmap), cedExp$Sample),]

all(cedExp$Sample == row.names(sampleUmap))

str(cedExp)
#colCed <- brewer.pal(9, "GnBu")[as.numeric(cut(cedExp$Enrichment.Z.score[match(row.names(sampleUmap), cedExp$Sample)],breaks = 9))]

colCed <- rep("grey90", nrow(sampleUmap)) #lightblue1
#colCed[cedExp$Enrichment.Z.score > 0 & cedExp$FDR.5..significant] <- "orange"
#colCed[cedExp$Enrichment.Z.score > 0 & cedExp$Bonferroni.significant] <- "red4"


fdrZthreshold <- min(abs(cedExp$Enrichment.Z.score[cedExp$FDR.5..significant]))
bonfZscore <- min(abs(cedExp$Enrichment.Z.score[cedExp$Bonferroni.significant]))
maxZ <- max(cedExp$Enrichment.Z.score)
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


pdf(paste0("Umap/sampleUmap", name , ".pdf"), width = 14, height = 7.5, useDingbats = F, title = paste0("GeneNetwork UMAP, ", name, " enrichments"))

createUmap(sampleUmaporder, sampleCol)

posLegend = 0.86

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

mtext(paste0("Sample enrichment for: ", name ), side = 3, line = -4.2, col = "gray30")

dev.off()



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

