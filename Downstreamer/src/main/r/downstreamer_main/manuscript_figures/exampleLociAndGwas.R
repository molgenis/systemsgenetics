
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

library(readr)
cedGwas <- as.data.frame(read_delim("cedGwas.txt", delim = "\t", quote = ""))
str(cedGwas)
cedGwas$logP <- -log10(cedGwas$pValue)


cedDsGenes <- read.delim("cedGenes.txt", quote = "", stringsAsFactors = F)
str(cedDsGenes)

cedDsGenes$GWAS.gene.P.valueLog <- -log10(cedDsGenes$GWAS.gene.P.value)

#Rel locus2:60686829-61686829

maxLog10GeneP <- 16
chr=2
start=60686829
stop=61686829


cedGwasNfkb1 <- cedGwas[cedGwas$Chr==2 & cedGwas$Pos >= 60686829 & cedGwas$Pos <= 61686829 ,]

bonfGwasLogP <- -log10(0.05 / nrow(cedGwas))
bonfZscore <- 4.71


str(cedGwasNfkb1)

topPos <- cedGwasNfkb1$Pos[which.max(cedGwasNfkb1$logP)]

layout(matrix(1:3, ncol = 1))

plot.new()
plot.window(xlim = c(start,stop), ylim = c(0,maxLog10GeneP))
axis(1)
axis(2)
abline(v=topPos, col = "firebrick2", lwd = 2)
points(cedGwasNfkb1$Pos, cedGwasNfkb1$logP, bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))

abline(h=-log10(5e-8), col = "grey80", lwd = 1, lty = 2)

cedDsGenesNfkb1 <- cedDsGenes[cedDsGenes$Chromosome.Name=="2" & cedDsGenes$Gene.Start..bp. <= 61686829 & cedDsGenes$Gene.End..bp. >= 60686829, ]
str(cedDsGenesNfkb1)

plot.new()
plot.window(xlim = c(start,stop), ylim = c(0,maxLog10GeneP))
axis(1)
axis(2)
for(i in 1:nrow(cedDsGenesNfkb1)){
  lines(c(cedDsGenesNfkb1$Gene.Start..bp.[i],cedDsGenesNfkb1$Gene.End..bp.[i]), y = c(cedDsGenesNfkb1$GWAS.gene.P.valueLog[i],cedDsGenesNfkb1$GWAS.gene.P.valueLog[i]), lwd = 4)
}
abline(v=topPos, col = "firebrick2", lwd = 2)
abline(h=bonfGwasLogP, col = "grey80", lwd = 1, lty = 2)


maxDsZscore <-  max(cedDsGenesNfkb1$Enrichment.Z.score,na.rm = T)
minDsZscore <-  min(cedDsGenesNfkb1$Enrichment.Z.score,na.rm = T)




plot.new()
plot.window(xlim = c(start,stop), ylim = c(minDsZscore,maxDsZscore))
axis(1)
axis(2)

for(i in 1:nrow(cedDsGenesNfkb1)){
  lines(c(cedDsGenesNfkb1$Gene.Start..bp.[i],cedDsGenesNfkb1$Gene.End..bp.[i]), y = c(cedDsGenesNfkb1$Enrichment.Z.score[i],cedDsGenesNfkb1$Enrichment.Z.score[i]), lwd = 4)
}
abline(v=topPos, col = "firebrick2", lwd = 2)
abline(h=bonfZscore, col = "grey80", lwd = 1, lty = 2)








chr=4
start=103422486 - 250000
stop=103538459 + 250000




cedGwasNfkb1 <- cedGwas[cedGwas$Chr==chr & cedGwas$Pos >= start & cedGwas$Pos <= stop ,]


str(cedGwasNfkb1)

topPos <- cedGwasNfkb1$Pos[which.max(cedGwasNfkb1$logP)]

layout(matrix(1:3, ncol = 1))

plot.new()
plot.window(xlim = c(start,stop), ylim = c(0,maxLog10GeneP))
axis(1)
axis(2)
abline(v=topPos, col = "firebrick2", lwd = 2)
points(cedGwasNfkb1$Pos, cedGwasNfkb1$logP, bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))

abline(h=-log10(5e-8), col = "grey80", lwd = 1, lty = 2)

cedDsGenesNfkb1 <- cedDsGenes[cedDsGenes$Chromosome.Name==chr & cedDsGenes$Gene.Start..bp. <= stop & cedDsGenes$Gene.End..bp. >= start, ]

cedDsGenesNfkb1$GWAS.gene.P.valueLog[is.na(cedDsGenesNfkb1$GWAS.gene.P.valueLog)]<-0

str(cedDsGenesNfkb1)

plot.new()
plot.window(xlim = c(start,stop), ylim = c(0,maxLog10GeneP))
axis(1)
axis(2)
for(i in 1:nrow(cedDsGenesNfkb1)){
  lines(c(cedDsGenesNfkb1$Gene.Start..bp.[i],cedDsGenesNfkb1$Gene.End..bp.[i]), y = c(cedDsGenesNfkb1$GWAS.gene.P.valueLog[i],cedDsGenesNfkb1$GWAS.gene.P.valueLog[i]), lwd = 4)
}
abline(v=topPos, col = "firebrick2", lwd = 2)
abline(h=bonfGwasLogP, col = "grey80", lwd = 1, lty = 2)


maxDsZscore <-  max(cedDsGenesNfkb1$Enrichment.Z.score,na.rm = T)
minDsZscore <-  min(cedDsGenesNfkb1$Enrichment.Z.score,na.rm = T)


plot.new()
plot.window(xlim = c(start,stop), ylim = c(minDsZscore,maxDsZscore))
axis(1)
axis(2)

for(i in 1:nrow(cedDsGenesNfkb1)){
  lines(c(cedDsGenesNfkb1$Gene.Start..bp.[i],cedDsGenesNfkb1$Gene.End..bp.[i]), y = c(cedDsGenesNfkb1$Enrichment.Z.score[i],cedDsGenesNfkb1$Enrichment.Z.score[i]), lwd = 4)
}
abline(v=topPos, col = "firebrick2", lwd = 2)
abline(h=bonfZscore, col = "grey80", lwd = 1, lty = 2)





plot.new()
plot(cedGwas$logP, bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
