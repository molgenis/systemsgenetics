setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

source(paste0("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\Downstreamer\\src\\main\\r\\downstreamer_main/downstreamer_functions.r"))

traits <- read.delim("MetaBrain/traits.txt")


i <- 1

pdf("MetaBrain/withAndWithoutEqtls.pdf", height = 20, width = 10)
#png("MetaBrain/withAndWithoutEqtls.png", height = 2000, width = 1000)
layout(matrix(1:8, ncol =2))
par(pty="s")
for(i in 1:nrow(traits)){

  

trait <- traits[i, "trait"]
name <- traits[i, "name"]

enrichments <- read.depict2(paste0("MetaBrain/normal/",trait,"_enrichtments.xlsx"))$GenePrioritization_MetaBrain
enrichmentsIncEqtl <- read.depict2(paste0("MetaBrain/inceqt/",trait,"_enrichtments.xlsx"))$GenePrioritization_MetaBrain

enrichmentsBoth <- merge(enrichments, enrichmentsIncEqtl,  "Gene.ID" , suffixes= c("Normal", "incEqtl"))

maxZ <- max(range(enrichmentsBoth$Enrichment.Z.scoreNormal, enrichmentsBoth$Enrichment.Z.scoreincEqtl))
r <- cor(enrichmentsBoth$Enrichment.Z.scoreNormal, enrichmentsBoth$Enrichment.Z.scoreincEqtl)
plot(enrichmentsBoth$Enrichment.Z.scoreNormal, enrichmentsBoth$Enrichment.Z.scoreincEqtl, bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5), asp = 1, xlab = "Key gene score without eqtl information",  ylab = "Key gene score without eqtl information", xlim = c(-maxZ,maxZ), ylim = c(-maxZ,maxZ), main = name)
mtext(paste0("Pearson r: ", signif(r,2)))
}
dev.off()
