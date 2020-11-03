source("../../downstreamer_main/downstreamer_functions.r")



path        <- "~/Desktop/collaborations/lude_cancer_gwas/"
files       <- list.files(path, pattern = "*.xlsx")
data        <- lapply(paste0(path, files), read.depict2)
names(data) <- files



plot.xy <- function(x.name, y.name, data) {
  p1 <- ggplot(data=df.plot, aes_string(x=x.name, y=y.name)) +
    geom_point(alpha=0.15) +
    geom_abline(slope=1, intercept=0, col="blue", lwd=1.25) +
    geom_smooth(method='lm', formula= y~x, col="red") +
    coord_fixed()
  return(theme.nature(p1))
}


type  <- "Coregulation"

order   <- rownames(data[["breast_cancer_triple_negative_2020_32424353_hg19_enrichtments_exHla.xlsx"]][[type]])
df.plot <- data.frame(BRC=data[["breast_cancer_triple_negative_2020_32424353_hg19_enrichtments_exHla.xlsx"]][[type]][order,]$Enrichment.Z.score,
                      CRC=data[["colorectal_cancer_2018_30104761_hg19_enrichtments_exHla.xlsx"]][[type]][order,]$Enrichment.Z.score,
                      EC=data[["endometrial_cancer_2018_30093612_hg19_enrichtments_exHla.xlsx"]][[type]][order,]$Enrichment.Z.score,
                      PC=data[["prostate_cancer_2018_29892016_hg19_enrichtments_exHla.xlsx"]][[type]][order,]$Enrichment.Z.score)#,
#                     kb150=kb150[[trait]][[type]][order,]$Enrichment.Z.score)

p1 <- plot.xy("BRC", "CRC", df.plot)
p2 <- plot.xy("kb50", "kb25", df.plot)
p3 <- plot.xy("kb50", "kb100", df.plot)
p4 <- plot.xy("kb50", "kb150", df.plot)


grid.arrange(grobs=list(p1, p2, p3, p4), ncol=4)
grid.arrange(grobs=list(p1, p2, p3), ncol=3)