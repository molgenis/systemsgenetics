smartSeq <- read.delim("smartseqSamples.txt")[,1]
str(smartSeq)


defaultCol <- adjustcolor("grey", alpha.f = 0.3)
pcsAndMeta$colSmartseq <- defaultCol
pcsAndMeta$colSmartseq[pcsAndMeta$Row.names %in% smartSeq] <- "darkslategray"
pcsAndMeta$colSmartseq[pcsAndMeta$Row.names %in% qseq] <- "orangered"
pcsAndMeta$colSmartseq[outliersPc1 == "TRUE" ] <- "darkblue"
plotOrderSmartseq <- order((pcsAndMeta$colSmartseq != defaultCol) + 1)

rpng()
plot(pcsAndMeta[plotOrderSmartseq,"PC_1"], pcsAndMeta[plotOrderSmartseq,"PC_2"], col = pcsAndMeta$colSmartseq[plotOrderSmartseq], cex = 0.4, main = "Smartseq")
dev.off()

rpng()
plot(pcsAndMeta[plotOrderSmartseq,"PC_1"], pcsAndMeta[plotOrderSmartseq,"PC_6"], col = pcsAndMeta$colSmartseq[plotOrderSmartseq], cex = 0.4, main = "Quantseq and Smartseq", pch =16)

library(pROC)

smartSeqClass <- as.factor(pcsAndMeta$Row.names %in% smartSeq)
table(smartSeqClass)
dim(pcsAndMeta)
smartSeqAuc <- apply(pcsAndMeta[,2:100],2,function(x){
  tryCatch(
    {
      #wilcox.test(x ~ smartSeqClass)$p.value
      as.numeric(auc(response = smartSeqClass, predictor = x))
    },
    error=function(cond){return(1)}
  )
})
sort(smartSeqAuc)
str(pcsAndMeta[,2])
boxplot(pcsAndMeta[,2]~smartSeqClass)

boxplot(log10(pcsAndMeta[,"sra.sample_spots"]) ~ smartSeqClass )

sum(pcsAndMeta[,"sra.sample_spots"] < 10000000, na.rm = T)


library(vioplot)
vioplot(log10(pcsAndMeta[,"sra.sample_spots"]) ~ smartSeqClass)


plot(pcsAndMeta[plotOrderSmartseq,"PC_6"], log10(pcsAndMeta[plotOrderSmartseq,"recount_qc.star.number_of_splices:_total"]), col = pcsAndMeta$colSmartseq[plotOrderSmartseq], cex = 0.6)




defaultCol <- adjustcolor("grey", alpha.f = 0.3)
sum(!is.na(pcsAndMeta[,"recount_qc.star.number_of_splices:_total"]) & pcsAndMeta[,"recount_qc.star.number_of_splices:_total"] < 150000)

pcsAndMeta$colSmartseq <- defaultCol
pcsAndMeta$colSmartseq[!is.na(pcsAndMeta[,"recount_qc.star.number_of_splices:_total"]) & pcsAndMeta[,"recount_qc.star.number_of_splices:_total"] < 150000] <- "aquamarine2"
table(pcsAndMeta$colSmartseq)
plotOrderSmartseq <- order((pcsAndMeta$colSmartseq != defaultCol) + 1)

plot(pcsAndMeta[plotOrderSmartseq,"PC_6"], pcsAndMeta[plotOrderSmartseq,"PC_1"], col = pcsAndMeta$colSmartseq[plotOrderSmartseq], cex = 0.4, main = "recount_qc.star.number_of_splices:_total < 150000")



defaultCol <- adjustcolor("grey", alpha.f = 0.5)
pcsAndMeta$colSmartseq <- defaultCol
pcsAndMeta$colSmartseq[!is.na(pcsAndMeta[,"recount_seq_qc.%a"]) & pcsAndMeta[,"recount_seq_qc.%a"] < 20] <- "darkslateblue"
table(pcsAndMeta$colSmartseq)
plotOrderSmartseq <- order((pcsAndMeta$colSmartseq != defaultCol) + 1)

plot(pcsAndMeta[plotOrderSmartseq,"PC_6"], pcsAndMeta[plotOrderSmartseq,"PC_1"], col = pcsAndMeta$colSmartseq[plotOrderSmartseq], cex = 0.6, main = "recount_seq_qc.%c < 20")






pcNoCenter <- read.delim("Components.txt", sep = ",", row.names = 1)
pcNoCenter <- merge(pcNoCenter, combinedMeta, all.x = T, by = 0)
rownames(pcNoCenter) <- pcNoCenter$Row.names

defaultCol <- adjustcolor("grey", alpha.f = 0.3)
pcNoCenter$colSmartseq <- defaultCol
pcNoCenter$colSmartseq[rownames(pcNoCenter) %in% smartSeq] <- "darkslategray"
table(pcNoCenter$colSmartseq)
plotOrderSmartseq <- order((pcNoCenter$colSmartseq != defaultCol) + 1)

plot(pcNoCenter[plotOrderSmartseq,"X0"], pcNoCenter[plotOrderSmartseq,"X1"], col = pcNoCenter$colSmartseq[plotOrderSmartseq], cex = 0.4, main = "Smartseq")


plot(pcNoCenter[plotOrderSmartseq,"X0"], pcNoCenter[plotOrderSmartseq,"X1"], col = adjustcolor("grey", alpha.f = 0.2), cex = 0.3)









defaultCol <- adjustcolor("grey", alpha.f = 0.6)
pcNoCenter$col <- defaultCol

tissueAndCol <- tolower(pcNoCenter[,"Tissue"]) %in% tolower(tissueCol$PlotClass)

pcNoCenter$col[tissueAndCol] <- tissueCol$col[match(tolower(pcNoCenter[tissueAndCol,"Tissue"]), tolower(tissueCol$PlotClass))]


tissue2AndCol <- tolower(pcNoCenter[,"Tissue2"]) %in% tolower(tissueCol$PlotClass)
sum(tissue2AndCol)
pcNoCenter$col[tissue2AndCol] <- tissueCol$col[match(tolower(pcNoCenter[tissue2AndCol,"Tissue2"]), tolower(tissueCol$PlotClass))]



sum(is.na(tolower(pcNoCenter[,"Tissue"]) %in% tolower(tissueCol$PlotClass)))

#pcNoCenter$col <- tissueCol$col[match(tolower(pcNoCenter[,"Tissue"]), tolower(tissueCol$PlotClass), nomatch = nrow(tissueCol))]

plotOrder <- order((pcNoCenter$col != defaultCol) + 1)

plot(pcNoCenter[plotOrder,"X0"], pcNoCenter[plotOrder,"X1"], col = pcNoCenter$col[plotOrder], cex = 0.4)


pcNoCenter$gtexCol <- defaultCol
pcNoCenter$gtexCol[pcNoCenter$Cohort == "GTEx" ] <- "goldenrod3"
pcNoCenter$gtexCol[pcNoCenter$Cohort == "TCGA" ] <- "cyan1"

plotOrder <- order((pcNoCenter$gtexCol != defaultCol) + 1)
plot(pcNoCenter[plotOrder,"X0"], pcNoCenter[plotOrder,"X1"], col = pcNoCenter$gtexCol[plotOrder], cex = 0.4)




toExclude <- 
  (!is.na(pcsAndMeta[,"recount_qc.star.number_of_splices:_total"]) & pcsAndMeta[,"recount_qc.star.number_of_splices:_total"] < 15000) |
  (!is.na(pcsAndMeta[,"recount_seq_qc.%a"]) & pcsAndMeta[,"recount_seq_qc.%a"] < 20) |
  (!is.na(pcsAndMeta[,"recount_seq_qc.%t"]) & pcsAndMeta[,"recount_seq_qc.%t"] < 20) |
  (!is.na(pcsAndMeta[,"recount_seq_qc.%c"]) & pcsAndMeta[,"recount_seq_qc.%c"] < 20) |
  (!is.na(pcsAndMeta[,"recount_seq_qc.%g"]) & pcsAndMeta[,"recount_seq_qc.%g"] < 20)
sum(toExclude)


samplesToKeep <- pcsAndMeta$Row.names[!toExclude]
length(samplesToKeep) + sum(toExclude) == nrow(pcs) 

write.table(samplesToKeep, file = "samplesToKeep.txt", row.names = F, quote = F)


boxplot(pcsAndMeta[,"PC_1"])

outliersPc1 <- as.factor(pcsAndMeta[,"PC_2"] >= 120)
table(outliersPc1)
library(pROC)
outlierAuc <- sapply(colnames(pcsAndMeta),function(x){
    tryCatch(
      {
        #wilcox.test(x ~ smartSeqClass)$p.value
        as.numeric(auc(response = outliersPc1, predictor = pcsAndMeta[,x]))
      },
      error=function(cond){return(NA)}
    )
  })
sort(outlierAuc)





auc(response = outliersPc1, predictor = pcsAndMeta[,"PC_2"])
