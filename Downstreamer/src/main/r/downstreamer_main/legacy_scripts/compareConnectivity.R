setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")

connectivity <- read.delim("./connectivity.txt")

str(connectivity)


gnomadDs <- read.delim("./max_coregulation_zscore_and_gnomad_metrics.tsv", stringsAsFactors = F)
str(gnomadDs)

connectivityGnomad <- merge(connectivity, gnomadDs, by.x = "X.", by.y = "gene")
str(connectivityGnomad)

layout(matrix(1:6, byrow = T, ncol =2))
plot(connectivityGnomad$sumChi2, connectivityGnomad$max_ds_zscore, xlab = "sumChi2 co-regulation z-scores", ylab = "Max DS z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X20.0, connectivityGnomad$max_ds_zscore, xlab = "Count co-regulation z-scores >= 20", ylab = "Max DS z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X10.0, connectivityGnomad$max_ds_zscore, xlab = "Count co-regulation z-scores >= 10", ylab = "Max DS z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X5.0, connectivityGnomad$max_ds_zscore, xlab = "Count co-regulation z-scores >= 5", ylab = "Max DS z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X2.0, connectivityGnomad$max_ds_zscore, xlab = "Count co-regulation z-scores >= 2", ylab = "Max DS z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X1.0, connectivityGnomad$max_ds_zscore, xlab = "Count co-regulation z-scores >= 1", ylab = "Max DS z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))

layout(matrix(1:6, byrow = T, ncol =2))
plot(connectivityGnomad$sumChi2, connectivityGnomad$lof_z, xlab = "sumChi2 co-regulation z-scores", ylab = "LOF Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X20.0, connectivityGnomad$lof_z, xlab = "Count co-regulation z-scores >= 20", ylab = "LOF Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X10.0, connectivityGnomad$lof_z, xlab = "Count co-regulation z-scores >= 10", ylab = "LOF Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X5.0, connectivityGnomad$lof_z, xlab = "Count co-regulation z-scores >= 5", ylab = "LOF Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X2.0, connectivityGnomad$lof_z, xlab = "Count co-regulation z-scores >= 2", ylab = "LOF Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X1.0, connectivityGnomad$lof_z, xlab = "Count co-regulation z-scores >= 1", ylab = "LOF Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))

layout(matrix(1:6, byrow = T, ncol =2))
plot(connectivityGnomad$sumChi2, connectivityGnomad$mis_z, xlab = "sumChi2 co-regulation z-scores", ylab = "Missense Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X20.0, connectivityGnomad$mis_z, xlab = "Count co-regulation z-scores >= 20", ylab = "Missense Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X10.0, connectivityGnomad$mis_z, xlab = "Count co-regulation z-scores >= 10", ylab = "Missense Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X5.0, connectivityGnomad$mis_z, xlab = "Count co-regulation z-scores >= 5", ylab = "Missense Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X2.0, connectivityGnomad$mis_z, xlab = "Count co-regulation z-scores >= 2", ylab = "Missense Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))
plot(connectivityGnomad$X1.0, connectivityGnomad$mis_z, xlab = "Count co-regulation z-scores >= 1", ylab = "Missense Z-score", bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.5))


library(hexbin)
bin<-hexbin(connectivityGnomad$sumChi2, connectivityGnomad$max_ds_zscore, xbins=50)
plot(bin, main="") 

cor.test(connectivityGnomad$sumChi2, connectivityGnomad$lof_z)
cor.test(connectivityGnomad$max_ds_zscore, connectivityGnomad$lof_z)

cor.test(connectivityGnomad$sumChi2, connectivityGnomad$mis_z)
  

write.table(connectivityGnomad, file = "connectivityGnomadMaxDsZ.txt", sep = "\t", quote = F, row.names = F)
