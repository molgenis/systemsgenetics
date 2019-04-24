
setwd("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\genetica-libraries\\src\\test\\resources\\")

testMatrix <- read.delim("testMatrix.txt", row.names = 1)

write.table(cor(testMatrix), file = "testMatrixColumnCorMatrix.txt", quote = F, sep = "\t", col.names = NA)
write.table(cov(testMatrix), file = "testMatrixColumnCovMatrix.txt", quote = F, sep = "\t", col.names = NA)

write.table(scale(testMatrix), file = "testMatrixColumnScaledMatrix.txt", quote = F, sep = "\t", col.names = NA)



library("psych")

corPvalues <- corr.test(testMatrix, adjust = "none")$p

corZscores <- apply(corPvalues, 1:2, function(p){qnorm(1 - (p/2))})

corDirections <- apply(cor(testMatrix), 1:2, function(x){if(x<0){return(-1)}else{return(1)}})

corZscoresDirected <- corZscores * corDirections

write.table(corZscoresDirected, file = "testMatrixColumnCorZscoreMatrix.txt", quote = F, sep = "\t", col.names = NA)


library(weights)

weights <- runif (100)
weights = weights + (1 - mean(weights))

weightedCorRes <- wtd.cor(testMatrix, weight = weights, mean1 = FALSE)
weightedCor <- weightedCorRes$correlation

write.table(weights, file = "randomWeights.txt", quote = F, sep = "\t", col.names = NA)
write.table(weightedCor, file = "testMatrixColumnWeightedCorMatrix.txt", quote = F, sep = "\t", col.names = NA)
