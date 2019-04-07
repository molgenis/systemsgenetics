
setwd("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\genetica-libraries\\src\\test\\resources\\")

testMatrix <- read.delim("testMatrix.txt", row.names = 1)

write.table(cor(testMatrix), file = "testMatrixColumnCorMatrix.txt", quote = F, sep = "\t", col.names = NA)
write.table(cov(testMatrix), file = "testMatrixColumnCovMatrix.txt", quote = F, sep = "\t", col.names = NA)

write.table(scale(testMatrix), file = "testMatrixColumnScaledMatrix.txt", quote = F, sep = "\t", col.names = NA)



library("psych")

corPvalues <- corr.test(testMatrix, adjust = "none")$p

corZscores <- apply(corPvalues, 1:2, function(p){qnorm(1 - (p/2))})

write.table(corZscores, file = "testMatrixColumnCorZscoreMatrix.txt", quote = F, sep = "\t", col.names = NA)
