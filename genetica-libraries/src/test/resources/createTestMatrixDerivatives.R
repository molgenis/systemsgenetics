
setwd("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\genetica-libraries\\src\\test\\resources\\")

testMatrix <- read.delim("testMatrix.txt", row.names = 1)

write.table(cor(testMatrix), file = "testMatrixColumnCorMatrix.txt", quote = F, sep = "\t", col.names = NA)
write.table(cov(testMatrix), file = "testMatrixColumnCovMatrix.txt", quote = F, sep = "\t", col.names = NA)

write.table(scale(testMatrix), file = "testMatrixColumnScaledMatrix.txt", quote = F, sep = "\t", col.names = NA)



