
setwd("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\genetica-libraries\\src\\test\\resources\\")


write.table(cor(finalx), file = "testMatrixColumnCorMatrix.txt", quote = F, sep = "\t", col.names = NA)
write.table(cov(finalx), file = "testMatrixColumnCovMatrix.txt", quote = F, sep = "\t", col.names = NA)



