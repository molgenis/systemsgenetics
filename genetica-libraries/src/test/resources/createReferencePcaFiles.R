A <- matrix(c(1,2,3,4,5,6), ncol = 2)


pca <- prcomp(A)
str(pca)
View(pca)

C <- t(t(A) - apply(t(A), 1, mean))


V <- cov(C) 


eigenRes <- eigen(V)

t(t(eigenRes$vectors) %*% t(C))
pca$x



setwd("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\genetica-libraries\\src\\test\\resources\\")

testMatrix <- read.delim("testMatrix.txt", row.names = 1)

pca <- prcomp(testMatrix)


C <- t(t(testMatrix) - apply(t(testMatrix), 1, mean))


V <- cov(C) 

eigenRes <- eigen(V)
