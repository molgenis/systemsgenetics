

#
## Create geneZscoresNullGwas matrix used to find correlated genes genes x traits
#

#Code to generate matrix from: http://r.789695.n4.nabble.com/Re-Generate-a-serie-of-new-vars-that-correlate-with-existingvar-td822187.html

# generates ndistr vectors of same mean and sd, with various cor.coeffs
# input :
#         x1 : a vector
#         ndistr : number of distributions
#         coefs : vector o ndistr correl. coeffs

CorelSets<-function(x1= rnorm(100, 15, 5),ndistr=3, newc){
  
  # x2, x3, and x4 in a matrix, these will be modified to meet the criteria
  x234 <- scale(matrix( rnorm(ndistr*length(x1)), ncol=ndistr ))
  
  # put all into 1 matrix for simplicity
  x1234 <- cbind(scale(x1),x234)
  
  # find the current correlation matrix
  c1 <- var(x1234)
  
  # cholesky decomposition to get independence
  chol1 <- solve(chol(c1))
  
  newx <-  x1234 %*% chol1
  
  # check that we have independence and x1 unchanged
  zapsmall(cor(newx))
  all.equal( x1234[,1], newx[,1] )
  
  # create new correlation structure
  
  chol2 <- chol(newc)
  
  finalx <- newx %*% chol2 * sd(x1) + mean(x1)
  pairs(finalx)
  CorelSets<-finalx
} 


newc<-diag(8)
newc[1,2] <- 0.9
newc[1,3] <- 0.6
newc[1,4] <- 0.8
newc[2,3] <- 0.7
newc[2,4] <- 0.9
newc[3,4] <- 0.5

newc[5,6] <- 0.5
newc[5,7] <- 0.6
newc[6,7] <- 0.95

# newc[1,8] <- 0.95
# newc[2,8] <- 0.6
# newc[3,8] <- 0.3
# newc[4,8] <- 0.6
 newc[5,8] <- -0.91
 newc[6,8] <- -0.7
 newc[7,8] <- -0.8


newc[lower.tri(newc)] <- t(newc)[lower.tri(newc)]

chol(newc)

newc
heatmap(newc, scale= "none", Rowv = NA, Colv = NA)
heatmap((abs(newc) >= 0.9)+0, scale= "none", Rowv = NA, Colv = NA)

sum(abs(newc) >= 0.9)

geneZscoresNullGwas <- t(CorelSets(x1= rnorm(10000, 15, 5),ndistr=7, newc))
zapsmall(cor(geneZscoresNullGwas))

sum(abs(cor(geneZscoresNullGwas)) >= 0.9)

str(geneZscoresNullGwas)


write.table(geneZscoresNullGwas, file = "C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\DEPICT2\\src\\test\\resources\\geneZscoresNullGwas.txt", quote = F, sep = "\t", col.names = NA)

collapse1 <- matrix(nrow = 4, ncol = 10000)
rownames(collapse1) <- c("1_2_4", "3", "5_8", "6_7")

collapse1[1,] <- apply(geneZscoresNullGwas[c(1,2,4),], 2, mean) 
collapse1[2,] <- geneZscoresNullGwas[3,]
collapse1[3,] <- apply(geneZscoresNullGwas[c(5,8),], 2, mean) 
collapse1[4,] <- apply(geneZscoresNullGwas[c(6,7),], 2, mean) 

write.table(collapse1, file = "C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\DEPICT2\\src\\test\\resources\\mergedMetaGenes1.txt", quote = F, sep = "\t", col.names = NA)


collapse2 <- matrix(nrow = 4, ncol = 10000)
rownames(collapse2) <- c("1_2_4", "3", "5_8", "6_7")

collapse2[1,] <- apply(geneZscoresNullGwas[c(1,2,4),], 2, function(x){sum(x)/sqrt(length(x))}) 
collapse2[2,] <- geneZscoresNullGwas[3,]
collapse2[3,] <- apply(geneZscoresNullGwas[c(5,8),], 2, function(x){sum(x)/sqrt(length(x))})
collapse2[4,] <- apply(geneZscoresNullGwas[c(6,7),], 2, function(x){sum(x)/sqrt(length(x))})

  
write.table(collapse2, file = "C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\DEPICT2\\src\\test\\resources\\mergedMetaGenes2.txt", quote = F, sep = "\t", col.names = NA)
