#Code to generate matrix from: http://r.789695.n4.nabble.com/Re-Generate-a-serie-of-new-vars-that-correlate-with-existingvar-td822187.html


# create the initial x variable
x1 <- rnorm(100, 15, 5)

# x2, x3, and x4 in a matrix, these will be modified to meet the criteria
x234 <- scale(matrix( rnorm(300), ncol=3 ))

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

# create new correlation structure (zeros can be replaced with other rvals)
newc <- matrix( 
  c(1  , 0.4, 0.5, 0.6, 
    0.4, 1  , 0.003  , 0  ,
    0.5, 0.003  , 1  , 0.8  ,
    0.6, 0  , 0.8  , 1  ), ncol=4 )

# check that it is positive definite
eigen(newc)

chol2 <- chol(newc)

randomGwasScores <- t(newx %*% chol2 * sd(x1) + mean(x1))

# verify success
mean(x1)
colMeans(randomGwasScores)

sd(x1)
apply(randomGwasScores, 2, sd)

randomGwasScores[,3] <- randomGwasScores[,3]*-1

zapsmall(cor(randomGwasScores))

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(randomGwasScores, upper.panel = panel.cor)

all.equal(x1, randomGwasScores[,1])

str(randomGwasScores)

rownames(randomGwasScores) <- paste0("Gene", 1:4)
colnames(randomGwasScores) <- paste0("RandomGwas", 1:100)

randomGwasScoresCor <- cor(t(randomGwasScores))

randomGwasScoresCorInverse <- solve(randomGwasScoresCor)

gwasGeneScores <- matrix(c(0.3,0.4,0.4,0.6,0,0.1,0.3,0.1), ncol = 2, nrow = 4)

#gwasGeneScores <- matrix(c(0.3,0.4,0.4,0.6), ncol = 1, nrow = 4)

m2 <- apply(gwasGeneScores, 2, mean)
sd2 <- apply(gwasGeneScores, 2, sd)

gwasGeneScores2 <- gwasGeneScores

for(i in 1:2){
  gwasGeneScores2[,i] <- (gwasGeneScores2[,i] - m2[i]) / sd2[i]
}



pathwayGeneScores <- matrix(c(0.2,0.3,0.4,0.5,1,0.5,0.2,0.3,0.1,0.2,0.3,0.4), ncol = 3, nrow = 4)
#pathwayGeneScores <- matrix(c(0.2,0.3,0.4,0.5), ncol = 1, nrow = 4)


m3 <- apply(pathwayGeneScores, 2, mean)
sd3 <- apply(pathwayGeneScores, 2, sd)

pathwayGeneScores2 <- pathwayGeneScores

for(i in 1:3){
  pathwayGeneScores2[,i] <- (pathwayGeneScores2[,i] - m3[i]) / sd3[i]
}

sd(pathwayGeneScores2[,1])



t(gwasGeneScores2) %*% randomGwasScoresCorInverse

gwas1 <- as.matrix(gwasGeneScores2[,2])

solve(t(gwas1) %*% randomGwasScoresCorInverse %*% gwas1) %*% (t(gwas1) %*% randomGwasScoresCorInverse) %*% as.matrix(pathwayGeneScores2[,1])



gwasGeneScores2[,1] %*% randomGwasScoresCorInverse %*% pathwayGeneScores2[,1]



x <- (t(gwasGeneScores2) %*% randomGwasScoresCorInverse)[1,]
y <- gwasGeneScores2[,1]

sum(x*y)



setwd("C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\DEPICT2\\src\\test\\resources\\")

write.table(randomGwasScoresCorInverse, file = "invCorMatrix.txt", quote = F, sep = "\t", col.names = NA)
write.table(gwasGeneScores2, file = "gwasGeneScores.txt", quote = F, sep = "\t", col.names = NA)
write.table(pathwayGeneScores2, file = "pathwayGeneScores.txt", quote = F, sep = "\t", col.names = NA)

gls(gwas1 ~ pathwayGeneScores2[,1] + pathwayGeneScores2[,2] + pathwayGeneScores2[,3], correlation = randomGwasScoresCor)



fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary, correlation = corAR1(form = ~ 1 | Mare))
fm1
