
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

finalx <- newx %*% chol2 * sd(x1) + mean(x1)

# verify success
mean(x1)
colMeans(finalx)

sd(x1)
apply(finalx, 2, sd)

finalx[,3] <- finalx[,3]*-1

zapsmall(cor(finalx))

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

pairs(finalx, upper.panel = panel.cor)

all.equal(x1, finalx[,1])


write.table(finalx, file = "C:\\Users\\patri\\Documents\\GitHub\\systemsgenetics\\genetica-libraries\\src\\test\\resources\\testMatrix.txt", quote = F, sep = "\t", col.names = NA)
