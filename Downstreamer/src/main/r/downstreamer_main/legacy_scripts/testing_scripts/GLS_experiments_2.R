# See http://www.biostat.jhsph.edu/~iruczins/teaching/jf/ch5.pdf

data(longley)
g <- lm(Employed~Population,data=longley)
summary(g, cor=T)

# GLS
y     <- as.matrix(longley$Employed)
x     <- as.matrix(model.matrix(g)[,2])

# Construct sigma matrix
Sigma <- diag(16)
Sigma <- 0.31041 ^ abs(row(Sigma)-col(Sigma))

rownames(Sigma) <- 1:16
rownames(x) <- 1:16
rownames(y) <- 1:16

colnames(Sigma) <- 1:16
colnames(x) <- c("intercept", "population")
colnames(y) <- c("y")

write.table(Sigma, file="~/Desktop/sigma_test.tsv", col.names=T, row.names=T, quote=F, sep="\t")
write.table(x, file="~/Desktop/x_test.tsv", col.names=T, row.names=T, quote=F, sep="\t")
write.table(y, file="~/Desktop/y_test.tsv", col.names=T, row.names=T, quote=F, sep="\t")

Sigi <- solve(Sigma)
xtxi <- solve(t(x) %*% Sigi %*% x) #B1
beta <- xtxi %*% t(x) %*% Sigi %*% y


xtx    <- t(x) %*% Sigi %*% x 
beta2  <- (t(x) %*% Sigi %*% y) / xtx


# Calculate SE
res <- y - x %*% beta
sig <- sqrt(sum(res^2) / g$df)
se  <- sqrt(diag(xtxi)) * sig

se2 <- sig / sqrt(xtx)

# Calculate p
tstats <- abs(beta / se2)
2 * pt(tstats[1], df=16, lower=F)


library(nlme)
g <- gls(Employed ~ GNP + Population, correlation=corAR1(form=~Year),data=longley)
summary(g)
