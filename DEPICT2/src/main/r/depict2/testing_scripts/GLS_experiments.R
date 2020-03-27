# See http://www.biostat.jhsph.edu/~iruczins/teaching/jf/ch5.pdf

data(longley)
g <- lm(Employed~GNP+Population,data=longley)
summary(g, cor=T)


# GLS
x     <- model.matrix(g)

# Construct sigma matrix
Sigma <- diag(16)
#Sigma <- 0.31041^abs(row(Sigma)-col(Sigma))
Sigma <- 0.3662414^abs(row(Sigma)-col(Sigma))


Sigi <- solve(Sigma)
xtxi <- solve(t(x) %*% Sigi %*% x )
beta <- xtxi %*% t(x) %*% Sigi %*% longley$Empl

# Calculate SE
res <- longley$Empl - x %*% beta
sig <- sqrt(sum(res^2) / g$df)
se  <- sqrt(diag(xtxi))*sig

# Calculate p
tstats <- abs(beta[,1] / se)
2 * pt(tstats[3], df=16, lower=F)


library(nlme)
g <- gls(Employed ~ GNP + Population, correlation=corAR1(form=~Year),data=longley)
summary(g)
