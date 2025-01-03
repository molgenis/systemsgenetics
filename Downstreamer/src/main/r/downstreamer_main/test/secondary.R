

library(readr)

#change X1 in case of specified header for row names. Use ...1 if first colun has no ID
colTypes <- cols( .default = col_double(),  `...1` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/users/umcg-pdeelen/tp53genotypes_LL.txt.dosages.txt", delim = "\t", quote = "", col_types = colTypes)
dosages <- as.matrix(table_tmp[,-1])
rownames(dosages) <- table_tmp[,1][[1]]
rm(table_tmp)

dosages <- dosages[1000:1999,]

snp1Index <- match("rs17881035", rownames(dosages))
snp2Index <- match("rs8073498", rownames(dosages))
snp3Index <- match("rs9894946", rownames(dosages))



effect1 <- dosages[snp1Index,] * 0.42
effect2 <- dosages[snp2Index,] * 0.2
effect3 <- dosages[snp3Index,] * 0.1

set.seed(42);
noise <- rnorm(ncol(dosages), sd = 0.2)

exp <- effect1 + noise
exp <- effect1 + effect2 + noise
exp <- effect1 + effect2 + effect3 + noise


layout(matrix(1:3, nrow = 1))
plot(dosages[snp1Index,], exp, xlab = "rs17881035", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8, ylab = "Simulated expression")
plot(dosages[snp2Index,], exp, xlab = "rs8073498", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8, ylab = "Simulated expression")
plot(dosages[snp3Index,], exp, xlab = "rs9894946", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8, ylab = "Simulated expression")



cor(exp,dosages[snp1Index,])

summary(lm(exp ~ dosages[snp1Index,]))
cor.test(exp, predictBasedOnSnp)


summary(lm(exp ~ dosages[snp2Index,]))
summary(lm(exp ~ dosages[snp3Index,]))

summary(lm(exp ~ dosages[snp1Index,] + dosages[snp2Index,] + dosages[snp3Index,]))



summary(lm(exp ~ dosages[snp1Index,] + dosages[snp2Index,] + dosages[snp3Index,] + dosages[367,]))


cor.test(exp,dosages[snp1Index,])$p.value
cor.test(exp,dosages[snp2Index,])$p.value
cor.test(exp,dosages[snp3Index,])$p.value

rounds <- 10
pvalues <- matrix(NA, nrow = nrow(dosages), ncol = rounds)
betas <- matrix(NA, nrow = nrow(dosages), ncol = rounds)
roundStats <- matrix(NA, nrow = rounds, ncol = 3)
colnames(roundStats) <- c("variantIndex", "min p-value", "beta")

previousResiduals <- exp
for(r in 1:rounds){
  
  roundMinP = 1
  roundResiduals = NA
  
  for(v in 1:nrow(dosages)){
    
    fit <- lm(previousResiduals ~ dosages[v,])
    snpP <- summary(fit)$coefficients[2,4]
    snpBeta <- summary(fit)$coefficients[2,1]
    pvalues[v,r] <- snpP
    betas[v,r] <- snpBeta
    
    if(snpP < roundMinP){
      roundResiduals <- residuals(fit)
      roundMinP <- snpP
      
      roundStats[r,1] <- v
      roundStats[r,2] <- snpP
      roundStats[r,3] <- snpBeta
      
    }
    
  }
  
  print(paste("Round", r, "min pvalue", roundMinP))
  previousResiduals <- roundResiduals
  
}

round1Res <- data.frame(snp = rownames(dosages), beta = betas[,1], pvalue = pvalues[,1])
str(round1Res)
write.table(round1Res, file = "/home/umcg-pdeelen/simulatedEqtl.txt", sep = "\t", quote = F, row.names = F)


str(exp)
str(as.data.frame(t(dosages[roundStats[1:8,1],])))
summary(lm(exp ~ ., data = as.data.frame(t(dosages[roundStats[,1],]))))
predictBasedOnRounds <- predict(lm(exp ~ ., data = as.data.frame(t(dosages[roundStats[,1],]))))
cor.test(exp, predictBasedOnRounds)


plot(dosages[snp1Index,], previousResiduals)
plot(dosages[snp3Index,], previousResiduals)

cor.test(exp,dosages[snp1Index,])$p.value
cor.test(exp,dosages[snp2Index,])$p.value
cor.test(exp,dosages[snp3Index,])$p.value

roundStats


library("RColorBrewer")
snpColPallete <- brewer.pal(3, name = "Accent")

367 
cor(t(dosages[c(snp1Index, snp2Index, snp3Index, 367),]))

layout(matrix(1:3, nrow = 1))
plot(dosages[snp1Index,], dosages[snp2Index,], xlab = "rs17881035 (effect 1)", ylab = "rs8073498 (effect 2)", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)
plot(dosages[snp1Index,], dosages[367,], xlab = "rs17881035 (effect 1)", ylab = "rs1614984 (detected effect)", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)
plot(dosages[snp2Index,], dosages[367,], xlab = "rs8073498 (effect 2)", ylab = "rs1614984 (detected effect)", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)



snpCol <- rep(adjustcolor("grey", alpha.f = 0.5),nrow(dosages))
snpCol[snp1Index] <- snpColPallete[1]
snpCol[snp2Index] <- snpColPallete[2]
snpCol[snp3Index] <- snpColPallete[3]

snpSize <- rep(1, nrow(dosages))
snpSize[snp1Index] <- 2
snpSize[snp2Index] <- 2
snpSize[snp3Index] <- 2

layout(1)
plot(betas[,1], pch = 16, col=snpCol, cex = snpSize, ylab = "Fitted beta", main = "Round 1", las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(betas[,2], pch = 16, col=snpCol, cex = snpSize, ylab = "Fitted beta", main = "Round 2", las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(betas[,3], pch = 16, col=snpCol, cex = snpSize, ylab = "Fitted beta", main = "Round 3", las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(betas[,4], pch = 16, col=snpCol, cex = snpSize, ylab = "Fitted beta", main = "Round 4", las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(betas[,5], pch = 16, col=snpCol, cex = snpSize, ylab = "Fitted beta", main = "Round 5", las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(betas[,6], pch = 16, col=snpCol, cex = snpSize, ylab = "Fitted beta", main = "Round 5", las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))


plot(-log10(pvalues[,1]), main = "Round 1", ylab = "-log10(p-value)", pch = 16, col=snpCol, cex = snpSize, las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(-log10(pvalues[,2]), main = "Round 2", ylab = "-log10(p-value)", pch = 16, col=snpCol, cex = snpSize, las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(-log10(pvalues[,3]), main = "Round 3", ylab = "-log10(p-value)", pch = 16, col=snpCol, cex = snpSize, las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(-log10(pvalues[,4]), main = "Round 4", ylab = "-log10(p-value)", pch = 16, col=snpCol, cex = snpSize, las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(-log10(pvalues[,5]), main = "Round 5", ylab = "-log10(p-value)", pch = 16, col=snpCol, cex = snpSize, las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(-log10(pvalues[,6]), main = "Round 6", ylab = "-log10(p-value)", pch = 16, col=snpCol, cex = snpSize, las = 2)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))


#dosagesT <- t(dosages)
#https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
dosagesScale = dosages - rowMeans(dosages);
# Standardize each variable
dosagesScale = dosagesScale / sqrt(rowSums(dosagesScale^2));   


dosageSvd <- svd(t(dosagesScale), nu = 500, nv = 500)

x <- eigen(cor(t(dosages)))


str(dosageSvd)
str(x)

str(dosageSvd$u)

str(dosageSvd)
plot(dosageSvd$d)
head(dosageSvd$d^2)

eigenFit <- summary(lm(exp ~ ., data = as.data.frame(dosageSvd$u)))

eigenFit2Coef <- eigenFit$coefficients[-1,]
eigenFit2Coef[eigenFit2Coef[,4]>0.01,1] <- 0
eigenFit2Beta <- eigenFit2Coef[,1]
sum(eigenFit2Beta != 0)


eqtlSignal <-  dosageSvd$v %*% eigenFit2Beta
str(eqtlSignal)
str( diag(eigenFit$coefficients[2:501,1]) )
str( t(dosageSvd$v))

plot(eqtlSignal, pch = 16, col=snpCol, cex = snpSize, ylab = "Reconstructed regulation score", main = "Joint modeling of locus")
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

layout(matrix(1:3, nrow =1))
modeledExpression <-  dosageSvd$u %*% eigenFit2Beta
str(modeledExpression)
plot(exp, modeledExpression, main = "PCA method", xlab = "Simulated expression", ylab = "Predicted expression using PCA", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)
cor.test(exp, modeledExpression)

plot(exp, predictBasedOnSnp, main = "3 known top SNPs", xlab = "Simulated expression", ylab = "Predicted expression using 3 SNPs", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)


cor.test(exp, predictBasedOnRounds)
plot(exp, predictBasedOnRounds, main = "8 regression rounds", xlab = "Predicted expression using ", ylab = "Predicted expression using top indepdent SNPs", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)






str(dosageSvd)

reconstructedSvd <- dosageSvd$u %*% (diag(dosageSvd$d[1:500]) %*% t(dosageSvd$v))

x <- t(diag(1/dosageSvd$d[1:500]) %*% t(dosageSvd$v) %*% dosagesScale)

str(x)
str(dosageSvd$u)

plot(x[200,], dosageSvd$u[200,])

str(t(reconstructedSvd))
str(dosages)

plot(t(reconstructedSvd)[,1], dosages[,1])

plot(t(reconstructedSvd)[,1], dosagesScale[,1])


colTypes <- cols( .default = col_double(),  `...1` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/users/umcg-pdeelen/tp53genotypes_1000g.txt.dosages.txt", delim = "\t", quote = "", col_types = colTypes)
dosages1000g <- as.matrix(table_tmp[,-1])
rownames(dosages1000g) <- table_tmp[,1][[1]]
rm(table_tmp)


rownames(dosages)[!rownames(dosages) %in% rownames(dosages1000g)]
all(rownames(dosages) %in% rownames(dosages1000g))
dosages1000g <- dosages1000g[rownames(dosages),]

all(rownames(dosages1000g) == rownames(dosages))
 

#https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
dosages1000gScale = dosages1000g - rowMeans(dosages1000g);
# Standardize each variable
sdScale <-  sqrt(rowSums(dosages1000gScale^2))
dosages1000gScale = dosages1000gScale / sdScale;   


dosage1000gSvd <- svd(t(dosages1000gScale), nu = 50, nv = 50)

str(dosage1000gSvd)

x <- t(dosages1000gScale) %*% dosage1000gSvd$v %*% diag(1/dosage1000gSvd$d[1:50])



plot(x[1,], dosage1000gSvd$u[1,])



str(dosage1000gSvd$v)
plot(log(dosage1000gSvd))

dosagesScale2 = (dosages - rowMeans(dosages1000g)) / sdScale

llMapped <- t(dosagesScale2) %*% dosage1000gSvd$v %*% diag(1/dosage1000gSvd$d[1:50])

plot(llMapped[,1], dosageSvd$u[,1])

cor.test(llMapped[,1], dosageSvd$u[,1])
cor.test(llMapped[,1], exp)
cor.test(dosageSvd$u[,1], exp)

a <- cor( dosage1000gSvd$v ,  dosageSvd$v )
library(pheatmap)
pheatmap(a[1:50, 1:50], scale = "none", cluster_cols = F, cluster_rows = F)

eigenFitB <- summary(lm(exp ~ ., data = as.data.frame(llMapped)))

eigenFitB2Coef <- eigenFitB$coefficients[-1,]
eigenFitB2Coef[eigenFitB2Coef[,4]>0.01,1] <- 0
eigenFitB2Beta <- eigenFitB2Coef[,1]
sum(eigenFitB2Beta != 0)

layout(1)
eqtlSignal <-  llMapped %*% eigenFitB2Beta

plot(eqtlSignal, pch = 16, col=snpCol, cex = snpSize, ylab = "Reconstructed regulation score", main = "Joint modeling of locus")
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

layout(matrix(1:3, nrow =1))

modeledExpression <- llMapped %*% eigenFitB2Beta
plot(exp, modeledExpression, main = "PCA method ext ref", xlab = "Simulated expression", ylab = "Predicted expression using PCA of 1000G", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)
cor.test(exp, modeledExpression)


modeledExpression <-  dosageSvd$u %*% eigenFit2Beta
plot(exp, modeledExpression, main = "PCA method", xlab = "Simulated expression", ylab = "Predicted expression using PCA", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)
cor.test(exp, modeledExpression)

plot(exp, predictBasedOnSnp, main = "3 known top SNPs", xlab = "Simulated expression", ylab = "Predicted expression using 3 SNPs", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)


cor.test(exp, predictBasedOnRounds)





