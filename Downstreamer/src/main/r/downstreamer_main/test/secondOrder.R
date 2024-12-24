

colTypes <- cols( .default = col_double(),  `...1` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/users/umcg-pdeelen/nod2genotypes_1000g.txt.dosages.txt", delim = "\t", quote = "", col_types = colTypes)
dosages1000g <- as.matrix(table_tmp[,-1])
rownames(dosages1000g) <- table_tmp[,1][[1]]
rm(table_tmp)


#change X1 in case of specified header for row names. Use ...1 if first colun has no ID
colTypes <- cols( .default = col_double(),  `...1` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/users/umcg-pdeelen/nod2genotypes_LL.txt.dosages.txt", delim = "\t", quote = "", col_types = colTypes)
dosagesLL <- as.matrix(table_tmp[,-1])
rownames(dosagesLL) <- table_tmp[,1][[1]]
rm(table_tmp)

dosagesLL <- dosagesLL[rownames(dosagesLL) %in% rownames(dosages1000g),]
colTypes <- cols( .default = col_double(),  `...1` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/users/umcg-pdeelen/ikzf1genotypes_LLarray.txt.dosages.txt", delim = "\t", quote = "", col_types = colTypes)
dosagesLLa <- as.matrix(table_tmp[,-1])
rownames(dosagesLLa) <- table_tmp[,1][[1]]
rm(table_tmp)


dosagesLLa <- dosagesLLa[rownames(dosagesLLa) %in% rownames(dosages1000g),]

str(dosagesLL)

#change X1 in case of specified header for row names. Use ...1 if first colun has no ID
colTypes <- cols( .default = col_double(),  `...1` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/users/umcg-pdeelen/nod2genotypes_RS.txt.dosages.txt", delim = "\t", quote = "", col_types = colTypes)
dosagesRS <- as.matrix(table_tmp[,-1])
rownames(dosagesRS) <- table_tmp[,1][[1]]
rm(table_tmp)


colTypes <- cols( .default = col_double(),  `ID` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426/1_DataQC/out/LL/outputfolder_exp/exp_data_QCd/exp_data_preprocessed.txt", delim = "\t", quote = "", col_types = colTypes)
expLL <- as.matrix(table_tmp[,-1])
rownames(expLL) <- table_tmp[,1][[1]]
rm(table_tmp)



colTypes <- cols( .default = col_double(),  `ID` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426/1_DataQC/out/RS/outputfolder_exp/exp_data_QCd/exp_data_preprocessed.txt", delim = "\t", quote = "", col_types = colTypes)
expRS <- as.matrix(table_tmp[,-1])
rownames(expRS) <- table_tmp[,1][[1]]
rm(table_tmp)


all(rownames(expLL) == colnames(dosagesLL))
all(rownames(expRS) == colnames(dosagesRS))

dosagesLLa <- dosagesLLa[,match(rownames(expLL), colnames(dosagesLLa))]
all(rownames(expLL) == colnames(dosagesLLa))


str(dosagesLL)
str(dosages1000g)

#https://stackoverflow.com/questions/18964837/fast-correlation-in-r-using-c-and-parallelization/18965892#18965892
center <- rowMeans(dosages1000g);
dosages1000gScale = dosages1000g - center
scale <-  sqrt(rowSums(dosages1000gScale^2))
dosages1000gScale = dosages1000gScale / scale;   


dosage1000gSvd <- svd(t(dosages1000gScale), nu = 100, nv = 100)
rownames(dosage1000gSvd$v) <- rownames(dosages1000gScale)
compsToUse <- 50


plot(log(dosage1000gSvd$d[1:compsToUse]))
plot(dosage1000gSvd$d[1:compsToUse] * 100/ sum(dosage1000gSvd$d))
plot(cumsum(dosage1000gSvd$d[1:compsToUse] * 100/ sum(dosage1000gSvd$d)))


dosagesLLScale = (dosagesLL - center[rownames(dosagesLL)]) / scale[rownames(dosagesLL)]
dosagesLLaScale = (dosagesLLa - center[rownames(dosagesLLa)]) / scale[rownames(dosagesLLa)]

llMapped <- t(dosagesLLScale) %*% dosage1000gSvd$v[rownames(dosagesLL),1:compsToUse] %*% diag(1/dosage1000gSvd$d[1:compsToUse])
llaMapped <- t(dosagesLLaScale) %*% dosage1000gSvd$v[rownames(dosagesLLa),1:compsToUse] %*% diag(1/dosage1000gSvd$d[1:compsToUse])

plot(llaMapped[,70], llMapped[,70])

plot(diag(cor(llaMapped, llMapped)))
#tp53 	ENSG00000141510
#IKZF1  ENSG00000185811
#NOD2   ENSG00000167207
#STX3   ENSG00000166900

gene <- "ENSG00000167207"

eigenFitI <- summary(lm(expLL[,gene] ~ ., data = as.data.frame(llMapped[,1:compsToUse])))
eigenFitA <- summary(lm(expLL[,gene] ~ ., data = as.data.frame(llaMapped[,1:compsToUse])))
#eigenFit <- summary(lm(expLL[,gene] ~ ., data = as.data.frame(llaMapped[,17])))

eigenFitABetaCoef <- eigenFitA$coefficients[-1,]
eigenFitABetaCoef[eigenFitABetaCoef[,4]>0.01,1] <- 0
eigenFitABeta <- eigenFitABetaCoef[,1]
sum(eigenFitABeta != 0)


eigenFitIBetaCoef <- eigenFitI$coefficients[-1,]
#eigenFitIBetaCoef[eigenFitIBetaCoef[,4]>0.01,1] <- 0
eigenFitIBeta <- eigenFitIBetaCoef[,1]
sum(eigenFitIBeta != 0)

layout(1)
eqtlSignal <-  llMapped[,1:compsToUse] %*% eigenFitIBeta

plot(eqtlSignal, pch = 16, ylab = "Reconstructed regulation score", main = "Joint modeling of locus")

layout(matrix(1:3, nrow =1))

modeledExpression <-  llMapped[,1:compsToUse]  %*% eigenFitIBeta
plot(expLL[,gene], modeledExpression, main = "PCA method ext ref", xlab = "Simulated expression", ylab = "Predicted expression using PCA of 1000G", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)
cor.test(expLL[,gene], modeledExpression)

modeledExpression <-  llMapped[,1:50]  %*% eigenFitIBeta
plot(expLL[,gene], modeledExpression, main = "PCA method ext ref", xlab = "Simulated expression", ylab = "Predicted expression using PCA of 1000G", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)
cor.test(expLL[,gene], modeledExpression)


modeledExpression <-  dosageSvd$u %*% eigenFit2Beta
plot(exp, modeledExpression, main = "PCA method", xlab = "Simulated expression", ylab = "Predicted expression using PCA", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)
cor.test(exp, modeledExpression)

plot(exp, predictBasedOnSnp, main = "3 known top SNPs", xlab = "Simulated expression", ylab = "Predicted expression using 3 SNPs", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)


cor.test(exp, predictBasedOnRounds)






rounds <- 10
pvalues <- matrix(NA, nrow = nrow(dosagesLL), ncol = rounds)
betas <- matrix(NA, nrow = nrow(dosagesLL), ncol = rounds)
roundStats <- matrix(NA, nrow = rounds, ncol = 3)
colnames(roundStats) <- c("variantIndex", "min p-value", "beta")

previousResiduals <- expLL[,gene]
for(r in 1:rounds){
  
  roundMinP = 1
  roundResiduals = NA
  
  for(v in 1:nrow(dosagesLL)){
    
    fit <- lm(previousResiduals ~ dosagesLL[v,])
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
roundStats2 <- as.data.frame(roundStats)
roundStats2$snp <- rownames(dosagesLL)[roundStats2$variantIndex]
View(roundStats2)

summary(lm( expLL[,gene] ~ ., data = as.data.frame(t(dosagesLL[unique(roundStats[,1]),]))))
predictBasedOnRounds <- predict(lm( expLL[,gene] ~ ., data = as.data.frame(t(unique(dosagesLL[roundStats[,1],])))))
cor.test( expLL[,gene], predictBasedOnRounds)

str(as.data.frame(t(dosagesLL[roundStats[,1],])))

x <- as.data.frame(t(unique(dosagesLL[roundStats[,1],])))

colnames(x) <- make.names(colnames(x))
summary(lm(as.formula(
  paste('expLL[,gene] ~',paste('poly(',colnames(x),',2)',collapse = ' + '))
) , data = x))
y <- predict(lm(as.formula(
  paste('expLL[,gene] ~',paste('poly(',colnames(x),',2)',collapse = ' + '))
) , data = x))
cor.test( expLL[,gene], y)





center <- rowMeans(dosagesLL);
dosagesLLScale = dosagesLL - center
scale <-  sqrt(rowSums(dosagesLLScale^2))
dosagesLLScale = dosagesLLScale / scale;   


dosageLLSvd <- svd(t(dosagesLLScale), nu = 100, nv = 100)


eigenFit <- summary(lm(expLL[,gene] ~ ., data = as.data.frame(dosageLLSvd$u)))
eigenFitBCoef <- eigenFit$coefficients[-1,]
eigenFitBCoef[eigenFitBCoef[,4]>0.01,1] <- 0
eigenFitBBeta <- eigenFitBCoef[,1]
sum(eigenFitBBeta != 0)

modeledExpression <- dosageLLSvd$u %*% eigenFitBBeta
plot(expLL[,gene], modeledExpression, main = "PCA method ext ref", xlab = "Real expression", ylab = "Predicted expression using PCA of 1000G", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)
cor.test(expLL[,gene], modeledExpression)


cor.test(expLL[,gene], dosagesLLScale[590,])

str(dosageLLSvd$u)
tail(sort(abs(cor(dosagesLLScale[590,], dosageLLSvd$u))))
which.max(abs(cor(dosagesLLScale[590,], dosageLLSvd$u)))

cor.test(expLL[,gene], dosageLLSvd$u[,11])
cor.test(dosageLLSvd$u[,11], dosagesLLScale[590,])

tail(sort(abs(cor(dosageLLSvd$v[,11], dosage1000gSvd$v[rownames(dosagesLL),]))))
which.max(cor(dosageLLSvd$v[,11], dosage1000gSvd$v[rownames(dosagesLL),]))


sort(abs(cor(dosageLLSvd$u[,11], llMapped)))
which.max(cor(dosageLLSvd$u[,11], llMapped))

plot(dosageLLSvd$u[,2], llMapped[,1])
cor.test(dosageLLSvd$u[,2], llMapped[,1])

tail(sort(abs(cor(expLL[,gene], llMapped))))

tail(sort(abs(cor(dosagesLLScale[590,], llMapped))))

which.max(abs(cor(dosagesLLScale[590,], llMapped)))

library(heatmap3)

heatmap3(cor(cbind(expLL[,gene], dosagesLLScale[590,], dosageLLSvd$u[,11], llMapped[,17])), scale = "none", Rowv = NA, Colv = NA, balanceColor = T)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(cbind(expLL[,gene], dosagesLLScale[590,], dosageLLSvd$u[,11], llaMapped[,17], llMapped[,17]),upper.panel = panel.cor)


cor.test(expLL[,gene], llaMapped[,17])
cor.test(expLL[,gene], llMapped[,17])
cor.test(expLL[,gene], dosagesLLScale[590,])







comps <- 50

library(glmnet)
cfit <- cv.glmnet(x = llMapped[,1:comps], y = expLL[,gene])
cfit

plot(cfit) 

str(cfit)

assess.glmnet(cfit, newx = llMapped[,1:comps], newy = expLL[,gene], keep = TRUE, alpha=1, lambda = "1se")

coef(cfit, s = "lambda.min")

predictionsTest <- predict(cfit, s = "lambda.min", newx = llMapped[,1:comps])

plot(expLL[,gene], predictionsTest)
cor.test(expLL[,gene], predictionsTest)







gene2 <- "ENSG00000166900"


plot(expLL[,gene], expLL[,gene2])
cor.test(expLL[,gene], expLL[,gene2])

eigenFitI2 <- summary(lm(expLL[,gene] ~ . * expLL[,gene2], data = as.data.frame(llMapped[,1:compsToUse])))




x <- as.data.frame(t(dosagesLL))
x$gene2 <- expLL[,gene2]
colnames(x) <- make.names(colnames(x))


v <- as.data.frame(t(dosagesRS))
v$gene2 <- expRS[,gene2]
colnames(v) <- make.names(colnames(v))



snps <- make.names(rownames(dosagesLL)[roundStats[1:3,1]])

summary(lm(as.formula(
  paste('expLL[,gene] ~',paste('poly(',snps,',2)',collapse = ' + '))
) , data = x))
y <- predict(lm(as.formula(
  paste('expLL[,gene] ~',paste('poly(',snps,',2)',collapse = ' + '))
) , data = x))
cor.test( expLL[,gene], y)

fit <- lm(as.formula(
  paste('expLL[,gene] ~', paste(snps, "gene2", sep = " * ", collapse = " + "))
) , data = x)
summary(fit)
y <- predict(fit)
cor.test( expLL[,gene], y)

paste('expLL[,gene] ~', paste(snps, "gene2", sep = " * ", collapse = " + "))
paste('expLL[,gene] ~', paste(snps, sep = " + ", collapse = " + "))
paste('expLL[,gene] ~ gene2 + rs2111235 + X16.50714603.C.CTGTG + rs1861757')
paste('expLL[,gene] ~ rs2111235 + X16.50714603.C.CTGTG + rs1861757')
paste('gene2 ~', paste(snps, sep = " * ", collapse = " + "))
paste('expLL[,gene] ~ gene2 * rs2111235 + gene2 * X16.50714603.C.CTGTG + gene2 * rs1861757')

paste(' expLL[,gene] ~ ',paste(paste('poly(',snps,',2)'), "gene2", sep = "*", collapse = ' + '))
paste(' expLL[,gene] ~ ',paste(paste('poly(',snps,',2)'), sep = "*", collapse = ' + '))

fit <- lm(as.formula(
  paste('expLL[,gene] ~', paste(snps, "gene2", sep = " * ", collapse = " + "))
) , data = x)
summary(fit)
y <- predict(fit)
cor.test( expLL[,gene], y)
plot( expLL[,gene], y)


fitRS <- lm(as.formula(
  paste('expRS[,gene] ~', paste(snps, "gene2", sep = " * ", collapse = " + "))
) , data = v)
summary(fitRS)
yRS <- predict(fitRS)
cor.test( expRS[,gene], yRS)
plot( expRS[,gene], yRS)


yv <- predict(fit, newdata = v)
cor.test( expRS[,gene], yv)
plot( expRS[,gene], yv)

fit2 <- lm(as.formula(
  paste(' expLL[,gene] ~ ',paste(paste('poly(',snps,',2)'), sep = "*", collapse = ' + '))
) , data = x)
summary(fit2)
y2 <- predict(fit2)
cor.test( expLL[,gene], y2)
plot( expLL[,gene], y2)


fit3 <- lm(as.formula(
  paste(' expLL[,gene] ~ ',paste(paste('poly(',snps,',2)'), "gene2", sep = "*", collapse = ' + '))
) , data = x)
summary(fit3)
y3 <- predict(fit3)
cor.test( expLL[,gene], y3)
plot( expLL[,gene], y3)

fit3RS <- lm(as.formula(
  paste(' expRS[,gene] ~ ',paste(paste('poly(',snps,',2)'), "gene2", sep = "*", collapse = ' + '))
) , data = v)
summary(fit3RS)
y3RS <- predict(fit3RS)
cor.test( expRS[,gene], y3RS)
plot( expRS[,gene], y3RS)

yv3 <- predict(fit3, newdata = v)
cor.test( expRS[,gene], yv3)
plot( expRS[,gene], yv3)


fit4 <- lm(as.formula(
  paste('expLL[,gene] ~', paste(snps, sep = " + ", collapse = " + "))
) , data = x)
summary(fit4)
y4 <- predict(fit4)
cor.test( expLL[,gene], y4)plot( expLL[,gene], y4)




fit5 <- lm(as.formula(
  paste(' expLL[,gene] ~ ',paste(paste('poly(',snps,',5)'), paste('poly(gene2,3)'), sep = "*", collapse = ' + '))
) , data = x)
summary(fit5)
y5 <- predict(fit5)
cor.test( expLL[,gene], y5)
plot( expLL[,gene], y5)


fit6 <- lm(as.formula(
  paste(' expLL[,gene] ~ poly(rs2111235,2) + poly(X16.50714603.C.CTGTG,2) + poly(rs1861757,2) + gene2 + gene2 * rs2111235  + gene2 * X16.50714603.C.CTGTG + gene2 * rs1861757')
) , data = x)
summary(fit6)
y6 <- predict(fit6)
cor.test( expLL[,gene], y6)
plot( expLL[,gene], y6)



plot(expLL[,gene], y, main = "Model C", xlab = "NOD2", ylab = "NOD2 described by model", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)

plot(expLL[,gene], y3, main = "Model D", xlab = "NOD2", ylab = "NOD2 described by model", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)

layout(matrix(1:2, nrow = 1))
plot(expRS[,gene], yv, main = "Model C validated on RS", xlab = "NOD2", ylab = "NOD2 described by model", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)

plot(expRS[,gene], yv3, main = "Model D validated on RS", xlab = "NOD2", ylab = "NOD2 described by model", pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.8)



anova(fit, fit3)
anova(fit, fit6)
anova(fitRS, fit3RS)

anova(fit4, fit2)


anova(fit3, fit5)

cor.test( y, expLL[,gene2])

z <- residuals(fit)
cor.test( z, expLL[,gene2])
plot( z, expLL[,gene2])

cor.test( expLL[,gene], expLL[,gene2])



