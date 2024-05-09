

library(readr)

#change X1 in case of specified header for row names. Use ...1 if first colun has no ID
colTypes <- cols( .default = col_double(),  `...1` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/users/umcg-pdeelen/tp53genotypes_LL.txt.dosages.txt", delim = "\t", quote = "", col_types = colTypes)
dosages <- as.matrix(table_tmp[,-1])
rownames(dosages) <- table_tmp[,1][[1]]
rm(table_tmp)



str(dosages)
"rs17881035" %in% rownames(dosages)


exp <- dosages["rs17881035",] * 0.42 + rnorm(ncol(dosages))

plot(dosages["rs17881035",], exp)

cor(exp,dosages["rs17881035",])

summary(lm(exp ~ dosages["rs17881035",]))$coefficients[2,4]

cor.test(exp,dosages["rs17881035",])

rounds <- 5
pvalues <- matrix(NA, nrow = nrow(dosages), ncol = rounds)
roundStats <- matrix(NA, nrow = 5, ncol = 3)
colnames(roundStats) <- c("variantIndex", "min p-value", "beta")

previousResiduals <- exp
for(r in 1:rounds){
  
  roundMinP = 1
  roundResiduals = NA
  
  for(v in 1:nrow(dosages)){
    
    fit <- lm(previousResiduals ~ dosages[v,])
    snpP <- summary(fit)$coefficients[2,4]
    pvalues[v,r] <- snpP
    
    if(snpP < roundMinP){
      roundResiduals <- residuals(fit)
      roundMinP <- snpP
      
      roundStats[r,1] <- v
      roundStats[r,2] <- snpP
      roundStats[r,3] <- summary(fit)$coefficients[2,1]
      
    }
    
  }
  
  print(paste("Round", r, "min pvalue", roundMinP))
  previousResiduals <- roundResiduals
  
}
roundStats

plot(-log10(pvalues[,2]), pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5))


