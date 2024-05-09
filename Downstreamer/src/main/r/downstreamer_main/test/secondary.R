

library(readr)

#change X1 in case of specified header for row names. Use ...1 if first colun has no ID
colTypes <- cols( .default = col_double(),  `...1` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/users/umcg-pdeelen/tp53genotypes_LL.txt.dosages.txt", delim = "\t", quote = "", col_types = colTypes)
dosages <- as.matrix(table_tmp[,-1])
rownames(dosages) <- table_tmp[,1][[1]]
rm(table_tmp)

dosages <- dosages[1000:1999,]


str(dosages)


effect1 <- dosages["rs17881035",] * 0.42
effect2 <- dosages["rs8073498",] * 0.2
effect3 <- 0 #dosages["rs1614984",] * 0.1

set.seed(42);
noise <- rnorm(ncol(dosages), sd = 0.2)
exp <- effect1 + noise


exp <- effect1 + effect2 + effect3 + noise

plot(dosages["rs17881035",], exp)
plot(dosages["rs1614984",], exp)


cor(exp,dosages["rs17881035",])

summary(lm(exp ~ dosages["rs17881035",]))
summary(lm(exp ~ dosages["rs8073498",]))
summary(lm(exp ~ dosages["rs1614984",]))

summary(lm(exp ~ dosages["rs17881035",] + dosages["rs8073498",] + dosages["rs1614984",]))

cor.test(exp,dosages["rs17881035",])$p.value
cor.test(exp,dosages["rs8073498",])$p.value
cor.test(exp,dosages["rs1614984",])$p.value

rounds <- 5
pvalues <- matrix(NA, nrow = nrow(dosages), ncol = rounds)
betas <- matrix(NA, nrow = nrow(dosages), ncol = rounds)
roundStats <- matrix(NA, nrow = 5, ncol = 3)
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
roundStats


library("RColorBrewer")
snpColPallete <- brewer.pal(3, name = "Accent")

snp1Index <- match("rs17881035", rownames(dosages))
snp2Index <- match("rs8073498", rownames(dosages))
snp3Index <- match("rs1614984", rownames(dosages))

snpCol <- rep(adjustcolor("grey", alpha.f = 0.5),nrow(dosages))
snpCol[snp1Index] <- snpColPallete[1]
snpCol[snp2Index] <- snpColPallete[2]
snpCol[snp3Index] <- snpColPallete[3]

snpSize <- rep(1, nrow(dosages))
snpSize[snp1Index] <- 2
snpSize[snp2Index] <- 2
snpSize[snp3Index] <- 2

plot(betas[,1], pch = 16, col=snpCol, cex = snpSize)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(betas[,2], pch = 16, col=snpCol, cex = snpSize)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(betas[,3], pch = 16, col=snpCol, cex = snpSize)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))


plot(-log10(pvalues[,1]), pch = 16, col=snpCol, cex = snpSize)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(-log10(pvalues[,2]), pch = 16, col=snpCol, cex = snpSize)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))

plot(-log10(pvalues[,3]), pch = 16, col=snpCol, cex = snpSize)
legend("topright", fill = snpColPallete, legend = c("Primairy (b=0.42)","Secondary (b=0.2)","Tertiary (b=0.1)"))
