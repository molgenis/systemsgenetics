#!/usr/bin/env Rscript
library(data.table)
library(simGWAS)
library(mvtnorm)
library(corpcor)

# Function definitions
#---------------------------------------------------------------
# Function to read haplotype block
gethap <- function(i) {
  y=fread(ldd$comm[i], skip="#CHROM")
  ha <- simGWAS:::vcf2haps(as.matrix(y[,-c(1:9)]))
  rownames(ha) <- paste0("pos",y$POS)
  t(ha)
}

# Funtion to calculate LD
cor2 <- function (x) {
  1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}

# Load arguments
#---------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
CHR <- args[1] # chromosome
RUN <- args[2]

# Hardcoded for now
NCSE <- NCTL <- N <- 300000 #500000 #10000 #number of cases and controls
NSIM <- 1 #number of simulations
NCV  <- 5  #number of causal variants
DIR  <- "/data/umcg-obakker/projects/simulating_gwas/"
MAF  <- 0.05
NLD  <- 0.5 # Number of LD blocks to simulate a causal effect for. If < 0 use all if > 0 < 1 takes a percentage of LD blocks

print(paste0('Number of cases and controls: ', N))
print(paste0('Number of simulations: ', NSIM))
print(paste0('Number of causal variants: ', NCV))
#print(paste0('Number of LD blocks to simulate: ', NLD))
print(paste0('Chromosome: ', CHR))
print(paste0('Working directory: ', DIR))
print(paste0('Considering only SNPs at MAF > ', MAF))
#---------------------------------------------------------------

# Identify LD and VCF and legend files (European)
# CURRENTLY HARDOCODED TO 1KG EUR AND PREDIFNED LD BLOCKS
file.vcf <- paste0(DIR, 'reference_1000G/EUR.chr', CHR, '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.recode.vcf.gz')
file.ld  <- paste0(DIR, 'reference_1000G/LD_blocks/EUR/fourier_ls-chr',CHR,'.bed')
file.snp <- paste0('zcat ', DIR, 'reference_1000G/1000GP_Phase3_chr',CHR,'.legend.gz')

# Read in SNP and LD blocks data
snp_info <- fread(file.snp)
snp_info <- snp_info[, c('id', 'position')]
ldd      <- fread(file.ld)


if (NLD > nrow(ldd)){
  NLD <- nrow(ldd)
}

if (NLD > 0 & NLD < 1) {
  NLD <- round(NLD * nrow(ldd))
}

print(paste0('Simulating causal variants in ', NLD, ' ld blocks'))
ldblocks.to.use <- sample(1:nrow(ldd), NLD)

#if (NLD > 0) {
#    # Randomly sample X LD blocks for this chromosome
#   ldd <- ldd[sample(1:nrow(ldd), NLD)]
#}

# Split VCF by LD blocks, define command to retrieve haplotyes from VCF. bi-allelic only MAF < x
ldd[,blocknum:=1:.N]
ldd[,comm:=paste0("bcftools view ",file.vcf,
                  " --min-af ", MAF, ":minor --max-alleles 2 --min-alleles 2 -r ",CHR,":",start,"-",stop," -Ov")]

# Make some empty vectors to store the output
simv     <- simz <- simbeta <- simpval <- allres <- vars <- vector("list", nrow(ldd))
measures <- c('expected_z', paste0('simulated_z_',1:NSIM), paste0('simulated_beta_',1:NSIM), paste0('simulated_pval_',1:NSIM))

for(i in 1:nrow(ldd)) {

  # Choose odds ratios simulate either a strong effect, or a very, very weak one
  #if (i %in% ldblocks.to.use) {
  #   #g1 <- sample(seq(from = 1.01, to = 1.4, by = 0.01), NCV, replace = T)
  #   g1 <- sample(rnorm(mean = 1, sd=0.05, n=10000), NCV, replace=T)
  #} else {
  #   g1 <- sample(rnorm(mean = 1, sd=0.005, n=10000), NCV, replace=T)
  #}

  # New strategy, have one combined distribution to sample from
  g1 <- sample(c(rnorm(mean = 1, sd=0.05, n=100000), rnorm(mean = 1, sd=0.012, n=100000)), NCV, replace=T)

  # Load the haplotypes for this LD block
  h <- gethap(i)
  print(paste0('Haplotypes for block ', i,' read in correctly'))

  # Select only variants that have some variance, probably redundant?
  use              <- apply(h, 2, var) > 0
  h                <- h[, use, drop = FALSE]

  # Calculte genotype frequencies?
  freq             <- as.data.frame(h+1)
  freq$Probability <- 1 / nrow(freq)

  # Calculate LD matrix for the locus
  LD               <- cor2(h)
  print('LD computed')

  # Select NCV random causal variants
  CV <- sample(colnames(h), NCV)
  FP <- make_GenoProbList(snps=colnames(h), W=CV, freq=freq)
  print('Genotype probabilities calculated')

  # Find expected Z-scores
  exp_z <- simGWAS:::est_statistic(N0=NCTL, # number of controls
                                   N1=NCSE, # number of cases
                                   snps=colnames(h), # column names in freq of SNPs for which Z scores should be generated
                                   W=CV, # causal variants, subset of snps
                                   gamma.W=log(g1), # odds ratios
                                   freq=freq, # reference haplotypes
                                   GenoProbList=FP) # FP above
  print('Z-scores calculated')

  # Simulate Z-scores
  simz[[i]] <- rmvnorm(n = NSIM, mean = exp_z, sigma = LD)
  print('Z-score simulations calculated')

  # Simulate accompanying beta values
  simv[[i]] <- simulated_vbeta(N0=NCTL, # number of controls
                               N1=NCSE, # number of cases
                               snps=colnames(h), # column names in freq of SNPs for which Z scores should be generated
                               W=CV, # causal variants, subset of snps
                               gamma.W=log(g1), # odds ratios
                               freq=freq, # reference haplotypes
                               GenoProbList=FP,
                               nrep=NSIM)

  # Calculate the beta's
  simbeta[[i]] <- simz[[i]] * sqrt(simv[[i]])

  # Replaced below with above as this should be the same.
  #simv[[i]]    <- 1 / tmpsimv
  #simbeta[[i]] <- simz[[i]] / sqrt(simv[[i]])

  # Calculate the pvalues
  simpval[[i]] <- pnorm(abs(simz[[i]]), lower.tail=FALSE)*2

  print('Beta simulations calculated')

  # Save results
  allres[[i]]           <- rbind(exp_z, simz[[i]], simbeta[[i]], simpval[[i]])
  colnames(allres[[i]]) <- colnames(h)
  allres[[i]]           <- as.data.frame(t(allres[[i]]))
  colnames(allres[[i]]) <- measures

  # Save causal variant positions
  vars[[i]]             <- as.data.frame(matrix(c(rep(i, length(CV)), gsub('pos', paste0(CHR,':'), CV), g1), ncol=3, nrow=length(CV), byrow=F))
  colnames(vars[[i]])   <- c("LD_block", "Variant", "OR")

}

# Combine results
dm <- do.call("rbind", allres)
dm$position <- gsub('pos', '',rownames(dm))

# Merge with SNP info to get rs IDs and allele info, reorder
dm           <- merge(dm, snp_info, by = 'position')
dm$chr       <- rep(CHR, nrow(dm))
dm$rs        <- gsub("(.*)\\:.*\\:(.*)\\:(.*)", "\\1", dm$id)
#dm$reference <- gsub("(.*)\\:.*\\:(.*)\\:(.*)", "\\2", dm$id)
#dm$effect    <- gsub("(.*)\\:.*\\:(.*)\\:(.*)", "\\3", dm$id)
#dm           <- dm[,c(ncol(dm), ncol(dm)-1, 2:(ncol(dm)-2))]

# Make one file with the causal variants per block
causal_variants <- do.call("rbind", vars)

dir.create(paste0(DIR, 'output/', RUN))

write.table(dm[,c("rs","simulated_pval_1")], paste0(DIR, 'output/', RUN ,'/simgwas_', NCV, '-cv_', NLD ,'-LD-blocks_chr', CHR, '.snppval.txt'),
            sep = '\t', quote = F, col.names = T, row.names = F)

write.table(dm, paste0(DIR, 'output/', RUN ,'/simgwas_', NCV, '-cv_', NLD ,'-LD-blocks_chr', CHR, '.assoc.txt'),
            sep = '\t', quote = F, col.names = T, row.names = F)

write.table(causal_variants, paste0(DIR, 'output/', RUN ,'/simgwas_', NCV, '-cv_', NLD ,'-LD-blocks_chr', CHR, '_causal-variants.txt'),
            sep = '\t', quote = F, col.names = T, row.names = F)
