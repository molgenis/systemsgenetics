

setwd("/groups/umcg-fg/tmp01/projects/downstreamer/depict2_bundle/summary_statistics/plaintext/pricePaper")


files <- list.files(pattern = ".gz")


file <- files[1]

for(file in files){
  summstats <- read.delim(gzfile(file))
  summstats$P <- 2*pnorm(-abs(summstats$Z))
  write.table(summstats, file = gzfile(paste0("pvalues/",file)), sep = "\t", quote = F, row.names = F)
}

names <- files
names <- sub("PASS_","", names)
names <- sub(".sumstats.gz","", names)
names <- sub("UKB_460K.[^_]+_","", names)


gwasList <- data.frame(names = names, rscol = 0, pcol = 6, files = paste0(getwd(),"/pvalues/",file))

write.table(gwasList, file = "pvalues/list.txt", sep = "\t", quote = F, row.names = F, col.names = F)



