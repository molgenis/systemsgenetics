#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)


remoter::client("localhost", port = 55506, password = "laberkak")



library(uwot)

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")


load(file = "Metadata/combinedMeta_2022_09_15.RData")

combinedMeta <- combinedMeta[!combinedMeta$exclude,]

samplesPasExpressionQc <- rownames(expPcs)

combinedMeta <- combinedMeta[rownames(combinedMeta) %in% samplesPasExpressionQc,]


sum()
tcgaAll <- combinedMeta[combinedMeta$study == "TCGA",]
colnames(tcgaAll)

tcgaCancer <- tcgaAll[tcgaAll$Cancer,]
colnames(tcgaCancer)

table(tcgaCancer$tcga.cgc_sample_sample_type)

tcgaCancer <- tcgaCancer[tcgaCancer$tcga.cgc_sample_sample_type == "Primary Tumor",]
table(tcgaCancer$tcga.gdc_cases.project.primary_site)

rpng(width = 1000, height = 1000)
par(mar = c(9,4,3,1))
barplot(table(tcgaCancer$tcga.gdc_cases.project.primary_site), las = 2, main = "TCGA primary tumors")
dev.off()

