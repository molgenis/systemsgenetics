#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)



remoter::client("localhost", port = 55503, password = "laberkak")

setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")


library(readr)


table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Metadata/metadata_sra1.txt", delim = "\t", quote = "", guess_max = 20000)
sraMeta1 <- as.data.frame(table_tmp[,-1])
rm(table_tmp)

table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Metadata/metadata_sra2.txt", delim = "\t", quote = "", guess_max = 20000)
sraMeta2 <- as.data.frame(table_tmp[,-1])
rm(table_tmp)

sraSharedCol <- intersect(colnames(sraMeta2), colnames(sraMeta1))
length(sraSharedCol)

sraMeta <- rbind(sraMeta1[,sraSharedCol], sraMeta2[,sraSharedCol])

#For some reason some runs are duplicated in the meta data file. 
#Quick inspection showed that they have the same values
#Solution exclude duplicate row
sraUniqueIds <- unique(sraMeta$external_id)
str(sraUniqueIds)
sraMeta <- sraMeta[ match(sraUniqueIds, sraMeta$external_id), ]
rownames(sraMeta) <- sraMeta$external_id



#extra columns in part 2
sraPart2Col <- colnames(sraMeta2)[!colnames(sraMeta2) %in% colnames(sraMeta1)]

sraMetaExtended <- sraMeta2[match(sraUniqueIds, sraMeta2$external_id),sraPart2Col]

str(sraMetaExtended)

sum(length(unique(sraMeta2$external_id)))
sum(length(unique(sraMeta1$external_id)))

sum(unique(length(sraMeta2$external_id)))
sum(unique(length(sraMeta1$external_id)))

#metadata_gtex
load(file = "/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Metadata/metadata_gtex.Rda", verbose = T)
rownames(metadata_gtex) <- metadata_gtex$external_id
metadata_gtex2 <- metadata_gtex[,c("gtex.smts", "gtex.smtsd")]
str(metadata_gtex2)

#metadata_tcga
load(file = "/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Metadata/metadata_tcga.Rda", verbose = T)
rownames(metadata_tcga) <- metadata_tcga$external_id
metadata_tcga2 <- metadata_tcga[,c("tcga.gdc_cases.project.primary_site", "tcga.cgc_sample_sample_type")]



#ARCH4 data
table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Metadata/metadataSRA.txt", delim = "\t", quote = "")
metadata_archs4 <- as.data.frame(table_tmp[,-1])
rownames(metadata_archs4) <- table_tmp[,3][[1]]
rm(table_tmp)
metadata_archs4_2 <- metadata_archs4[,c("Tissue", "CellType", "CellLine")]
colnames(metadata_archs4_2) <- c("archs4.Tissue", "archs4.CellType", "archs4.CellLine")

#GADO data
metadata_Gado <- read.delim("celllinesAndCancer/oldAnnotations/sampleAnnotations.txt")
rownames(metadata_Gado) <- metadata_Gado$Sample
metadata_Gado2 <- metadata_Gado[,c("CellLine", "TissueType", "CellType", "PlotClass")]
colnames(metadata_Gado2) <- paste0("gado.",colnames(metadata_Gado2))
gadoTissueCol <- read.delim("celllinesAndCancer/oldAnnotations/tissueCol5.txt")


#Kidney Network Annotaions
metadata_Kn <- read.delim("Metadata/KidneyNetwork.txt")
rownames(metadata_Kn) <- metadata_Kn$Sample
metadata_Kn2 <- metadata_Kn[,c("Origin", "Cell_type", "Cell_type_simplified", "Cell_type_manual")]
colnames(metadata_Kn2) <- paste0("KidneyNetwork.",colnames(metadata_Kn2))



allSamples <- c(rownames(metadata_gtex2), rownames(metadata_tcga2), rownames(sraMeta))

length(unique(allSamples)) == length(allSamples)

numberSamples = length(allSamples)
finalAnnotations <- data.frame(
  Tissue = rep("",numberSamples), 
  Tissue2 = rep("",numberSamples), 
  Cellline =  vector(mode = "logical", length = numberSamples), 
  CelllineName = rep("",numberSamples), 
  Cancer =  vector(mode = "logical", length = numberSamples),
  Cohort = rep("SRA",numberSamples),
  row.names = allSamples, stringsAsFactors = F)
finalAnnotations$Cellline = NA
finalAnnotations$Cancer = NA

finalAnnotations$Cohort[rownames(finalAnnotations) %in% rownames(metadata_gtex2)] <- "GTEx"
finalAnnotations$Cohort[rownames(finalAnnotations) %in% rownames(metadata_tcga2)] <- "TCGA"
table(finalAnnotations$Cohort, useNA = "always")
str(finalAnnotations)


dim(finalAnnotations)
dim(metadata_gtex2)


sum(rownames(finalAnnotations) %in% rownames(metadata_gtex2))

a <- merge(finalAnnotations, metadata_gtex2, all.x = T, by = 0)
row.names(a) <- a$Row.names

b <- merge(a, metadata_tcga2, all.x = T, by = 0)
row.names(b) <- b$Row.names

c <- merge(b, metadata_archs4_2, all.x = T, by = 0)
row.names(c) <- c$Row.names

d <- merge(c, metadata_Gado2, all.x = T, by = 0)
row.names(d) <- d$Row.names

e <- merge(d, sraMetaExtended, all.x = T, by = 0)
row.names(e) <- e$Row.names

f <- merge(e, sraMeta, all.x = T, by = 0)
row.names(f) <- f$Row.names

g <- merge(f, metadata_Kn2, all.x = T, by = 0)
row.names(g) <- g$Row.names

str(g)


combinedMeta <- g[,-c(1:7)]


head(rownames(combinedMeta))

combinedMeta$exclude <- FALSE
combinedMeta$fetal <- NA

#save(combinedMeta, file = "combinedMeta.RData")
#load(file = "combinedMeta.RData")

combinedMeta$Tissue[combinedMeta$Cohort == "GTEx"] <- combinedMeta$gtex.smts[combinedMeta$Cohort == "GTEx"]
combinedMeta$Tissue2[combinedMeta$Cohort == "GTEx"] <- combinedMeta$gtex.smtsd[combinedMeta$Cohort == "GTEx"]
combinedMeta$Cellline[combinedMeta$Cohort == "GTEx"] <- FALSE
combinedMeta$Cancer[combinedMeta$Cohort == "GTEx"] <- FALSE

gtexLcl <- combinedMeta$Cohort == "GTEx" & (!is.na(combinedMeta$gtex.smtsd) & combinedMeta$gtex.smtsd == "Cells - EBV-transformed lymphocytes")
combinedMeta$Cellline[gtexLcl] <- TRUE
combinedMeta$CelllineName[gtexLcl] <- "lcl"
combinedMeta$Tissue[gtexLcl] <- ""
combinedMeta$Tissue2[gtexLcl] <- ""

gtexCml <- combinedMeta$Cohort == "GTEx" & (!is.na(combinedMeta$gtex.smtsd) & combinedMeta$gtex.smtsd == "Cells - Leukemia cell line (CML)")
combinedMeta$Cellline[gtexCml] <- TRUE
combinedMeta$CelllineName[gtexCml] <- "cml"
combinedMeta$Tissue[gtexCml] <- ""
combinedMeta$Tissue2[gtexCml] <- ""


table(combinedMeta$gtex.smtsd[combinedMeta$Cohort == "GTEx"])

combinedMeta$Tissue[combinedMeta$Cohort == "TCGA"] <- combinedMeta$tcga.gdc_cases.project.primary_site[combinedMeta$Cohort == "TCGA"]
combinedMeta$Cellline[combinedMeta$Cohort == "TCGA"] <- FALSE
combinedMeta$Cancer[combinedMeta$Cohort == "TCGA"] <- TRUE #default for TCGA exception below
combinedMeta$Cancer[combinedMeta$Cohort == "TCGA" & combinedMeta$tcga.cgc_sample_sample_type == "Solid Tissue Normal"] <- FALSE


#Map GADO names to gtex names
combinedMeta$gado.TissueType[!is.na(combinedMeta$gado.TissueType) & combinedMeta$gado.TissueType != ""] <- gsub("(^[[:alpha:]])", "\\U\\1", combinedMeta$gado.TissueType[!is.na(combinedMeta$gado.TissueType) & combinedMeta$gado.TissueType != ""], perl=TRUE)#https://stackoverflow.com/questions/18509527/first-letter-to-upper-case
gadoTissueCol$PlotClass <- gsub("(^[[:alpha:]])", "\\U\\1", gadoTissueCol$PlotClass, perl=TRUE)#https://stackoverflow.com/questions/18509527/first-letter-to-upper-case
combinedMeta$gado.TissueType[!is.na(combinedMeta$gado.TissueType) & combinedMeta$gado.TissueType == "Adipose"] <- "Adipose Tissue"
gadoTissueCol$PlotClass[gadoTissueCol$PlotClass == "Adipose"] <- "Adipose Tissue"

#Fix
combinedMeta$gado.CellType[!is.na(combinedMeta$gado.CellType) & combinedMeta$gado.CellType == "acute myeloid leukemia"] <- "AML"

#Only annotations with a color are checked and highly realiable
gadoAnnotatedTissues <- !is.na(combinedMeta$gado.TissueType) & combinedMeta$gado.TissueType != "" & combinedMeta$gado.TissueType  %in% gadoTissueCol$PlotClass
combinedMeta$Tissue[gadoAnnotatedTissues] <- combinedMeta$gado.TissueType[gadoAnnotatedTissues]
combinedMeta$Cancer[gadoAnnotatedTissues] <- FALSE
combinedMeta$Cellline[gadoAnnotatedTissues] <- FALSE

gadoAnnotatedCelltypes <- !is.na(combinedMeta$gado.CellType) & combinedMeta$gado.CellType != "" & combinedMeta$gado.CellType  %in% gadoTissueCol$PlotClass
combinedMeta$Tissue2[gadoAnnotatedCelltypes] <- combinedMeta$gado.CellType[gadoAnnotatedCelltypes]
combinedMeta$Cancer[gadoAnnotatedCelltypes] <- FALSE
combinedMeta$Cellline[gadoAnnotatedCelltypes] <- FALSE

combinedMeta$Cancer[!is.na(combinedMeta$gado.CellType) & combinedMeta$gado.CellType == "AML"] <- TRUE
combinedMeta$Cancer[!is.na(combinedMeta$gado.CellType) & combinedMeta$gado.CellType == "DLBCL"] <- TRUE


gadoAnnotatedCelllines <- !is.na(combinedMeta$gado.CellLine) & combinedMeta$gado.CellLine != "" & tolower(combinedMeta$gado.CellLine)  %in% tolower(gadoTissueCol$PlotClass)
combinedMeta$CelllineName[gadoAnnotatedCelllines] <- combinedMeta$gado.CellLine[gadoAnnotatedCelllines]
combinedMeta$Cancer[gadoAnnotatedCelllines] <- FALSE
combinedMeta$Cellline[gadoAnnotatedCelllines] <- TRUE

#Some manual stuff for big studies

combinedMeta$CelllineName[!is.na(combinedMeta$study) & combinedMeta$study == "SRP166108"] <- "HepaRG"
combinedMeta$Cellline[!is.na(combinedMeta$study) & combinedMeta$study == "SRP166108"] <- TRUE

combinedMeta$exclude[!is.na(combinedMeta$study) & combinedMeta$study == "SRP025982"] <- TRUE


#TODO SRP013565 encode data search for table

#TODO ERP020977 We differentiated macrophages from induced pluripotent stem cells in 86 unrelated, healthy individuals derived by the Human Induced Pluripotent Stem Cells Initiative (HIPSCI),

combinedMeta$Tissue[!is.na(combinedMeta$study) & combinedMeta$study == "SRP192714"] <- "Blood"
combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "SRP192714"] <- "Whole Blood"

combinedMeta$Cellline[!is.na(combinedMeta$study) & combinedMeta$study == "SRP186687"] <- TRUE

combinedMeta$Tissue[!is.na(combinedMeta$study) & combinedMeta$study == "SRP092402"] <- "Blood"
combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "SRP092402"] <- "Whole Blood"

combinedMeta$Tissue[!is.na(combinedMeta$study) & combinedMeta$study == "SRP116272"] <- "Blood"
combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "SRP116272"] <- "Whole Blood"

combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "ERP007111"] <- "iPSC"

combinedMeta$Cellline[!is.na(combinedMeta$study) & combinedMeta$study == "ERP012914"] <- TRUE
combinedMeta$CelllineName[!is.na(combinedMeta$study) & combinedMeta$study == "ERP012914"] <- "HAP1"


combinedMeta$Tissue[!is.na(combinedMeta$study) & combinedMeta$study == "SRP187978"] <- "Liver" 

#TODO SRP074349 check for table. Combination of cancer and control of differnt tissues

combinedMeta$Tissue[!is.na(combinedMeta$study) & combinedMeta$study == "ERP016243"] <- "Brain" 
combinedMeta$fetal[!is.na(combinedMeta$study) & combinedMeta$study == "ERP016243"] <- TRUE 

combinedMeta$Cellline[!is.na(combinedMeta$study) & combinedMeta$study == "ERP020478"] <- TRUE
combinedMeta$CelllineName[!is.na(combinedMeta$study) & combinedMeta$study == "ERP020478"] <- "hela"

combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "SRP056295"] <- "AML"
combinedMeta$Cancer[!is.na(combinedMeta$study) & combinedMeta$study == "SRP056295"] <- TRUE

sum(rownames(combinedMeta)[!is.na(combinedMeta$study) & combinedMeta$study == "ERP107748"]%in% rownames(pcs))

combinedMeta$Tissue1[!is.na(combinedMeta$study) & combinedMeta$study == "SRP151763"] <- "Eye"
combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "SRP151763"] <- "Retina"


combinedMeta$Cellline[!is.na(combinedMeta$study) & combinedMeta$study == "SRP102077"] <- TRUE
combinedMeta$CelllineName[!is.na(combinedMeta$study) & combinedMeta$study == "SRP102077"] <- "H-STS"

combinedMeta$Tissue[!is.na(combinedMeta$study) & combinedMeta$study == "SRP162411"] <- "Blood"
combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "SRP162411"] <- "Whole Blood"



table(combinedMeta$Tissue)
table(combinedMeta$Tissue2)
table(combinedMeta$Cancer)

table(combinedMeta$Tissue, combinedMeta$Cancer)

table(combinedMeta$Tissue[combinedMeta$Cohort == "TCGA"], combinedMeta$Cancer[combinedMeta$Cohort == "TCGA"])


sum(combinedMeta$gado.TissueType %in% combinedMeta$gado.PlotClass)

sum((!is.na(combinedMeta$gado.TissueType) & combinedMeta$gado.TissueType != ""))

sum( )
unique(combinedMeta$gado.TissueType)

table(combinedMeta$gado.PlotClass, useNA = "a")

combinedMeta$Tissue[combinedMeta$Cohort == "GSA"]



sort(table(combinedMeta$study[rownames(combinedMeta) %in% rownames(pcs)]))
