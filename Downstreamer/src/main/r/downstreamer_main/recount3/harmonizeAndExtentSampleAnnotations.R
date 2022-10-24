#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)


remoter::client("localhost", port = 55501, password = "laberkak")


#save.image("tmp2.RData")


  
setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Recount3\\")
#load("tmp2.RData")


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


#Mahmoud annotations
load("Recount3_QC_2ndRun/SRA_Studies_Annotations_Patrick/Annotations.rda", verbose = T)
str(Annotations)
rownames(Annotations) <- Annotations$SampleID
mahmoudAnnotations <- Annotations[,-(1:2)]
#write.table(Annotations, sep = "\t", quote = F, col.names = NA, file = "tmp.txt")

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
finalAnnotations$Fetal <- NA

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
str(combinedMeta)

#now fillin the gtex and gcta recount meta data. 

tmp <- metadata_gtex[,colnames(metadata_gtex) %in% sraSharedCol]
combinedMeta[rownames(tmp),colnames(tmp)] <- tmp

tmp <- metadata_tcga[,colnames(metadata_tcga) %in% sraSharedCol]
combinedMeta[rownames(tmp),colnames(tmp)] <- tmp

rm(tmp)

#set study make column uniform
combinedMeta$study[combinedMeta$Cohort == "GTEx"] <- "GTEx"
combinedMeta$study[combinedMeta$Cohort == "TCGA"] <- "TCGA"

combinedMeta$exclude <- FALSE

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

gtexFibroblasts <- combinedMeta$Cohort == "GTEx" & (!is.na(combinedMeta$gtex.smtsd) & combinedMeta$gtex.smtsd == "Cells - Cultured fibroblasts")
combinedMeta$Cellline[gtexFibroblasts] <- TRUE
combinedMeta$CelllineName[gtexFibroblasts] <- "Fibroblasts"
combinedMeta$Tissue[gtexFibroblasts] <- ""
combinedMeta$Tissue2[gtexFibroblasts] <- ""


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


combinedMeta$Tissue[!is.na(combinedMeta$study) & combinedMeta$study == "SRP192714"] <- "Blood"
combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "SRP192714"] <- "Whole Blood"
combinedMeta$Cellline[!is.na(combinedMeta$study) & combinedMeta$study == "SRP192714"] <- FALSE
combinedMeta$Cancer[!is.na(combinedMeta$study) & combinedMeta$study == "SRP192714"] <- FALSE

combinedMeta$Cellline[!is.na(combinedMeta$study) & combinedMeta$study == "SRP186687"] <- TRUE

combinedMeta$Tissue[!is.na(combinedMeta$study) & combinedMeta$study == "SRP092402"] <- "Blood"
combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "SRP092402"] <- "Whole Blood"
combinedMeta$Cellline[!is.na(combinedMeta$study) & combinedMeta$study == "SRP092402"] <- FALSE
combinedMeta$Cancer[!is.na(combinedMeta$study) & combinedMeta$study == "SRP092402"] <- FALSE


combinedMeta$Tissue[combinedMeta$study == "SRP116272" & grepl("source_name;;T cells", combinedMeta$sra.sample_attributes)] <- "Blood"
combinedMeta$Tissue2[combinedMeta$study == "SRP116272" & grepl("source_name;;T cells", combinedMeta$sra.sample_attributes)] <- "T-cells"
combinedMeta$Cellline[combinedMeta$study == "SRP116272" & grepl("source_name;;T cells", combinedMeta$sra.sample_attributes)] <- FALSE
combinedMeta$Cancer[combinedMeta$study == "SRP116272" & grepl("source_name;;T cells", combinedMeta$sra.sample_attributes)] <- FALSE

combinedMeta$Tissue[combinedMeta$study == "SRP116272" & grepl("source_name;;Monocytes", combinedMeta$sra.sample_attributes)] <- "Blood"
combinedMeta$Tissue2[combinedMeta$study == "SRP116272" & grepl("source_name;;Monocytes", combinedMeta$sra.sample_attributes)] <- "Monocytes"
combinedMeta$Cellline[combinedMeta$study == "SRP116272" & grepl("source_name;;Monocytes", combinedMeta$sra.sample_attributes)] <- FALSE
combinedMeta$Cancer[combinedMeta$study == "SRP116272" & grepl("source_name;;Monocytes", combinedMeta$sra.sample_attributes)] <- FALSE

combinedMeta$Tissue[combinedMeta$study == "SRP061932"] <- ""
combinedMeta$Tissue2[combinedMeta$study == "SRP061932"] <- ""
combinedMeta$Cellline[combinedMeta$study == "SRP061932"] <- FALSE
combinedMeta$Cancer[combinedMeta$study == "SRP061932"] <- FALSE



combinedMeta$Tissue[combinedMeta$study == "SRP047323"] <- ""
combinedMeta$Tissue2[combinedMeta$study == "SRP047323"] <- ""
combinedMeta$Cellline[combinedMeta$study == "SRP047323"] <- FALSE
combinedMeta$Cancer[combinedMeta$study == "SRP047323"] <- FALSE



combinedMeta$CelllineName[combinedMeta$study == "ERP001942"] <- "lcl"
combinedMeta$Cellline[combinedMeta$study == "ERP001942"] <- TRUE


combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "ERP007111"] <- "iPSC"

combinedMeta$Cellline[!is.na(combinedMeta$study) & combinedMeta$study == "ERP012914"] <- TRUE
combinedMeta$CelllineName[!is.na(combinedMeta$study) & combinedMeta$study == "ERP012914"] <- "HAP1"


combinedMeta$Tissue[!is.na(combinedMeta$study) & combinedMeta$study == "SRP151763"] <- "Eye"
combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "SRP151763"] <- "Retina"
combinedMeta$Cellline[combinedMeta$study == "SRP151763"] <- FALSE
combinedMeta$Cancer[combinedMeta$study == "SRP151763"] <- FALSE



combinedMeta$Tissue[!is.na(combinedMeta$study) & combinedMeta$study == "SRP162411"] <- "Blood"
combinedMeta$Tissue2[!is.na(combinedMeta$study) & combinedMeta$study == "SRP162411"] <- "Whole Blood"

studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP102542"
combinedMeta$Tissue[studySamples] <- "Muscle"
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$Cancer[studySamples] <- FALSE


studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP150311"
combinedMeta$Tissue[studySamples] <- "Muscle"
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$Cancer[studySamples] <- FALSE


studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP162873"
combinedMeta$Tissue[studySamples] <- "Muscle"
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$Cancer[studySamples] <- FALSE



studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP163524"
combinedMeta$Tissue[studySamples] <- "Muscle"
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$Cancer[studySamples] <- FALSE


studySamples <- !is.na(combinedMeta$study) & combinedMeta$study %in% c("SRP006676", "SRP071758", "SRP081599", "SRP086078", "SRP119923")
combinedMeta$Tissue[studySamples] <- "Airway Epithelial"
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$Cancer[studySamples] <- FALSE


samples <- combinedMeta$study == "SRP188219" & grepl("left atrial appendage", combinedMeta$sra.sample_attributes)
combinedMeta$Tissue[samples] <- "Heart"
combinedMeta$Tissue2[samples] <- "Left atrial appendage"
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- FALSE


samples <- combinedMeta$study == "SRP188219" & grepl("right atrial appendage", combinedMeta$sra.sample_attributes)
combinedMeta$Tissue[samples] <- "Heart"
combinedMeta$Tissue2[samples] <- "Right atrial appendage"
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- FALSE



studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRA755613"
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- TRUE
combinedMeta$CelllineName[studySamples] <- "iPSC"
combinedMeta$Cancer[studySamples] <- FALSE


studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRA755626"
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- TRUE
combinedMeta$CelllineName[studySamples] <- "iPSC"
combinedMeta$Cancer[studySamples] <- FALSE




studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP148659"
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- TRUE
combinedMeta$CelllineName[studySamples] <- "iPSC"
combinedMeta$Cancer[studySamples] <- FALSE

#Put to missing annotations unclear
studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP009316"
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- NA
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- NA


studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP021509"
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- NA
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- NA


#Unsure how to classify airway smooth muscle
studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP043162"
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- NA
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- NA

studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP052896"
combinedMeta$Cancer[studySamples] <- TRUE




#Some organoids and cancers
studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP058722"
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- NA
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- NA




#DPN and Tamoxifen treatments of parathyroid adenoma cells have cancer CNV profile
studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP012167"
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- NA
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- NA



studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP019936"
combinedMeta$Cancer[studySamples] <- TRUE



combinedMeta["SRR5341594", "sra.sample_title"] <- "Human differentiating macrophage"



studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "ERP011411"
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- NA
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- NA


studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP010166"
combinedMeta$Cancer[studySamples] <- TRUE




studySamples <- !is.na(combinedMeta$study) & combinedMeta$study == "SRP012656"
combinedMeta$Cancer[studySamples] <- TRUE


samples <- combinedMeta$study == "ERP006077" & grepl("Primary Prostate Tumour", combinedMeta$sra.sample_attributes)
combinedMeta$Tissue[samples] <- "Prostate"
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE

samples <- combinedMeta$study == "ERP006077" & grepl("Matched Adjacent", combinedMeta$sra.sample_attributes)
combinedMeta$Tissue[samples] <- "Prostate"
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- FALSE

samples <- combinedMeta$study == "ERP006077"
combinedMeta$Tissue[samples] <- "Pancreas"
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE

samples <- combinedMeta$study == "SRP058587"
combinedMeta$Tissue[samples] <- "Breast"
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE

samples <- combinedMeta$study == "SRP062332"
combinedMeta$Tissue[samples] <- ""
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- TRUE
combinedMeta$Cancer[samples] <- NA



samples <- combinedMeta$study == "SRP030401"
combinedMeta$Tissue[samples] <- "Breast"
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE



samples <- combinedMeta$study == "SRP028344"
combinedMeta$Tissue[samples] <- ""
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- TRUE
combinedMeta$Cancer[samples] <- NA




samples <- combinedMeta$study == "SRP073061" & grepl("Tumor", combinedMeta$sra.experiment_title)
combinedMeta$Tissue[samples] <- "Breast"
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE




samples <- combinedMeta$study == "SRP028346"
combinedMeta$Tissue[samples] <- ""
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- TRUE
combinedMeta$Cancer[samples] <- NA


samples <- combinedMeta$study == "SRP058571"
combinedMeta$Tissue[samples] <- ""
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- TRUE
combinedMeta$Cancer[samples] <- NA

studySamples <- combinedMeta$study %in% c("SRP014027", "SRP006575", "SRP071932", "ERP004617", "SRP034592","SRP049695")
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- NA
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- NA



samples <- combinedMeta$study == "SRP066596"
combinedMeta$Tissue[samples] <- ""
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- NA
combinedMeta$Cancer[samples] <- TRUE



samples <- combinedMeta$study == "SRP049648"
combinedMeta$Tissue[samples] <- ""
combinedMeta$Tissue2[samples] <- ""
combinedMeta$Cellline[samples] <- TRUE
combinedMeta$Cancer[samples] <- NA




studySamples <- combinedMeta$study == "SRP039694"
combinedMeta$Cancer[studySamples] <- TRUE


studySamples <- combinedMeta$study == "SRP066260"
combinedMeta$Cancer[studySamples] <- TRUE



#mahmoudAnnotations

colnames(mahmoudAnnotations)[colnames(mahmoudAnnotations) == "Cell_Line"] <- "Cellline"
colnames(mahmoudAnnotations)[colnames(mahmoudAnnotations) == "Cell_Line_Name"] <- "CelllineName"

all(colnames(mahmoudAnnotations) %in% colnames(combinedMeta))
all(rownames(mahmoudAnnotations) %in% rownames(combinedMeta))


combinedMeta[rownames(mahmoudAnnotations),colnames(mahmoudAnnotations)] <- mahmoudAnnotations


#All cellline to false for all tissues
combinedMeta$Cancer[combinedMeta$Cellline] <- FALSE

tmp <- !is.na(combinedMeta$Tissue) & combinedMeta$Tissue == "Cervix Uteri"
combinedMeta$Tissue[tmp] <- "Uterus"
combinedMeta$Tissue2[tmp] <- "Cervix"

tmp <- !is.na(combinedMeta$Tissue) & combinedMeta$Tissue == "Cervix"
combinedMeta$Tissue[tmp] <- "Uterus"
combinedMeta$Tissue2[tmp] <- "Cervix"

tmp <- !is.na(combinedMeta$Tissue) & combinedMeta$Tissue == "Lymph node"
combinedMeta$Tissue[tmp] <- "Lymph Nodes"

tmp <- !is.na(combinedMeta$Tissue) & combinedMeta$Tissue == "Bone marrow"
combinedMeta$Tissue[tmp] <- "Bone Marrow"

tmp <- !is.na(combinedMeta$Tissue) & combinedMeta$Tissue == "Colorectal"
combinedMeta$Tissue[tmp] <- "Colon"


tmp <- !is.na(combinedMeta$Tissue) & combinedMeta$Tissue == "Whole blood"
combinedMeta$Tissue2[tmp] <- "Whole Blood"
combinedMeta$Tissue[tmp] <- "Blood"



#Below are tissues2 fixes by Mahmoud
# annotations already present in Tissue are removed from Tissue2
#duplicated are harmonized
#set rare annotations to NA
### All parts of the basal ganglia (including substantia nigra) were annotated as basal ganglia
###brain fragements was set to NA
###Retina needs to have Eye as Tissue
###Sample annotated as both brain & stomach was annotated as NA

#Adipose Tissue 
# Tissue2 includes "Adipose - Subcutaneous" & "Adipose - Visceral (Omentum)"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Adipose - Subcutaneous"),2]= "Subcutaneous"
combinedMeta[!is.na(combinedMeta$Tissue2) & combinedMeta$Tissue2== "Adipose - Visceral (Omentum)",2]= "Visceral"
combinedMeta[!is.na(combinedMeta$Tissue) & combinedMeta$Tissue== "Adipose Tissue",1]= "Adipose"

# Adrena Gland
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Adrenal Gland"),2]= NA

#AML
# keep as is

# Arteries
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Artery - Aorta"),2]= "Aorta"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Artery - Coronary"),2]= "Coronary"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Artery - Tibial"),2]= "Tibial"

#B-cells
# keep as is

#basal ganglion
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "basal ganglion"),2]= "Basal Ganglia"

#Bladder
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Bladder"),2]= NA

#Brain (keep as GTEX)****
#Check for brain cortex vs cortex vs cerebral cortex
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Amygdala"),2]= "Amygdala"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Anterior cingulate cortex (BA24)"),2]= "Anterior cingulate cortex"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Caudate (basal ganglia)"),2]= "Caudate (basal ganglia)"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Cerebellar Hemisphere"),2]= "Cerebellar Hemisphere"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Cerebellum"),2]= "Cerebellum"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Cortex"),2]= "Cortex"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Frontal Cortex (BA9)"),2]= "Frontal Cortex"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Hippocampus"),2]= "Hippocampus"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Hypothalamus"),2]= "Hypothalamus"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Nucleus accumbens (basal ganglia)"),2]= "Nucleus accumbens (basal ganglia)"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Putamen (basal ganglia)"),2]= "Putamen (basal ganglia)"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Spinal cord (cervical c-1)"),2]= "Spinal Cord (cervical c-1)"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Brain - Substantia nigra"),2]= "Substantia nigra"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "brain fragment"),2]= NA

#Breast
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Breast - Mammary Tissue"),2]= "Mammary Tissue"

#CD34+
# Keep as is

#Cultured fibroblasts
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Cells - Cultured fibroblasts"),2]= "Cultured Fibroblasts"

#cerebellum
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "cerebellum"),2]= "Cerebellum"

# cerebral cortex 
#recheck****
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "cerebral cortex"),2]= "Cortex"

#Cervix
#keep as is

#choroid plexus
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "choroid plexus"),2]= NA

#Colon
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Colon - Sigmoid"),2]= "Sigmoid"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Colon - Transverse"),2]= "Transverse"

#diencephalon & diencephalon and midbrain
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "diencephalon"),2]= "Diencephalon"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "diencephalon and midbrain"),2]= NA

#DLBCL
#keep as is 

#Esophagus
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Esophagus - Gastroesophageal Junction"),2]= "Gastroesophageal Junction"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Esophagus - Mucosa"),2]= "Mucosa"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Esophagus - Muscularis"),2]= "Muscularis"

#Fallopian Tube
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Fallopian Tube"),2]= NA

#forebrain
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "forebrain"),2]= NA
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "forebrain and midbrain"),2]= NA
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "forebrain fragment"),2]= NA

#Heart
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Heart - Atrial Appendage"),2]= "Atrial Appendage"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Heart - Left Ventricle"),2]= "Left Ventricle"

#hindbrain
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "hindbrain"),2]= NA
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "hindbrain fragment"),2]= NA
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "hindbrain without cerebellum"),2]= NA
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "hippocampus"),2]= "Hippocampus"

#Kidney
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Kidney - Cortex"),2]= "Cortex"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Kidney - Medulla"),2]= "Medulla"

#Liver
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Liver"),2]= NA

#Lung
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Lung"),2]= NA

#medulla oblongata
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "medulla oblongata"),2]= "Medulla Oblongata"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "midbrain"),2]= "Midbrain"

#Minor salivary gland
#keep as is

#Monocytes
#keep as is

#Muscle-skeletal
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Muscle - Skeletal"),2]= "Skeletal"

#Nerve
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Nerve - Tibial"),2]= "Tibial"

#NK-cells
#keep as is

#Ovary
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Ovary"),2]= NA

#Pancreas
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Pancreas"),2]= NA

#PBMC
#keep as is

#Pituitary
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Pituitary"),2]= NA

#pituitary and diencephalon & pons
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "pituitary and diencephalon"),2]= NA
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "pons"),2]= NA

#prostate
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Prostate"),2]= NA

#skin
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Skin - Not Sun Exposed (Suprapubic)"),2]= "Suprapubic"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Skin - Sun Exposed (Lower leg)"),2]= "Lower Leg"

# Small Intesine
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Small Intestine - Terminal Ileum"),2]= "Terminal Ileum"

#spinal cord
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "spinal cord"),2]= ""

#Spleen 
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Spleen"),2]= NA

#Stomach 
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Stomach"),2]= NA
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "stomach"),2]= NA
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "stomach"),1]= NA

#T-cells
#Keep as is

#telencephalon
#too general
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "telencephalon"),2]= NA

#temporal lobe
#The temporal lobe is one of the four major lobes of the cerebral cortex
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "temporal lobe"),2]= "Cortex"

#Testis
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Testis"),2]=NA

#Thyroid
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Thyroid"),2]=NA

#Uterus
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Uterus"),2]=NA

#Vagina
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "Vagina"),2]=NA

#whole blood
#keep as is

#remove iPSCs from 
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "iPSC") ,4]="iPSC"
combinedMeta[!is.na(combinedMeta$Tissue2) & (combinedMeta$Tissue2== "iPSC") ,2]=""

#remove NA problem
combinedMeta$Tissue2[is.na(combinedMeta$Tissue2)]<- ""
combinedMeta$Tissue[is.na(combinedMeta$Tissue)]<- ""







#Harmonizing Cell Line Names for samples in recount3

# A549
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "a549"),4]= "A549"

#H-STS NET 
#Keep as is

#H1
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "h1"),4]= "H1"

#H9
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "h9"),4]= "H9"

#HAP1
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "hap1"),4]= "HAP1"

#HCT116
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "hct116"),4]= "HCT116"

#Hek293
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "hek293"),4]= "HEK293"

#HeLa
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "hela"),4]= "HeLa"

#HepaRG
#keep as is

#hepg2
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "hepg2"),4]= "HepG2"

#ipsc
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "ipsc"),4]= "iPSC"

#K562
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "k562"),4]= "K562"

#LCLs
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "lcl"),4]= "LCL"

#lcl_s4u_capturing
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "lcl_s4u_capturing"),4]= "LCL_S4U_Capturing"

#MCF10A
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "mcf10a"),4]= "MCF10A"

#MCF7
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "mcf7"),4]= "MCF7"

#MDA231
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "mda231"),4]= "MDA231"

#T47D
combinedMeta[!is.na(combinedMeta$CelllineName) & (combinedMeta$CelllineName== "t47d"),4]= "T47D"


#Fix SRP045234 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045234"),1]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045234"),2]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045234"),3]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045234"),4]= "iPSC"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045234"),5]= FALSE
#Fix SRP007525 Annotaions
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP007525"),1]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP007525"),2]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP007525"),3]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP007525"),4]= "OCI-LY1"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP007525"),5]= FALSE
#Fix SRP027358 & SRP032926
combinedMeta$Fetal[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP027358")]= TRUE
combinedMeta$Fetal[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP032926")]= TRUE

#Fix SRP026537 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP026537"),1]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP026537"),2]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP026537"),3]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP026537"),4]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP026537"),5]= FALSE

#Fix SRP049063 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP049063"),1]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP049063"),2]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP049063"),3]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP049063"),4]= "HT-29"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP049063"),5]= FALSE
#Fix SRP053034 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP053034"),1]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP053034"),2]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP053034"),3]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP053034"),4]= "RPE-1"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP053034"),5]= FALSE
#Fix SRP056197 Annoations
samples <- combinedMeta$study == "SRP056197" & grepl("Bone marrow", combinedMeta$sra.sample_attributes)
combinedMeta$Tissue[samples] <- "Bone Marrow"
combinedMeta$Tissue2[samples] <- "AML"
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE
samples <- combinedMeta$study == "SRP056197" & grepl("Heparinised blood", combinedMeta$sra.sample_attributes)
combinedMeta$Tissue[samples] <- "Blood"
combinedMeta$Tissue2[samples] <- "AML"
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE


#Fix SRP013565 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP013565"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP013565"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP013565"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP013565"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP013565"),"Cancer"]= FALSE

#Fix ERP008682 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP008682"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP008682"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP008682"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP008682"),"CelllineName"]= "H9"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP008682"),"Cancer"]= FALSE
#Fix SRP033646 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP033646"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP033646"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP033646"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP033646"),"CelllineName"]= "Caco2"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP033646"),"Cancer"]= FALSE
#Fix SRP027383 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP027383"),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP027383"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP027383"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP027383"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP027383"),"Cancer"]= TRUE
#Fix SRP050003 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050003"),"Tissue"]= "Liver"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050003"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050003"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050003"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050003") & (grepl("non-tumoral", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050003") & (grepl("carcinoma", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= TRUE
#Fix SRP073253 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP073253"),"Tissue"]= "Kidney"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP073253"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP073253"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP073253"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP073253"),"Cancer"]= TRUE
#Fix SRP069235 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP069235"),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP069235"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP069235"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP069235"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP069235"),"Cancer"]= TRUE
#Fix SRP074425 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074425"),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074425"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074425"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074425"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074425"),"Cancer"]= TRUE
#Fix SRP044668 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP044668"),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP044668"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP044668"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP044668"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP044668") & (grepl("non-neoplastic", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP044668") & (grepl("glioma", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= TRUE
#Fix SRP009123 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009123"),"Tissue"]= "Liver"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009123"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009123"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009123"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009123") & (grepl("non-tumor", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009123") & (grepl("carcinoma", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= TRUE
#Fix SRP041094 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP041094"),"Tissue"]= "Prostate"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP041094"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP041094"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP041094"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP041094"),"Cancer"]= TRUE
#Fix SRP040998 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP040998"),"Tissue"]= "Liver"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP040998"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP040998"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP040998"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP040998"),"Cancer"]= NA
#Fix SRP052056 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP052056"),"Tissue"]= "Thyroid"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP052056"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP052056"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP052056"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP052056") & (grepl("healthy", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP052056") & (grepl("carcinoma", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= TRUE
#Fix SRP029880 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP029880"),"Tissue"]= "Colon"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP029880"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP029880"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP029880"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP029880"),"Cancer"]= TRUE
#Fix SRP056696 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP056696"),"Tissue"]= "Liver"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP056696"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP056696"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP056696"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP056696") & (grepl("Normal", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP056696") & (grepl("Tumor", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= TRUE
#Fix SRP066794 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP066794"),"Tissue"]= "Lung"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP066794"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP066794"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP066794"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP066794"),"Cancer"]= TRUE
#Fix SRP149374 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP149374"),"Tissue"]= "Bone Marrow"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP149374"),"Tissue2"]= "CD34+"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP149374"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP149374"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP149374") & (grepl("Healthy", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP149374") & (grepl("Myelodysplastic", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= TRUE
#Fix SRP019250 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019250"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019250"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019250"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019250") & (grepl("HEK", combinedMeta$sra.sample_attributes, ignore.case=T)),"CelllineName"]= "HEK"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019250") & (grepl("LCL", combinedMeta$sra.sample_attributes, ignore.case=T)),"CelllineName"]= "LCL"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019250"),"Cancer"]= FALSE
#Fix SRP074349 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074349"),"Tissue"]= "Lung"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074349"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074349"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074349"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074349") & (grepl("NSCLC", combinedMeta$sra.sample_attributes, ignore.case=T))==F,"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP074349") & (grepl("NSCLC", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= TRUE


#Fix SRP009067 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009067"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009067"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009067"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009067"),"CelllineName"]= "LCL"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009067"),"Cancer"]= FALSE
#Fix SRP007885 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP007885"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP007885"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP007885"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP007885"),"CelllineName"]= "LCL"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP007885"),"Cancer"]= FALSE
#Fix SRP018218 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP018218") & (grepl("cell line", combinedMeta$sra.sample_attributes, ignore.case=T))==F,"Tissue"]= "Pancreas"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP018218") & (grepl("cell line", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP018218") & (grepl("cell line", combinedMeta$sra.sample_attributes, ignore.case=T))==F,"Tissue2"]= "Stellate Cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP018218") & (grepl("cell line", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP018218") & (grepl("cell line", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP018218") & (grepl("cell line", combinedMeta$sra.sample_attributes, ignore.case=T))==F,"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP018218"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP018218"),"Cancer"]= TRUE
#Fix SRP019275 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019275"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019275"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019275"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019275"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP019275"),"Cancer"]= FALSE
#Fix SRP042186 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042186"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042186"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042186"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042186"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042186"),"Cancer"]= FALSE
#Fix SRP042620 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042620") & (grepl("cell line", combinedMeta$sra.sample_attributes, ignore.case=T))==F,"Tissue"]= "Breast"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042620") & (grepl("cell line", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042620"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042620") & (grepl("cell line", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042620") & (grepl("cell line", combinedMeta$sra.sample_attributes, ignore.case=T))==F,"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042620"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042620") & (grepl("ER+", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042620") & (grepl("Triple Negative", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042620") & (grepl("Uninvolved", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP042620") & (grepl("No Known", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= FALSE


#Fix ERP010142 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP010142"),"Tissue"]= "Breast"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP010142"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP010142"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP010142"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP010142"),"Cancer"]= TRUE
#Fix SRP026600 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP026600"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP026600"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP026600"),"Cellline"]= NA
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP026600"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP026600"),"Cancer"]= NA
#Fix SRP028336 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336") & (grepl("muscle", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Muscle"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336") & (grepl("kidney", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Kidney"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336") & (grepl("prefrontal cortex", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336") & (grepl("cerebellar Cortex", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336") & (grepl("primary visual cortex", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336") & (grepl("muscle", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336") & (grepl("kidney", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336") & (grepl("prefrontal cortex", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Prefrontal Cortex"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336") & (grepl("cerebellar Cortex", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Cerebellar Cortex"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336") & (grepl("primary visual cortex", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Primary Visual Cortex"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP028336"),"Cancer"]= FALSE
#Fix SRP009029 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009029"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009029"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009029"),"Cellline"]= NA
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009029"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009029"),"Cancer"]= NA
#Fix SRP006912 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP006912"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP006912"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP006912"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP006912"),"CelllineName"]= "HK-2"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP006912"),"Cancer"]= FALSE
#Fix SRP055444 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055444"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055444"),"Tissue2"]= "CLL"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055444"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055444"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055444"),"Cancer"]= TRUE
#Fix SRP022942 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP022942"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP022942"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP022942"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP022942"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP022942"),"Cancer"]= FALSE
#Fix ERP012180 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP012180"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP012180"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP012180"),"Cellline"]= NA
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP012180"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP012180"),"Cancer"]= NA
#Fix SRP058717 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP058717"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP058717"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP058717"),"Cellline"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP058717"),"CelllineName"]= "HT-29"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP058717"),"Cancer"]= FALSE
#Fix SRP012568 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP012568"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP012568"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP012568"),"Cellline"]= NA
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP012568"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP012568"),"Cancer"]= NA
#Fix SRP065146 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP065146"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP065146"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP065146"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP065146"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP065146"),"Cancer"]= FALSE
#Fix ERP012188 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP012188"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP012188"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP012188"),"Cellline"]= NA
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP012188"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP012188"),"Cancer"]= NA
#Fix SRP055390 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055390"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055390"),"Tissue2"]= "CLL"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055390"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055390"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055390") & (grepl("normal", combinedMeta$sra.sample_attributes, ignore.case=T))==F,"Cancer"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055390") & (grepl("normal", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055390") & (grepl("normal", combinedMeta$sra.sample_attributes, ignore.case=T))==F,"Tissue2"]="CLL"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055390") & (grepl("normal", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "B-cells"


#Fix SRP036145 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP036145"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP036145"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP036145"),"Cellline"]= NA
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP036145"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP036145"),"Cancer"]= NA
#Fix SRP050533 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050533"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050533"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050533"),"Cellline"]= NA
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050533"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050533"),"Cancer"]= NA






#Fix ERP016243
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"Fetal"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"Tissue2"]= ""

table(paste0(combinedMeta$Tissue, " - ", combinedMeta$Tissue2, " - ", combinedMeta$Fetal)[combinedMeta$study == "ERP109002"])

#Fix SRP078234 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078234") & (grepl("hindbrain", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078234") & (grepl("hindbrain", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Hindbrain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078234") & (grepl("spinal cord", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078234") & (grepl("spinal cord", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Spinal Cord (cervical c-1)"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078234"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078234"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078234"),"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078234"),"Fetal"]= TRUE


#Fix SRP041044 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP041044"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP041044"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP041044"),"Cellline"]= NA
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP041044"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP041044"),"Cancer"]= NA
#Fix SRP050260 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050260"),"Tissue"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050260"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050260"),"Cellline"]= NA
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050260"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP050260"),"Cancer"]= NA
#Fix SRP076099 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076099"),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076099"),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076099"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076099"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076099"),"Cancer"]= NA
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076099"),"Fetal"]= TRUE
# Fix ERP115010 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP115010"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP115010"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP115010"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP115010"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP115010"),"Cancer"]= FALSE
# Fix SRP221482 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP221482"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP221482"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP221482"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP221482"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP221482"),"Cancer"]= FALSE
# Fix SRP059039 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP059039"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP059039"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP059039"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP059039"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP059039"),"Cancer"]= FALSE
# Fix SRP059172 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP059172"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP059172"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP059172"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP059172"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP059172"),"Cancer"]= FALSE
# Fix SRP062966 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP062966"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP062966"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP062966"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP062966"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP062966"),"Cancer"]= FALSE
# Fix SRP081605 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP081605"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP081605"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP081605"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP081605"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP081605"),"Cancer"]= FALSE
# Fix SRP103772 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP103772"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP103772"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP103772"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP103772"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP103772"),"Cancer"]= FALSE
# Fix SRP132939 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP132939"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP132939"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP132939"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP132939"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP132939"),"Cancer"]= FALSE
# Fix SRP136938 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP136938"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP136938"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP136938"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP136938"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP136938"),"Cancer"]= FALSE
# Fix SRP150552 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"Cancer"]= FALSE
# Fix SRP174223 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP174223"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP174223"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP174223"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP174223"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP174223"),"Cancer"]= FALSE
# Fix SRP174638 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP174638"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP174638"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP174638"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP174638"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP174638"),"Cancer"]= FALSE
# Fix SRP150552 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP150552"),"Cancer"]= FALSE




samples <- combinedMeta$study == "SRP033266" & grepl("Bone marrow", combinedMeta$sra.sample_attributes)
combinedMeta$Tissue[samples] <- "Bone Marrow"
combinedMeta$Tissue2[samples] <- "AML"
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE
samples <- combinedMeta$study == "SRP033266" & grepl("Heparinised blood", combinedMeta$sra.sample_attributes)
combinedMeta$Tissue[samples] <- "Blood"
combinedMeta$Tissue2[samples] <- "AML"
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE


samples <- combinedMeta$study == "SRP048759" & grepl("Bone marrow", combinedMeta$sra.sample_attributes)
combinedMeta$Tissue[samples] <- "Bone Marrow"
combinedMeta$Tissue2[samples] <- "AML"
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE
samples <- combinedMeta$study == "SRP048759" & grepl("Heparinised blood", combinedMeta$sra.sample_attributes)
combinedMeta$Tissue[samples] <- "Blood"
combinedMeta$Tissue2[samples] <- "AML"
combinedMeta$Cellline[samples] <- FALSE
combinedMeta$Cancer[samples] <- TRUE

#Fix SRP045500 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500") & (grepl("wholw blood", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500") & (grepl("CD4", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "T-cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500") & (grepl("CD8", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "T-cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500") & (grepl("Neutrophils", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Neutrophils"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500") & (grepl("Monocytes", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Monocytes"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500") & (grepl("B-cells", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "B-cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500") & (grepl("NK", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "NK-cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP045500"),"Cancer"]= FALSE

#Fix SRP076719 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076719") & (grepl("pbmc", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076719") & (grepl("ln", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Blood"#in other cases we also put all t-cell to blood regardless of source
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076719"),"Tissue2"]= "T-cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076719"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076719"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP076719"),"Cancer"]= FALSE

#Fix SRP051688 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688") & (grepl("T cells", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "T-cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688") & (grepl("B cells", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "B-cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688") & (grepl("monocytes", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Monocytes"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688") & (grepl("NK cells", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "NK-cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688") & (grepl("PBMC", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "PBMC"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688") & (grepl("myeloid DC", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Dendritic cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688") & (grepl("neutrophils", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Neutrophils"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP051688"),"Cancer"]= FALSE

#Fix SRP078912 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078912"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078912") & (grepl("T cells", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "T-cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078912") & (grepl("Monocyte", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Monocytes"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078912"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078912"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP078912"),"Cancer"]= FALSE

#Fix SRP110609 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP110609"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP110609") & (grepl("lymphocytes", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP110609") & (grepl("Monocyte", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Monocytes"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP110609"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP110609"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP110609"),"Cancer"]= FALSE

#Fix SRP158943 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP158943"),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP158943") & (grepl("cll", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "CLL" #It doesn't state further clasification 
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP158943") & (grepl("B cells", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "B-cells"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP158943"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP158943"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP158943") & (grepl("cll", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= TRUE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP158943") & (grepl("B cells", combinedMeta$sra.sample_attributes, ignore.case=T)),"Cancer"]= FALSE

#Fix ERP104864 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP104864") & (grepl("blood", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP104864") & (grepl("synovium", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Synovium"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP104864") & (grepl("blood", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Whole Blood"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP104864") & (grepl("synovium", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP104864"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP104864"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP104864"),"Cancer"]= FALSE

#Fix intestine samples
#ERP000546
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP000546") & (combinedMeta$Tissue=="Intestine"),"Tissue"]= "Colon"
#ERP003613
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP003613") & (combinedMeta$Tissue=="Intestine") & (grepl("colon", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Colon"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP003613") & (combinedMeta$Tissue=="Intestine") & (grepl("smallintestine", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Small Intestine"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP003613") & (combinedMeta$Tissue=="Intestine") & (grepl("duodenum", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Small Intestine"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP003613") & (combinedMeta$Tissue=="Intestine") & (grepl("duodenum", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Duodenum"
#ERP006650
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP006650") & (combinedMeta$Tissue=="Intestine"),"Tissue"]= "Colon"
#SRP039090
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP039090") & (combinedMeta$Tissue=="Intestine") & (grepl("Small Intestine", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Small Intestine"
#SRP043391
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP043391") & (combinedMeta$Tissue=="Intestine") & (grepl("colon", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Colon"
#SRP048801
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP048801") & (combinedMeta$Tissue=="Intestine") & (grepl("ileum", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Small Intestine"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP048801") & (combinedMeta$Tissue=="Intestine") & (grepl("ileum", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= "Ileum"
#SRP055438
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP055438") & (combinedMeta$Tissue=="Intestine"),"Tissue"]= "Colon"
#SRP056520
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP056520") & (combinedMeta$Tissue=="Intestine"),"Tissue"]= "Colon"
#SRP006900
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP006900") & (combinedMeta$Tissue=="Intestine"),"Tissue"]= "Colon"
#SRP063496
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP063496") & (combinedMeta$Tissue=="Intestine"),"Tissue"]= "Colon"
#SRP000941
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP000941") & (combinedMeta$Tissue=="Intestine") & (grepl("colon", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Colon"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP000941") & (combinedMeta$Tissue=="Intestine") & (grepl("small intestine", combinedMeta$sra.sample_attributes, ignore.case=T)),"Tissue"]= "Small Intestine"
#SRP021221
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP021221") & (TissucombinedMetaes$Tissue=="Intestine"),"Tissue"]= "Colon"
#SRP009386
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "SRP009386") & (combinedMeta$Tissue=="Intestine"),"Tissue"]= "Colon"
#exclude SRP048804 (Cell line)
combinedMeta=combinedMeta[!combinedMeta$study=="SRP048804",]
#exclude the remaining sample of Intestine
combinedMeta=combinedMeta[!combinedMeta$Tissue=="Intestine",]

#Fix ERP109002 Annotations
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Heart", combinedMeta$sra.library_name, ignore.case=T)),"Tissue"]= "Heart"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Heart", combinedMeta$sra.library_name, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Kidney", combinedMeta$sra.library_name, ignore.case=T)),"Tissue"]= "Kidney"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Kidney", combinedMeta$sra.library_name, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Testis", combinedMeta$sra.library_name, ignore.case=T)),"Tissue"]= "Testis"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Testis", combinedMeta$sra.library_name, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Liver", combinedMeta$sra.library_name, ignore.case=T)),"Tissue"]= "Liver"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Liver", combinedMeta$sra.library_name, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Brain", combinedMeta$sra.library_name, ignore.case=T)),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Brain", combinedMeta$sra.library_name, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Forebrain", combinedMeta$sra.experiment_attributes, ignore.case=T)),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Forebrain", combinedMeta$sra.experiment_attributes, ignore.case=T)),"Tissue2"]= "Forebrain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Hindbrain", combinedMeta$sra.experiment_attributes, ignore.case=T)),"Tissue"]= "Brain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Hindbrain", combinedMeta$sra.experiment_attributes, ignore.case=T)),"Tissue2"]= "Hindbrain"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Ovary", combinedMeta$sra.library_name, ignore.case=T)),"Tissue"]= "Ovary"
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("Ovary", combinedMeta$sra.library_name, ignore.case=T)),"Tissue2"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002"),"Cellline"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002"),"CelllineName"]= ""
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002"),"Cancer"]= FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002"),"Fetal"] <- FALSE
combinedMeta[!is.na(combinedMeta$study) & (combinedMeta$study== "ERP109002") & (grepl("embryo", combinedMeta$sra.sample_attributes, ignore.case=T)),"Fetal"]= TRUE


studySamples <- combinedMeta$study %in% c("SRP105369", "ERP006121", "SRP062144")
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- "AML"
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- TRUE



studySamples <- combinedMeta$study %in% c("SRP221351", "SRP110313", "SRP115151", "SRP133278", "SRP156583", "SRP201603")
combinedMeta$Tissue[studySamples] <- "Blood"
combinedMeta$Tissue2[studySamples] <- "B-cells"
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- FALSE


studySamples <- combinedMeta$study %in% c("ERP107715", "ERP111116", "SRP092158", "SRP133442", "SRP065795", "SRP119636")
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- NA
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- TRUE

studySamples <- combinedMeta$study %in% c("ERP109703", "SRP100686", "SRP161505")
combinedMeta$Tissue[studySamples] <- "Blood"
combinedMeta$Tissue2[studySamples] <- "CLL"
combinedMeta$Cellline[studySamples] <- NA
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- TRUE



studySamples <- combinedMeta$study == "SRP123604"
combinedMeta$Tissue[studySamples] <- "Colon"
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- TRUE



studySamples <- combinedMeta$study %in% c("ERP113862", "ERP002323", "ERP114921", "SRP051368", "SRP097893", "SRP101856", "SRP151577")
combinedMeta$Tissue[studySamples] <- "Blood"
combinedMeta$Tissue2[studySamples] <- "Dendritic cells"
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- FALSE



studySamples <- combinedMeta$study %in% c("SRP056733", "SRP062278", "SRP064515", "SRP074274", "SRP076097", "SRP095287", "SRP103821", "SRP109107", "SRP110187", "SRP118741", "SRP118760", "SRP145599", "SRP190161", "SRP218274", "SRP155941")
combinedMeta$Tissue[studySamples] <- "Blood"
combinedMeta$Tissue2[studySamples] <- "Macrophages"
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- FALSE


studySamples <- combinedMeta$study %in% c("ERP020977", "ERP022909")
combinedMeta$Tissue[studySamples] <- "Blood"
combinedMeta$Tissue2[studySamples] <- "Macrophages-iPSC"
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- FALSE



studySamples <- combinedMeta$study %in% c("ERP014531", "SRP041826", "SRP058953", "SRP096201", "SRP113586", "SRP192825", "SRP173842", "SRP045352", "SRP055514", "SRP069333", "SRP101726")
combinedMeta$Tissue[studySamples] <- "Blood"
combinedMeta$Tissue2[studySamples] <- "Monocytes"
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- FALSE


studySamples <- combinedMeta$study %in% c("SRP150456")
combinedMeta$Tissue[studySamples] <- "Nasal Lavage"
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- FALSE



studySamples <- combinedMeta$study %in% c("SRP102104", "SRP162654", "SRP042596", "SRP049605", "SRP074736", "SRP090282", "SRP125882", "SRP140711", "SRP162023", "SRP168421", "SRP201023", "SRP212077", "SRP140558")
combinedMeta$Tissue[studySamples] <- "Blood"
combinedMeta$Tissue2[studySamples] <- "PBMC"
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- FALSE


studySamples <- combinedMeta$study %in% c("SRP072980", "SRP169062", "SRP086613", "SRP092010", "ERP105662", "SRP093990", "SRP215282", "SRP032926", "SRP053186", "SRP059057", "SRP098715", "SRP101784", "SRP117629", "SRP140710", "SRP155217", "SRP158900", "SRP192607")
combinedMeta$Tissue[studySamples] <- "Blood"
combinedMeta$Tissue2[studySamples] <- "T-cells"
combinedMeta$Cellline[studySamples] <- FALSE
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- FALSE

#SRP081020 is mis-anotated in SRA as PBMC, paper and clustering both state wholeblood
studySamples <- combinedMeta$study %in% c("ERP114104", "SRP051848", "SRP056784", "SRP071965", "SRP077975", "SRP081020", "SRP098758", "SRP113245", "SRP126580", "SRP126582", "SRP126583", "SRP136057", "SRP144583", "SRP150872", "SRP214077", "SRP056443")
combinedMeta$Tissue[studySamples] <- ""
combinedMeta$Tissue2[studySamples] <- ""
combinedMeta$Cellline[studySamples] <- NA
combinedMeta$CelllineName[studySamples] <- ""
combinedMeta$Cancer[studySamples] <- NA






combinedMeta$Cellline[!is.na(combinedMeta$CelllineName)&combinedMeta$CelllineName=="iPSC"] <- TRUE
combinedMeta$Cancer[!is.na(combinedMeta$Cellline) & combinedMeta$Cellline] <- NA


table(paste0(combinedMeta$Tissue, " - ", combinedMeta$Tissue2, " - ", combinedMeta$Fetal)[combinedMeta$study == "ERP109002"])

#Exclude spike in
combinedMeta$exclude[combinedMeta$study == "SRP041955"] <- TRUE

(x <- table(paste0(combinedMeta$Tissue, " - ", combinedMeta$Tissue2),combinedMeta$Cancer))
write.table(x, file = "test.txt", row.names = T, col.names = NA, quote = F, sep = "\t")


table(paste0(combinedMeta$Tissue, " - ", combinedMeta$Tissue2),combinedMeta$Cancer)

table(combinedMeta$Tissue)
table(combinedMeta$Tissue2)

table(combinedMeta$Tissue, combinedMeta$Cellline)

table(combinedMeta$Cancer)

table(combinedMeta$Tissue, combinedMeta$Cancer)

table(combinedMeta$Tissue[combinedMeta$Cohort == "TCGA"], combinedMeta$Cancer[combinedMeta$Cohort == "TCGA"])


sum(combinedMeta$gado.TissueType %in% combinedMeta$gado.PlotClass)

sum((!is.na(combinedMeta$gado.TissueType) & combinedMeta$gado.TissueType != ""))

sum( )
unique(combinedMeta$gado.TissueType)

table(combinedMeta$gado.PlotClass, useNA = "a")

combinedMeta$Tissue[combinedMeta$Cohort == "GSA"]



#save(combinedMeta, file = "combinedMeta_2022_09_15.RData")

load(file = "combinedMeta_2022_08_19.RData")

pcsAndMeta <- merge(expPcs[,1:100], combinedMeta, by = 0, all.x = T)
dim(pcsAndMeta)
str(combinedMeta)

tissueCol <- read.delim("Recount3_QC_2ndRun/SRA_Studies_Annotations_Patrick/Annotations_color2.txt", row.names = 1)

sum(unique(pcsAndMeta[,"Tissue"]) %in% tissueCol$PlotClass)
sum(unique(pcsAndMeta[,"Tissue2"]) %in% tissueCol$PlotClass)

x <- unique(pcsAndMeta[,"Tissue2"])
x[!x %in% tissueCol$PlotClass]


defaultCol <- adjustcolor("grey", alpha.f = 0.6)
pcsAndMeta$col <- defaultCol

tissueAndCol <- pcsAndMeta[,"Tissue"] %in% tissueCol$PlotClass

pcsAndMeta$col[tissueAndCol] <-  adjustcolor(tissueCol$col[match(pcsAndMeta[tissueAndCol,"Tissue"], tissueCol$PlotClass)], alpha.f = 0.6)


tissue2AndCol <- pcsAndMeta[,"Tissue2"] %in% tissueCol$PlotClass
sum(tissue2AndCol)
pcsAndMeta$col[tissue2AndCol] <-  adjustcolor(tissueCol$col[match(pcsAndMeta[tissue2AndCol,"Tissue2"], tissueCol$PlotClass)], alpha.f = 0.6)

table(pcsAndMeta[pcsAndMeta[,"PC_2"] >= 0,"Tissue2"])


sum(is.na(tolower(pcsAndMeta[,"Tissue"]) %in% tolower(tisueCol$PlotClass)))

#pcsAndMeta$col <- tissueCol$col[match(tolower(pcsAndMeta[,"Tissue"]), tolower(tissueCol$PlotClass), nomatch = nrow(tissueCol))]

plotOrder <- order((pcsAndMeta$col != defaultCol) + 1)

rpng(width = 800, height = 800)
#pdf(file = "test.pdf")
plot(pcsAndMeta[plotOrder,"PC_1"], pcsAndMeta[plotOrder,"PC_2"], col = pcsAndMeta$col[plotOrder], cex = 0.3, pch = 16)
dev.off()


rpng(width = 800, height = 800)
#pdf(file = "test.pdf")
plot(pcsAndMeta[plotOrder,"PC_3"], pcsAndMeta[plotOrder,"PC_7"], col = pcsAndMeta$col[plotOrder], cex = 0.3, pch = 16)
dev.off()


#rpng(width = 800, height = 800)
png("tissues.png",width = 2000, height = 2000)
pairs(pcsAndMeta[plotOrder,paste0("PC_",1:5)], col = pcsAndMeta$col[plotOrder], cex = 0.4, upper.panel = NULL, pch = 16)
dev.off()

png("tissues2.png",width = 2000, height = 2000)
pairs(pcsAndMeta[plotOrder,paste0("PC_",6:10)], col = pcsAndMeta$col[plotOrder], cex = 0.4, upper.panel = NULL)
dev.off()

png("tissues3.png",width = 2000, height = 2000)
pairs(pcsAndMeta[plotOrder,paste0("PC_",11:15)], col = pcsAndMeta$col[plotOrder], cex = 0.4, upper.panel = NULL)
dev.off()

defaultCol <- adjustcolor("grey", alpha.f = 0.3)
pcsAndMeta$colCelline <- defaultCol
pcsAndMeta$colCelline[!is.na(pcsAndMeta[,"Cellline"]) & pcsAndMeta[,"Cellline"]] <-  adjustcolor("magenta", alpha.f = 0.3)
pcsAndMeta$colCelline[!is.na(pcsAndMeta[,"Cellline"]) & !pcsAndMeta[,"Cellline"]] <-  adjustcolor("royalblue1", alpha.f = 0.3)
pcsAndMeta$colCelline[!is.na(pcsAndMeta[,"Cancer"]) & pcsAndMeta[,"Cancer"]] <-  adjustcolor("forestgreen", alpha.f = 0.3)
plotOrder <- order((pcsAndMeta$colCelline != defaultCol) + 1)


pcsAndMeta$cellineTissueCancer <- "Unkown"
pcsAndMeta$cellineTissueCancer[!is.na(pcsAndMeta[,"Cellline"]) & pcsAndMeta[,"Cellline"]] <- "Cellline"
pcsAndMeta$cellineTissueCancer[!is.na(pcsAndMeta[,"Cellline"]) & !pcsAndMeta[,"Cellline"]] <- "Tissue"
pcsAndMeta$cellineTissueCancer[!is.na(pcsAndMeta[,"Cancer"]) & pcsAndMeta[,"Cancer"]] <- "Cancer"

pcsAndMeta$cellineTissueCancer <- factor(pcsAndMeta$cellineTissueCancer, levels = c("Tissue", "Cancer", "Cellline", "Unkown"))

table(pcsAndMeta$cellineTissueCancer, useNA = "always")
  
rpng(width = 800, height = 800)
plot(pcsAndMeta[plotOrder,"PC_1"], pcsAndMeta[plotOrder,"PC_2"], col = pcsAndMeta$colCelline[plotOrder], cex = 0.4, pch = 16)
dev.off()

rpng(width = 800, height = 800)
plot(pcsAndMeta[plotOrder,"PC_3"], pcsAndMeta[plotOrder,"PC_75"], col = pcsAndMeta$colCelline[plotOrder], cex = 0.4, pch = 16)
dev.off()

for(i in c(1,3:100)){
  png(paste0("cellinePlots/pc",i,".png"),width = 1000, height = 1000)
  #rpng()
  plot(pcsAndMeta[plotOrder,"PC_2"], pcsAndMeta[plotOrder,paste0("PC_",i)], col = pcsAndMeta$colCelline[plotOrder], cex = 1, pch = 16, xlab = "PC2", ylab = paste0("PC", i))
  dev.off()
}

library(vioplot)

for(i in 1:100){
png(paste0("cellinePlots2/pc",i,".png"),width = 500, height = 500)
vioplot( pcsAndMeta[,paste0("PC_",i)] ~ pcsAndMeta$cellineTissueCancer, col = c(adjustcolor("royalblue1", alpha.f = 0.3), adjustcolor("forestgreen", alpha.f = 0.3), adjustcolor("magenta", alpha.f = 0.3), defaultCol))
dev.off()
}
table(paste0(combinedMeta$Tissue, " - ", combinedMeta$Tissue2))

png("celllines_c.png",width = 2000, height = 2000)
pairs(pcsAndMeta[plotOrder,paste0("PC_",1:5, "_c")], col = pcsAndMeta$colCelline[plotOrder], cex = 0.4, upper.panel = NULL, pch = 16)
dev.off()


png("celllines2.png",width = 2000, height = 2000)
pairs(pcsAndMeta[plotOrder,paste0("PC_",6:10, "")], col = pcsAndMeta$colCelline[plotOrder], cex = 0.4, upper.panel = NULL, pch = 16)
dev.off()


defaultCol <- adjustcolor("grey", alpha.f = 0.3)
pcsAndMeta$colCancer <- defaultCol
pcsAndMeta$colCancer[!is.na(pcsAndMeta[,"Cancer"]) & pcsAndMeta[,"Cancer"]] <- adjustcolor("chartreuse1", alpha.f = 0.6)
plotOrder <- order((pcsAndMeta$colCancer != defaultCol) + 1)

rpng(width = 800, height = 800)
plot(pcsAndMeta[plotOrder,"PC_1"], pcsAndMeta[plotOrder,"PC_2"], col = pcsAndMeta$colCancer[plotOrder], cex = 0.4)
dev.off()

png("cancers.png",width = 2000, height = 2000)
pairs(pcsAndMeta[plotOrder,paste0("PC_",1:5)], col = pcsAndMeta$colCancer[plotOrder], cex = 0.4, upper.panel = NULL, pch = 16)
dev.off()

library(pROC)
cancerAuc <- apply(pcsAndMeta[,paste0("PC_",1:100)], 2, function(x){as.numeric(auc(response = pcsAndMeta$Cancer, predictor = x))})
sort(cancerAuc)

rpng(width = 800, height = 800)
plot(pcsAndMeta[plotOrder,"PC_33"], pcsAndMeta[plotOrder,"PC_9"], col = pcsAndMeta$col[plotOrder], cex = 0.4)
dev.off()


library(pROC)
cancerAuc <- apply(pcsAndMeta[,paste0("PC_",1:100)], 2, function(x){as.numeric(auc(response = pcsAndMeta$Cancer, predictor = x))})
sort(cancerAuc)

celllineAuc <- apply(pcsAndMeta[,paste0("PC_",1:100)], 2, function(x){as.numeric(auc(response = pcsAndMeta$Cellline, predictor = x))})
sort(celllineAuc)


rpng(width = 800, height = 800)
plot(pcsAndMeta[plotOrder,"PC_33"], pcsAndMeta[plotOrder,"PC_9"], col = pcsAndMeta$colCancer[plotOrder], cex = 0.4)
dev.off()

rpng()
pairs(pcsAndMeta[plotOrder,paste0("PC_",c(33,32,9,10,21))], col = pcsAndMeta$colCancer[plotOrder], cex = 0.4, upper.panel = NULL)
dev.off()


rpng()
pairs(pcsAndMeta[plotOrder,paste0("PC_",c(3,2,27,6,20))], col = pcsAndMeta$colCelline[plotOrder], cex = 0.4, upper.panel = NULL)
dev.off()

combinedMeta$sra.sample_title <- gsub("\"", "", combinedMeta$sra.sample_title)

tmp <- merge(combinedMeta,pcs[,1:100], by = 0, all.y = T)
dim(tmp)
write.table(tmp, file = "tmpAnnotations.txt", sep = "\t", quote = FALSE, col.names = NA)

qseq <- read.delim("quantseqSamples.txt")[,1]
str(qseq)


defaultCol <- adjustcolor("grey", alpha.f = 0.3)
pcsAndMeta$colQseq <- defaultCol
pcsAndMeta$colQseq[pcsAndMeta$Row.names %in% qseq] <- "orangered"
plotOrderQseq <- order((pcsAndMeta$colQseq != defaultCol) + 1)

plot(pcsAndMeta[plotOrderQseq,"PC_1"], pcsAndMeta[plotOrderQseq,"PC_2"], col = pcsAndMeta$colQseq[plotOrderQseq], cex = 0.4, main = "Quantseq")

plot(pcsAndMeta[plotOrderQseq,"PC_6"], pcsAndMeta[plotOrderQseq,"PC_1"], col = pcsAndMeta$colQseq[plotOrderQseq], cex = 0.4, main = "Quantseq")


table(pcsAndMeta$sra.library_layout)


numColumns <- unlist(lapply(combinedMeta, is.numeric))


combinedMetaMatrix <- as.matrix(combinedMeta[,numColumns])

library(pROC)

qseqClass <- as.factor(rownames(combinedMeta) %in% qseq)
table(qseqClass)
dim(combinedMetaMatrix)
qseqPValues <- apply(combinedMetaMatrix,2,function(x){
  tryCatch(
    {
      #wilcox.test(x ~ qseqClass)$p.value
      as.numeric(auc(response = qseqClass, predictor = x))
    },
    error=function(cond){return(1)}
  )
})
sort(qseqPValues)


boxplot(combinedMetaMatrix[,"recount_qc.bc_frag.kallisto_mean_length"] ~ qseqClass )








plot(pcsAndMeta[plotOrderQseq,"recount_seq_qc.%c"], log10(pcsAndMeta[plotOrderQseq,"recount_qc.star.number_of_splices:_total"]), col = pcsAndMeta$colQseq[plotOrderQseq], cex = 0.6)
10^5.2
abline(h=log10(150000))
abline(v=60)
log10(10^5)

plot(log10(pcsAndMeta[plotOrderQseq,"recount_qc.star.number_of_splices:_total"]), pcsAndMeta[plotOrderQseq,"recount_qc.bc_frag.kallisto_mean_length"], col = pcsAndMeta$colQseq[plotOrderQseq], cex = 0.4)

plot(pcsAndMeta[plotOrderQseq,"PC_6"], log10(pcsAndMeta[plotOrderQseq,"recount_qc.star.number_of_splices:_total"]), col = pcsAndMeta$colQseq[plotOrderQseq], cex = 0.4)


pc1Cor <- cor(pcsAndMeta[,"PC_1"], pcsAndMeta[,numColumns], use = "pairwise.complete.obs")
sort(pc1Cor[1,])

pc6Cor <- apply(combinedMetaMatrix[pcsAndMeta$Row.names,],2,function(x){
  tryCatch(
    {
      #wilcox.test(x ~ qseqClass)$p.value
      cor(pcsAndMeta[,"PC_6"], x, use = "pairwise.complete.obs")
    },
    error=function(cond){return(0)},
    warning=function(cond){return(0)}
  )
})
sort(pc6Cor^2)

load("testPcPAtrickFrist100.RData", verbose = T)
colnames(expPcs)
str(expPcs)
colnames(expPcs) <- paste0("PC_", 1:ncol(expPcs))
pcsAndMeta <- merge(expPcs, combinedMeta, by = 0, all.x = T)
dim(pcsAndMeta)



load("gadoPca.RData", verbose = T)
colnames(expGadoPcsSub)
str(expGadoPcsSub)
colnames(expGadoPcsSub) <- paste0("PC_", 1:ncol(expGadoPcsSub))
pcsAndMeta <- merge(expGadoPcsSub, combinedMeta, by = 0, all.x = T)
dim(pcsAndMeta)



table(pcsAndMeta$CelllineName)
pcsAndMeta$Cellline[grepl("s4u", pcsAndMeta$CelllineName)]






