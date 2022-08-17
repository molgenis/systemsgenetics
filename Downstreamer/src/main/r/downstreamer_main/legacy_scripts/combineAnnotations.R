#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)


remoter::client("localhost", port = 55501, password = "laberkak")


#save.image("tmp.RData")



setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")
setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Recount3\\")
load("tmp.RData")


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


#now fillin the gtex and gcta recount meta data. 

tmp <- metadata_gtex[,colnames(metadata_gtex) %in% sraSharedCol]
combinedMeta[rownames(tmp),colnames(tmp)] <- tmp

tmp <- metadata_tcga[,colnames(metadata_tcga) %in% sraSharedCol]
combinedMeta[rownames(tmp),colnames(tmp)] <- tmp

rm(tmp)

#set study make column uniform
combinedMeta$study[combinedMeta$Cohort == "GTEx"] <- "GTEx"
combinedMeta$study[combinedMeta$Cohort == "TCGA"] <- "TCGA"


head(rownames(combinedMeta))

combinedMeta$Fetal <- NA

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


combinedMeta["SRR5341594", "sra.sample_title"] <- "Human differentiating macrophage"


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


combinedMeta$Cellline[!is.na(combinedMeta$CelllineName)&combinedMeta$CelllineName=="iPSC"] <- TRUE







(x <- table(paste0(combinedMeta$Tissue, " - ", combinedMeta$Tissue2),combinedMeta$Cancer))
str(x)
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




#save(combinedMeta, file = "combinedMeta_2022_08_16.RData")

load(file = "combinedMeta_2022_08_16.RData")

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






