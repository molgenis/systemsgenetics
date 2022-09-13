#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)


remoter::client("localhost", port = 55501, password = "laberkak")



library(uwot)

setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Recount3\\")
setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

tissueCol <- read.delim("umap/col.txt", row.names = 1, na.strings = "")

load(file = "DataForPredictions.RData")

#load(file = "combinedMeta_2022_08_30.RData", verbose = T)
#str(combinedMeta)
#updatedAnnotations <- combinedMeta[,c("Tissue",	"Tissue2",	"Cellline",	"CelllineName",	"Cancer",	"Cohort", "Fetal")]

#all(rownames(pcsAndMeta) %in% rownames(updatedAnnotations))
#updatedAnnotations <- updatedAnnotations[rownames(pcsAndMeta),]
#all(rownames(pcsAndMeta) == rownames(updatedAnnotations))

#pcsAndMeta[,colnames(updatedAnnotations)] <- updatedAnnotations

#pcsAndMeta$selectedSamples <- !pcsAndMeta$excludeBasedOnPredictionCellline2 & !pcsAndMeta$excludeBasedOnPredictionCancer & !(!is.na(pcsAndMeta$Cancer) & pcsAndMeta$Cancer) & !(!is.na(pcsAndMeta$Cellline) & pcsAndMeta$Cellline)

table(pcsAndMeta$selectedSamples, useNA = "a")


clusterAnnotations <- read.delim("umap/annotationsBasedOnOldUmap.txt", row.names = 1)
pcsAndMeta <- merge(pcsAndMeta, clusterAnnotations, by = 0, all.x = T)
rownames(pcsAndMeta) <- pcsAndMeta$Row.names
table(pcsAndMeta$ClusterAnnotation)




#pcsAndMeta[!is.na(pcsAndMeta$study) & (pcsAndMeta$study== "ERP104864") & (grepl("synovium", pcsAndMeta$sra.sample_attributes, ignore.case=T)),"Tissue2"]= ""


tissueSamples <- pcsAndMeta[pcsAndMeta$selectedSamples,]

tissueSamples$class <- tissueSamples$Tissue

hasT2 <- tissueSamples$Tissue2 != ""
tissueSamples$class[hasT2] <- paste0(tissueSamples$class[hasT2], "-", tissueSamples$Tissue2[hasT2])

isFetal <- !is.na(tissueSamples$Fetal) & tissueSamples$Fetal
tissueSamples$class[isFetal] <- paste0(tissueSamples$class[isFetal], "-Fetal")

noTbutCluster <- tissueSamples$class == "" & !is.na(tissueSamples$ClusterAnnotation)
table(noTbutCluster, useNA = "a")
tissueSamples$class[noTbutCluster] <- tissueSamples$ClusterAnnotation[noTbutCluster]

table(tissueSamples$class)
write.table(table(tissueSamples$class, useNA = "always"), file = "umap/tissues.txt", sep = "\t", quote = F, row.names = F)

str(tissueSamples)



mapping <- read.delim("umap/tissuesMapping.txt")
str(mapping)

all(tissueSamples$class %in% mapping$Class)


tissueSamples$umapFactor <- as.factor(mapping$ClassificationClass[match(tissueSamples$class, mapping$Class)])

table(tissueSamples$umapFactor, useNA = "always")


defaultCol <- adjustcolor("grey", alpha.f = 0.6)
tissueCol <- read.delim("umap/col.txt", row.names = 1)


tissueSamples$TissueCol <- defaultCol
sum(unique(tissueSamples$umapFactor) %in% rownames(tissueCol))
sum(tissueSamples$umapFactor %in% rownames(tissueCol))
tissueSamples$TissueCol[tissueSamples$umapFactor %in% rownames(tissueCol)] <- adjustcolor(tissueCol[as.character(tissueSamples$umapFactor[tissueSamples$umapFactor %in% rownames(tissueCol)]),1], alpha.f = 0.5)
#tissueSamples$TissueCol[tissueSamples$umapFactor %in% rownames(tissueCol)] <- tissueCol[as.character(tissueSamples$umapFactor[tissueSamples$umapFactor %in% rownames(tissueCol)]),1]
table(tissueSamples$TissueCol, useNA = "a")

tissueSamples$plotOrderTissues <- order(tissueSamples$TissueCol != defaultCol)


#, n_threads = 22

compsToUseForUmap <- compsToUse
init <- as.matrix(tissueSamples[,paste0("PC_",1:2)])
umapInput <- as.matrix(tissueSamples[,paste0("PC_",1:compsToUseForUmap)])

sampleUmap <- umap(
  umapInput, 
  n_epochs = 1000, 
  init = init, 
  n_neighbors = 500, 
  min_dist = 2, init_sdev = 1e-4, learning_rate = 1, 
  spread = 15, 
  bandwidth = 10,
  scale = "scale",
  local_connectivity = 1,
  metric = "correlation")


rownames(sampleUmap) <- rownames(tissueSamples)
colnames(sampleUmap) <- c("UMAP1", "UMAP2")
umapAndMeta <- merge(sampleUmap, tissueSamples, by = 0)
dim(umapAndMeta)





rpng()

par(mar = c(3,5,0.1,0.1), xpd = NA)
plot(umapAndMeta[umapAndMeta$plotOrderTissues,"UMAP1"], umapAndMeta[umapAndMeta$plotOrderTissues,"UMAP2"], col = umapAndMeta$TissueCol[umapAndMeta$plotOrderTissues], cex = 0.2, pch = 16)

dev.off()


locator(n =2, type = "l")
cluster1 <- locator(n =2, type = "l")
cluster2 <- locator(n =2, type = "l")


write.table(umapAndMeta,file = "umaptest.txt", sep = "\t", quote = F, col.names = NA)

#save.image( file="umap_tmp.RData")
#load("umap_tmp.RData") 

rpng()

par(mar = c(3,5,0.1,0.1), xpd = NA)
plot(umapAndMeta[plotOrder,"UMAP1"], umapAndMeta[plotOrder,"UMAP2"], col = umapAndMeta$TissueCol[plotOrder], cex = 0.8, pch = 16, xlim = c(-25,25), ylim = c(-25,25))

dev.off()



#png(file = "umaptest.png", width = 1600, height = 800)

pdf(file = "umaptest.pdf", width = 16, height = 8)
#rpng()

layout(matrix(1:2,ncol = 2))

par(mar = c(3,5,0.1,0.1), xpd = NA)
plot(umapAndMeta[plotOrderTissues,"UMAP1"], umapAndMeta[plotOrderTissues,"UMAP2"], col = umapAndMeta$TissueCol[plotOrderTissues], cex = 0.4, pch = 16)

par(mar = c(0,0,0,0), xpd = NA)
plot.new()
plot.window(xlim = 0:1, ylim = 0:1)
legend("center", fill = tissueCol[,1], legend = row.names(tissueCol), bty = "n", ncol = 2,cex = 0.7)


dev.off()







#smartseq plots

someSmartSeqStudies <- read.delim("selectionSmartseqStudies.txt", header = F)[,1]
str(someSmartSeqStudies)

someSmartSeqSamples <- read.delim("smartseqSamples.txt", header = T)[,1]
str(someSmartSeqSamples)

umapAndMeta$smartseqcol <- defaultCol
umapAndMeta$smartseqcol[umapAndMeta$study %in% someSmartSeqStudies] <- "pink"
umapAndMeta$smartseqcol[umapAndMeta$Row.names %in% someSmartSeqSamples] <- "pink"

umapAndMeta$plotOrdersq <- order(umapAndMeta$smartseqcol != defaultCol)


par(mar = c(3,5,0.1,0.1), xpd = NA)
plot(umapAndMeta[umapAndMeta$plotOrderTissues,"UMAP1"], umapAndMeta[umapAndMeta$plotOrderTissues,"UMAP2"], col = umapAndMeta$smartseqcol[umapAndMeta$plotOrderTissues], cex = 0.2, pch = 16)


#instestine clusters

umapAndMeta$intestineCluster <- ""
umapAndMeta$intestineCluster[umapAndMeta$UMAP1 >= cluster1$x[1] & umapAndMeta$UMAP1 <= cluster1$x[2] & umapAndMeta$UMAP2 >= cluster1$y[1] & umapAndMeta$UMAP2 <= cluster1$y[2]] <- "c1"
umapAndMeta$intestineCluster[umapAndMeta$UMAP1 >= cluster2$x[1] & umapAndMeta$UMAP1 <= cluster2$x[2] & umapAndMeta$UMAP2 >= cluster2$y[1] & umapAndMeta$UMAP2 <= cluster2$y[2]] <- "c2"
table(umapAndMeta$intestineCluster)

table(factor(umapAndMeta$umapFactor[umapAndMeta$intestineCluster=="c1"]))
table(factor(umapAndMeta$umapFactor[umapAndMeta$intestineCluster=="c2"]))

table(factor(umapAndMeta$class[umapAndMeta$intestineCluster=="c1"]))
table(factor(umapAndMeta$class[umapAndMeta$intestineCluster=="c2"]))

a <- as.data.frame(table(paste(umapAndMeta$Cohort, umapAndMeta$class)[umapAndMeta$intestineCluster=="c1"]))
b <- as.data.frame(table(paste(umapAndMeta$Cohort, umapAndMeta$class)[umapAndMeta$intestineCluster=="c2"]))


table(paste(umapAndMeta$Cohort, umapAndMeta$class)[umapAndMeta$intestineCluster!=""], umapAndMeta$intestineCluster[umapAndMeta$intestineCluster!=""])

str(a)
c <- merge(a,b,by = 0, all = T)
c

load("metadata_gtex.Rda", verbose = T)
View(metadata_gtex)


gtexTansverse <- umapAndMeta[umapAndMeta$study == "GTEx" & umapAndMeta$Tissue2 == "Transverse" & umapAndMeta$intestineCluster != "",]

rownames(gtexTansverse) <- gtexTansverse$Row.names

rownames(metadata_gtex) <- metadata_gtex$external_id

dim(gtexTansverse)
gtexTansverse <- merge(gtexTansverse, metadata_gtex[,!colnames(metadata_gtex) %in% colnames(gtexTansverse)], by = 0)
dim(gtexTansverse)

table(gtexTansverse$gtex.smatsscr, gtexTansverse$intestineCluster)

fisher.test(table(gtexTansverse$gtex.smatsscr, gtexTansverse$intestineCluster))
grep("MHBCTINF", colnames(gtexTansverse), ignore.case = T)




numCols <- colnames(gtexTansverse)[unlist(lapply(gtexTansverse, is.numeric))  ]

colName <-  "sra.paired_nominal_length"
clusterCompare <- sapply(numCols, function(colName){
  #print(colName)
  if(!all(is.na(gtexTansverse[,colName])) & sd(gtexTansverse[,colName], na.rm =T) > 0  ){
    t.test(gtexTansverse[,colName] ~ gtexTansverse$intestineCluster)$p.value
  }
  
})
clusterCompare <- unlist(clusterCompare)
clusterCompare2 <- clusterCompare[grep("PC_", names(clusterCompare), invert = T)]
sort(clusterCompare2, decreasing = T)
boxplot(gtexTansverse$`recount_qc.aligned_reads%.chrx` ~ gtexTansverse$intestineCluster)
boxplot(gtexTansverse$`recount_qc.aligned_reads%.chrx` ~ paste0(gtexTansverse$intestineCluster, "_",gtexTansverse$gtex.sex))
boxplot(gtexTansverse$`recount_qc.aligned_reads%.chrm` ~ gtexTansverse$intestineCluster)
boxplot(gtexTansverse$`` ~ gtexTansverse$intestineCluster)

boxplot(gtexTansverse$`recount_qc.star.number_of_reads_unmapped:_other_both` ~ gtexTansverse$intestineCluster)

boxplot(gtexTansverse$`gtex.smtsisch` ~ gtexTansverse$intestineCluster)
boxplot(gtexTansverse$`CnvAutoCor` ~ gtexTansverse$intestineCluster)

#save(gtexTansverse, file =  "gtexTansverse.RData")
load("gtexTansverse.RData")


str(row.names(gtexTansverse))
str(gtexTansverse$Row.names)
str(exp)
expgT <- exp[,gtexTansverse$Row.names]
save(expgT, file = "expgT.RData")
load( "expgT.RData")


colnames(expgT)
expgT <- t(expgT)
all(rownames(expgT) == gtexTansverse$Row.names)

x <- expgT[,1]

diffExp <- apply(expgT, 2, function(x){
  t.test(x ~gtexTansverse$intestineCluster)$statistic
})
hist(-log10(diffExp))
names(diffExp)[order(diffExp)[1:100]]
cat(sub("\\..+","",names(diffExp)[order(diffExp, decreasing = T)[1:200]]), sep = "\n")

load("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/SRA_Studies_Annotations_Patrick/Fibroblasts.rda", verbose = T)
str(fibroblasts)

load("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/Recount3_QC_2ndRun/SRA_Studies_Annotations_Patrick/BloodVessels.rda", verbose = T)




minSamplesTraining <- 50
maxFractionOfStudy <- 0.8


umapAndMetaClassified <- umapAndMeta[!is.na(umapAndMeta$umapFactor),]
umapAndMetaClassified$training <- FALSE

tissueClass <- levels(umapAndMetaClassified$umapFactor)[2]
study <- "GTEx"

set.seed(42)
#for each tissue slecect samples for training
for(tissueClass in levels(umapAndMetaClassified$umapFactor)){
  thisTissueSamples <- umapAndMetaClassified$umapFactor==tissueClass
  studiesForThisTissue <- unique(umapAndMetaClassified$study[thisTissueSamples])
  numberOfStudies <- length(studiesForThisTissue)
  numberOfSamplesPerStudy <- ceiling(minSamplesTraining / numberOfStudies)
  print(paste(tissueClass, length(studiesForThisTissue), numberOfSamplesPerStudy, sep = " - "))
  #for each studies put samples to training or test
  for(study in studiesForThisTissue){
    
    thisTissueAndStudySamples <- thisTissueSamples & umapAndMetaClassified$study == study
    thisTissueAndStudySamplesCount <- sum(thisTissueAndStudySamples)
    
    #Don't select more samples from study then the study has and also no more then set fraction. Do floor to put studies with single sample to testset
    potentialMax <- floor(thisTissueAndStudySamplesCount * maxFractionOfStudy)
    numberTrainingSamplesThisStudy <- if(potentialMax > numberOfSamplesPerStudy) numberOfSamplesPerStudy  else potentialMax
    if(numberTrainingSamplesThisStudy > 0){
        #The which will get all indices for the samples of this study-tissue combination. These are then samples for the samples used for training
        trainingSamplesThisStudy <- sample(which(thisTissueAndStudySamples), numberTrainingSamplesThisStudy)
        #Set selected to TRUE
        umapAndMetaClassified$training[trainingSamplesThisStudy] <- TRUE
    }
    
    
   #print(paste0(thisTissueAndStudySamplesCount, " - ", numberTrainingSamplesThisStudy))
  }
  
}

sum(umapAndMetaClassified$training)

umapAndMetaClassifiedTraining <- umapAndMetaClassified[umapAndMetaClassified$training,]
table(umapAndMetaClassifiedTraining$umapFactor)
umapAndMetaClassifiedTest <- umapAndMetaClassified[!umapAndMetaClassified$training,]
dim(umapAndMetaClassifiedTest)



cfit <- cv.glmnet(x = as.matrix(umapAndMetaClassifiedTraining[,paste0("PC_",1:compsToUse)]), y = umapAndMetaClassifiedTraining$umapFactor, family = "multinomial", type.measure = "class", alpha=1, nlambda=100)
best_lambda <- cfit$lambda.min
cfit



fibTraining <- fibroblasts$Row.names
bvTraining <-  bloodVessels$Row.names   
bvTraining
unique(bvTraining)
