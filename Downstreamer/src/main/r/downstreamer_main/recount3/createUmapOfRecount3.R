#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55001, password = "laberkak", sync = T)


#remoter::client("localhost", port = 55001, password = "laberkak")



library(uwot)

setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Recount3\\")
setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")


#load(file = "DataForPredictions.RData")
#rownames(pcsAndMeta) <- pcsAndMeta$Row.names


load("tissuePredictions/samplesWithPrediction_16_09_22_noOutliers.RData", verbose = T)
tissueCol <- read.delim("umap/col.txt", row.names = 1, na.strings = "")
rownames(tissueCol)

tissueCol2 <- tissueCol$Col
names(tissueCol2) <- rownames(tissueCol)

all(samplesWithPredictionNoOutliers$predictedTissue %in% names(tissueCol2))


samplesWithPredictionNoOutliers$col <- tissueCol2[samplesWithPredictionNoOutliers$predictedTissue]


load(paste0("CombinedHealthyTissue/combinedHealthyTissue_PCA.RData"))


all(rownames(combinedHealtyTissuePca$expPcs) == rownames(samplesWithPredictionNoOutliers))


str(samplesWithPredictionNoOutliers)
dim(combinedHealtyTissuePca$expPcs)

combinedHealtyTissuePca$expPcs <- combinedHealtyTissuePca$expPcs[,1:50]
str(combinedHealtyTissuePca$expPcs )

nnData <- umap(combinedHealtyTissuePca$expPcs, ret_nn = T)
nnDataCorrelation <- umap(combinedHealtyTissuePca$expPcs, ret_nn = T, metric = "correlation")
nnDataCosine <- umap(combinedHealtyTissuePca$expPcs, ret_nn = T, metric = "cosine")

str(nnData)

#nn <- nnDataCorrelation$nn[[1]]
#nn <- nnData$nn[[1]]
#nn <- nnDataCosine$nn[[1]]
init <- combinedHealtyTissuePca$expPcs[,1:2]

sampleUmap <- umap(X = NULL, nn_method = nn)
plot(sampleUmap[,1],sampleUmap[,2], pch = 16, col=adjustcolor(samplesWithPredictionNoOutliers$col, alpha.f = 0.5), bty="n", xlab = "UMAP-1", ylab = "UMAP-2", cex = 0.7)



sampleUmap <- umap(X = NULL, nn_method = nn, bandwidth = 10, n_epochs = 500)
plot(sampleUmap[,1],sampleUmap[,2], pch = 16, col=adjustcolor(samplesWithPredictionNoOutliers$col, alpha.f = 0.5), bty="n", xlab = "UMAP-1", ylab = "UMAP-2", cex = 0.7)


sampleUmap <- umap(X = NULL, nn_method = nnDataCorrelation$nn[[1]], init = init, 
                    bandwidth = 500,n_epochs = 200,learning_rate = 10,
                    n_neighbors = 50,local_connectivity = 1)#, n_neighbors = 50,n_epochs = 1000,learning_rate = 10,local_connectivity = 50,bandwidth = 1000, repulsion_strength = 0.1
rownames(sampleUmap) <- rownames(combinedHealtyTissuePca)
colnames(sampleUmap) <- c("UMAP1", "UMAP2")

plot(sampleUmap[,1],sampleUmap[,2], pch = 16, col=adjustcolor(samplesWithPredictionNoOutliers$col, alpha.f = 0.3), bty="n", xlab = "UMAP-1", ylab = "UMAP-2", cex = 0.7)




sampleUmap <- umap(X = NULL, nn_method = nnDataCorrelation$nn[[1]], init = init, 
                   bandwidth = 100,n_epochs = 200,learning_rate = 10,
                   n_neighbors = 10,local_connectivity = 10, spread = 10)#, n_neighbors = 50,n_epochs = 1000,learning_rate = 10,local_connectivity = 50,bandwidth = 1000, repulsion_strength = 0.1
rownames(sampleUmap) <- rownames(combinedHealtyTissuePca)
colnames(sampleUmap) <- c("UMAP1", "UMAP2")


par(pty="s")
plot(sampleUmap[,1],sampleUmap[,2], pch = 16, col=adjustcolor(samplesWithPredictionNoOutliers$col, alpha.f = 0.3), bty="n", xlab = "UMAP-1", ylab = "UMAP-2", cex = 0.7)



tissueCenters <- aggregate(sampleUmap[,1:2], by = list(samplesWithPredictionNoOutliers$annotatedTissue), FUN = median)
tissueCenters <- merge(tissueCenters, tissueCol, by.x = "Group.1", by.y = 0)
str(tissueCenters)



#points(tissueCenters$UMAP1, tissueCenters$UMAP2, col = tissueCenters$Col, cex = 6)

pdf("umap/umapWithLables.pdf", width = 10, height = 10, useDingbats = F)
par(pty="s", xpd = NA)
plot(sampleUmap[,1],sampleUmap[,2], pch = 16, col=adjustcolor(samplesWithPredictionNoOutliers$col, alpha.f = 0.3), bty="n", xlab = "UMAP-1", ylab = "UMAP-2", cex = 0.7)
text(tissueCenters$UMAP1, tissueCenters$UMAP2, tissueCenters$Group.1, col = tissueCenters$Col, adj = c(0.5, 0.5))
dev.off()

save.image("umap/umapSession.RData")







umapSweep

for(b in c(1,10,100,10000)){
  
  for(nnb in c(10, 100,1000)){
    
    for(lc in c(1,10,100)){
      
      for(ns in c(5)){
        
        for(rs in c(0.01,0.1,1,10)){
          
          for(lr in c(1,10)){
            
            sampleUmap <- umap(X = NULL, nn_method = nn, init = init, 
                               bandwidth = b,n_epochs = 1000,learning_rate = lr,
                               n_neighbors = nnb,local_connectivity = lc, repulsion_strength = rs, negative_sample_rate = ns )#, n_neighbors = 50,n_epochs = 1000,learning_rate = 10,local_connectivity = 50,bandwidth = 1000, repulsion_strength = 0.1
            rownames(sampleUmap) <- rownames(combinedHealtyTissuePca)
            colnames(sampleUmap) <- c("UMAP1", "UMAP2")
            
            
            pngFile <- paste0("umapSweep/", "b" , b, "nnb" , nnb, "lc" , lc, "ns", ns, "rs" , rs, "lr" , lr, ".png")
            
            if(!file.exists(pngFile)){
              
              png(file = pngFile, height = 1000, width = 1000)
              plot(sampleUmap[,1],sampleUmap[,2], pch = 16, col=adjustcolor(samplesWithPredictionNoOutliers$col, alpha.f = 0.5), bty="n", xlab = "UMAP-1", ylab = "UMAP-2", cex = 0.7)
              dev.off()
              
              
            }
            
            
            
            
          }
          
        }
        
      }
      
      
    }
    
  }
  
}




colnamesToUpdate <- colnames(pcsAndMeta)[colnames(pcsAndMeta) %in% colnames(combinedMeta)]
all(rownames(pcsAndMeta) %in% rownames(combinedMeta))
pcsAndMeta[,colnamesToUpdate] <- combinedMeta[rownames(pcsAndMeta),colnamesToUpdate]


table(pcsAndMeta$selectedSamples, useNA = "a")


clusterAnnotations <- read.delim("umap/annotationsBasedOnOldUmap.txt", row.names = 1)
samplesWithClusterAnnotation <- rownames(pcsAndMeta)[rownames(pcsAndMeta) %in% rownames(clusterAnnotations)]

pcsAndMeta$ClusterAnnotation <- NA
pcsAndMeta[samplesWithClusterAnnotation, "ClusterAnnotation"] <- clusterAnnotations[samplesWithClusterAnnotation,"ClusterAnnotation"]
table(pcsAndMeta$ClusterAnnotation, useNA = "a")

tissueSamples <- pcsAndMeta[pcsAndMeta$selectedSamples,]

tissueSamples$class <- tissueSamples$Tissue


hasT2 <- tissueSamples$Tissue2 != ""
tissueSamples$class[hasT2] <- paste0(tissueSamples$Tissue[hasT2], "-", tissueSamples$Tissue2[hasT2])
table(tissueSamples$class)
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
tissueSamples$class[!tissueSamples$class %in% mapping$Class]

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
  min_dist = 1, init_sdev = 1e-4, learning_rate = 2, 
  spread = 20, 
  bandwidth = 10,
  scale = "scale",
  local_connectivity = 10,
  repulsion_strength = 0.5,
  metric = "correlation")


rownames(sampleUmap) <- rownames(tissueSamples)
colnames(sampleUmap) <- c("UMAP1", "UMAP2")
save(sampleUmap, file = "umap/sampleUmap6.RData")

#load(file = "umap/sampleUmap3.RData")




umapAndMeta <- merge(sampleUmap, tissueSamples, by = 0)
rownames(umapAndMeta) <- umapAndMeta$Row.names
dim(umapAndMeta)




rpng()

par(mar = c(3,5,0.1,0.1), xpd = NA)
plot(umapAndMeta[umapAndMeta$plotOrderTissues,"UMAP1"], umapAndMeta[umapAndMeta$plotOrderTissues,"UMAP2"], col = umapAndMeta$TissueCol[umapAndMeta$plotOrderTissues], cex = 0.2, pch = 16)

dev.off()

plot(umapAndMeta[umapAndMeta$plotOrderTissues,"UMAP1"], umapAndMeta[umapAndMeta$plotOrderTissues,"UMAP2"], col = umapAndMeta$TissueCol[umapAndMeta$plotOrderTissues], cex = 0.2, pch = 16, xlim = c(-100,100), ylim = c(-100,100))

  

locator(n =2, type = "l")
cluster1 <- locator(n =2, type = "l")
cluster2 <- locator(n =2, type = "l")


write.table(umapAndMeta[,!grepl("PC_",colnames(umapAndMeta))],file = "umaptest.txt", sep = "\t", quote = F, col.names = NA)
#save(umapAndMeta, file = "umaptest.RData")
#load("umaptest.RData")

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
plot(umapAndMeta[umapAndMeta$plotOrderTissues,"UMAP1"], umapAndMeta[umapAndMeta$plotOrderTissues,"UMAP2"], col = umapAndMeta$TissueCol[umapAndMeta$plotOrderTissues], cex = 0.2, pch = 16)

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

#Take only samples that have an annotation
umapAndMetaClassified <- umapAndMeta[!is.na(umapAndMeta$umapFactor),]
#First put all in test, algorithm will put some 
umapAndMetaClassified$training <- FALSE

tissueClass <- levels(umapAndMetaClassified$umapFactor)[2]
study <- "GTEx"

set.seed(42)
#for each tissue select samples for training
for(tissueClass in levels(umapAndMetaClassified$umapFactor)){
  thisTissueSamples <- umapAndMetaClassified$umapFactor==tissueClass
  studiesForThisTissue <- unique(umapAndMetaClassified$study[thisTissueSamples])
  numberOfStudies <- length(studiesForThisTissue)
  numberOfSamplesPerStudy <- ceiling(minSamplesTraining / numberOfStudies)
  print(paste(tissueClass, length(studiesForThisTissue), numberOfSamplesPerStudy, sep = " - "))
  #for each study put samples to training or test
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


library(glmnet)
cfit <- cv.glmnet(x = as.matrix(umapAndMetaClassifiedTraining[,paste0("PC_",1:compsToUse)]), y = umapAndMetaClassifiedTraining$umapFactor, family = "multinomial", type.measure = "class")
cfit

rpng()
plot(cfit) 
dev.off()



assess.glmnet(cfit, newx = as.matrix(umapAndMetaClassifiedTest[,paste0("PC_",1:compsToUse)]), newy = umapAndMetaClassifiedTest$umapFactor, family = "multinomial", type.measure = "class", keep = TRUE, alpha=1, lambda = "1se")



predictionsTest <- predict(cfit, s = "lambda.1se", newx = as.matrix(umapAndMetaClassifiedTest[,paste0("PC_",1:compsToUse)]), type = "class")

predictionsTestScores <- predict(cfit, s = "lambda.1se", newx = as.matrix(umapAndMetaClassifiedTest[,paste0("PC_",1:compsToUse)]), type = "response")
predictionsTestScores <- predictionsTestScores[,,1]
umapAndMetaClassifiedTest$predictedTissueScore <- apply(predictionsTestScores, 1, max)

prop = 0.5

predictionsInTest <- sapply(seq(0,1,0.05), function(prop){  

  umapAndMetaClassifiedTest$predictedTissue <- predictionsTest[,1]
  
  
  umapAndMetaClassifiedTest$predictedTissue[umapAndMetaClassifiedTest$predictedTissueScore <= prop] <- NA
  
  umapAndMetaClassifiedTest$misclasified <- FALSE
  umapAndMetaClassifiedTest$misclasified[!is.na(umapAndMetaClassifiedTest$umapFactor) & !is.na(umapAndMetaClassifiedTest$predictedTissue) &  umapAndMetaClassifiedTest$umapFactor != umapAndMetaClassifiedTest$predictedTissue] <- TRUE
  errors <- sum(umapAndMetaClassifiedTest$misclasified )
  
  umapAndMetaClassifiedTest$notPredictedBack <- FALSE
  umapAndMetaClassifiedTest$notPredictedBack[!is.na(umapAndMetaClassifiedTest$umapFactor) & is.na(umapAndMetaClassifiedTest$predictedTissue) ] <- TRUE
  missed <- sum(umapAndMetaClassifiedTest$notPredictedBack)
  
  total <- nrow(umapAndMetaClassifiedTest)
  
  missedPercentage <- missed / total
  errorPercentage <- errors / total
  
  return(c("Threshold" = prop, "MissedPerc" = missedPercentage , "ErrorPerc" = errorPercentage ))
  
})
predictionsInTest

tissueClass <- levels(umapAndMetaClassified$umapFactor)[1]

predictionsInTestPerTissue <- lapply(levels(umapAndMetaClassified$umapFactor), function(tissueClass){
  predictionsInTestThisTissue <- sapply(seq(0,1,0.05), function(prop){  
    
    umapAndMetaClassifiedTestTissue <- umapAndMetaClassifiedTest[umapAndMetaClassifiedTest$umapFactor == tissueClass,]
    umapAndMetaClassifiedTestTissue$predictedTissue <- predictionsTest[umapAndMetaClassifiedTest$umapFactor == tissueClass,1]
    
    
    umapAndMetaClassifiedTestTissue$predictedTissue[umapAndMetaClassifiedTestTissue$predictedTissueScore <= prop] <- NA
    
    umapAndMetaClassifiedTestTissue$misclasified <- FALSE
    umapAndMetaClassifiedTestTissue$misclasified[!is.na(umapAndMetaClassifiedTestTissue$umapFactor) & !is.na(umapAndMetaClassifiedTestTissue$predictedTissue) &  umapAndMetaClassifiedTestTissue$umapFactor != umapAndMetaClassifiedTestTissue$predictedTissue] <- TRUE
    errors <- sum(umapAndMetaClassifiedTestTissue$misclasified )
    
    umapAndMetaClassifiedTestTissue$notPredictedBack <- FALSE
    umapAndMetaClassifiedTestTissue$notPredictedBack[!is.na(umapAndMetaClassifiedTestTissue$umapFactor) & is.na(umapAndMetaClassifiedTestTissue$predictedTissue) ] <- TRUE
    missed <- sum(umapAndMetaClassifiedTestTissue$notPredictedBack)
    
    total <- nrow(umapAndMetaClassifiedTestTissue)
    
    missedPercentage <- missed / total
    errorPercentage <- errors / total
    
    return(c("Threshold" = prop, "MissedPerc" = missedPercentage , "ErrorPerc" = errorPercentage ))
    
  })
  return(predictionsInTestThisTissue)
})
names(predictionsInTestPerTissue) <- levels(umapAndMetaClassified$umapFactor)
str(predictionsInTestPerTissue)

x <- sapply(predictionsInTestPerTissue, function(predictionsInTestThisTissue){
  return(predictionsInTestThisTissue[3,11])
})
sort(x)

predictionsInTest[2,15]

predictionsInTestPerTissue[["Whole Blood Fetal"]]

layout(matrix(1:2, nrow = 1))
plot(t(predictionsInTest[1:2,]), main = "Percentage classification missed in test dataset")
for(tissueClass in levels(umapAndMetaClassified$umapFactor)){
  predictionsInTestThisTissue <- predictionsInTestPerTissue[[tissueClass]]
  points(t(predictionsInTestThisTissue[1:2,]), type = "l",  col=adjustcolor("grey", alpha.f = 0.5))
}
plot(t(predictionsInTest[c(1,3),]), main = "Percentage wrong classification in test dataset")
sink <- sapply(predictionsInTestPerTissue, function(predictionsInTestThisTissue){
  points(t(predictionsInTestThisTissue[c(1,3),]), type = "l",  col=adjustcolor("grey", alpha.f = 0.5))
})




confusion <- confusion.glmnet(cfit, newx = as.matrix(umapAndMetaClassifiedTest[,paste0("PC_",1:compsToUse)]), newy = umapAndMetaClassifiedTest$umapFactor, family = "multinomial", type.measure = "class", keep = TRUE, alpha=1, lambda = "1se")
diag(confusion) <- 0

library(heatmap3)

rpng()
pdf("confusion.pdf", width = 12, height = 12)
heatmap3(confusion, Rowv = NA, Colv = NA, balanceColor =T, scale = "none")
dev.off()


predictions <- predict(cfit, s = "lambda.1se", newx = as.matrix(umapAndMeta[,paste0("PC_",1:compsToUse)]), type = "class")
umapAndMeta$predictedTissue <- predictions[,1]

predictionsScores <- predict(cfit, s = "lambda.1se", newx = as.matrix(umapAndMeta[,paste0("PC_",1:compsToUse)]), type = "response")
predictionsScores <- predictionsScores[,,1]
rownames(predictionsScores) <- umapAndMeta$Row.names
umapAndMeta$predictedTissueScore <- apply(predictionsScores, 1, max)

sum(umapAndMeta$predictedTissueScore <= 0.5)
umapAndMeta$predictedTissue[umapAndMeta$predictedTissueScore <= 0.5] <- NA


rpng()
hist(umapAndMeta$predictedTissueScore)
dev.off()

umapAndMeta$misclasified <- FALSE
umapAndMeta$misclasified[!is.na(umapAndMeta$umapFactor) & !is.na(umapAndMeta$predictedTissue) &  umapAndMeta$umapFactor != umapAndMeta$predictedTissue] <- TRUE
sum(umapAndMeta$misclasified )

umapAndMeta$notPredictedBack <- FALSE
umapAndMeta$notPredictedBack[!is.na(umapAndMeta$umapFactor) & is.na(umapAndMeta$predictedTissue) ] <- TRUE
sum(umapAndMeta$notPredictedBack)

sum(!is.na(umapAndMeta$predictedTissue))

length(unique((umapAndMeta$predictedTissue)))

sum(table((umapAndMeta$predictedTissue)) >= 1000)
hist(table((umapAndMeta$predictedTissue)), breaks =25)
barplot(table((umapAndMeta$predictedTissue)))

sort(table(umapAndMeta[umapAndMeta$misclasified, "umapFactor"]))
sort(table(umapAndMeta[umapAndMeta$notPredictedBack, "umapFactor"]))

tissueClass <- levels(umapAndMeta$umapFactor)[1]

pdf("tissuePrediction.pdf")
for(tissueClass in levels(umapAndMeta$umapFactor)){
  
  umapAndMeta$ThisTissueCol <- defaultCol
  umapAndMeta$ThisTissueCol[!is.na(umapAndMeta$umapFactor) & tissueClass == umapAndMeta$umapFactor & !umapAndMeta$misclasified] <- adjustcolor("forestgreen", alpha.f = 0.5)
  umapAndMeta$ThisTissueCol[!is.na(umapAndMeta$umapFactor) & tissueClass == umapAndMeta$umapFactor & umapAndMeta$notPredictedBack] <- adjustcolor("hotpink", alpha.f = 0.5)
  umapAndMeta$ThisTissueCol[!is.na(umapAndMeta$umapFactor) & tissueClass == umapAndMeta$umapFactor & umapAndMeta$misclasified] <- adjustcolor("violetred3", alpha.f = 0.5)
  umapAndMeta$ThisTissueCol[!is.na(umapAndMeta$umapFactor) & !is.na(umapAndMeta$predictedTissue) & tissueClass != umapAndMeta$umapFactor & tissueClass == umapAndMeta$predictedTissue] <- adjustcolor("orange1", alpha.f = 0.5)
  umapAndMeta$ThisTissueCol[is.na(umapAndMeta$umapFactor) & !is.na(umapAndMeta$predictedTissue) & tissueClass == umapAndMeta$predictedTissue] <- adjustcolor("dodgerblue1", alpha.f = 0.5)
    
  predictedBack <- sum(!is.na(umapAndMeta$umapFactor) & tissueClass == umapAndMeta$umapFactor & !umapAndMeta$misclasified)
  notPredictedBack <- sum(!is.na(umapAndMeta$umapFactor) & tissueClass == umapAndMeta$umapFactor & umapAndMeta$notPredictedBack) 
  predictedAsOther <-  sum(!is.na(umapAndMeta$umapFactor) & tissueClass == umapAndMeta$umapFactor & umapAndMeta$misclasified)
  otherPredicted <- sum(!is.na(umapAndMeta$umapFactor) & !is.na(umapAndMeta$predictedTissue) & tissueClass != umapAndMeta$umapFactor & tissueClass == umapAndMeta$predictedTissue)
  newPredicted <- sum(is.na(umapAndMeta$umapFactor) & !is.na(umapAndMeta$predictedTissue) & tissueClass == umapAndMeta$predictedTissue)
  
  table(umapAndMeta$ThisTissueCol, useNA = "a")
  
  umapAndMeta$plotOrderThisTissues <- order(umapAndMeta$ThisTissueCol != defaultCol)
  
  #rpng()
  layout(matrix(c(1,2,3), ncol = 1, byrow = T), heights = c(0.05,0.85,0.1))
  par(mar = c(0,0,0,0), xpd = NA)
  plot.new()
  plot.window(xlim = 0:1, ylim = 0:1)
  text(0.5,0.5,tissueClass, cex = 2 , font = 2)
  
  par(mar = c(5,5,0,0.1), xpd = NA)
  plot(umapAndMeta[umapAndMeta$plotOrderThisTissues,"UMAP1"], umapAndMeta[umapAndMeta$plotOrderThisTissues,"UMAP2"], col = umapAndMeta$ThisTissueCol[umapAndMeta$plotOrderThisTissues], cex = 0.2, pch = 16, bty="n", xlab = "UMAP-1", ylab = "UMAP-2")
  
  par(mar = c(0,0,0,0), xpd = NA)
  plot.new()
  plot.window(xlim = 0:1, ylim = 0:1)
  legend("center", fill = c(
        "forestgreen", 
        "hotpink",
        "violetred3",
        "orange1",
        "dodgerblue1"
      ),
    legend = c(
        paste0(tissueClass, " correctly predicted back (", predictedBack,")"),
        paste0(tissueClass, " not predicted back (", notPredictedBack,")"),
        paste0(tissueClass, " predicted as other (", predictedAsOther,")"),
        paste0("Other tissue predicted as ", tissueClass," (", otherPredicted,")"),
        paste0("Unkown predicted as ", tissueClass, " (", newPredicted,")")
      ), 
    bty = "n")
  
  
  
 #dev.off()
  
}
dev.off()

#save(umapAndMeta, file = "tissuePredictions/tissuePredictions_16_09_22.RData")
load("tissuePredictions/tissuePredictions_16_09_22.RData", verbose = T)

unique(umapAndMeta$predictedTissue)[!unique(umapAndMeta$predictedTissue) %in% rownames(tissueCol)]



clusterToExclude <- c("U2-OS", "Leukemia_blood-cell-line", "HAP1", "LNCaP")




samplesWithPrediction <- umapAndMeta[!is.na(umapAndMeta$predictedTissue) & !umapAndMeta$predictedTissue %in% clusterToExclude, c(
  "predictedTissue",
  "predictedTissueScore",
  "umapFactor",
  "misclasified",
  "study",
  "sra.library_layout"
)]
colnames(samplesWithPrediction)[3] <- "annotatedTissue"
str(samplesWithPrediction)
#save(samplesWithPrediction, file = "tissuePredictions/samplesWithPrediction_16_09_22.RData")

write.table(samplesWithPrediction, file = "samplesWithPrediction.txt")
load("tissuePredictions/samplesWithPrediction_16_09_22.RData")
str(samplesWithPrediction)
table(samplesWithPrediction$predictedTissue)

load(file = "umap/sampleUmap6.RData", verbose = T)


umapAndPredictions <- merge(samplesWithPrediction, sampleUmap, by = 0 )
rownames(umapAndPredictions) <- umapAndPredictions$Row.names


umapAndPredictions$TissuePredictedCol <- defaultCol
umapAndPredictions$TissuePredictedCol[umapAndPredictions$predictedTissue %in% rownames(tissueCol)] <- adjustcolor(tissueCol[as.character(umapAndPredictions$predictedTissue[umapAndPredictions$predictedTissue %in% rownames(tissueCol)]),1], alpha.f = 0.5)
umapAndPredictions$plotOrderTissuePredicted <- order(umapAndPredictions$TissuePredictedCol != defaultCol)

#rpng()

par(mar = c(3,3,0.1,0.1), xpd = NA)
plot(umapAndPredictions[umapAndPredictions$plotOrderTissuePredicted,"UMAP1"], umapAndPredictions[umapAndPredictions$plotOrderTissuePredicted,"UMAP2"], col = umapAndPredictions$TissuePredictedCol[umapAndPredictions$plotOrderTissuePredicted], cex = 0.2, pch = 16)

plot(umapAndPredictions[umapAndPredictions$plotOrderTissuePredicted,"UMAP1"], umapAndPredictions[umapAndPredictions$plotOrderTissuePredicted,"UMAP2"], col = umapAndPredictions$TissuePredictedCol[umapAndPredictions$plotOrderTissuePredicted], cex = 0.2, pch = 16, xlim = c(-100,70), ylim = c(-50,50))




#dev.off()

locator(n =2, type = "l")


pdf(file = "umapPredicted.pdf", width = 16, height = 8)
#rpng()

layout(matrix(1:2,ncol = 2))

par(mar = c(5,5,0.1,0.1), xpd = NA)
plot(umapAndPredictions[umapAndPredictions$plotOrderTissuePredicted,"UMAP1"], umapAndPredictions[umapAndPredictions$plotOrderTissuePredicted,"UMAP2"], col = umapAndPredictions$TissuePredictedCol[umapAndPredictions$plotOrderTissuePredicted], cex = 0.2, pch = 16, bty = "n", xlab = "UMAP-1", ylab = "UMAP-2")

par(mar = c(0,0,0,0), xpd = NA)
plot.new()
plot.window(xlim = 0:1, ylim = 0:1)
legend("center", fill = tissueCol[rownames(tissueCol) %in% umapAndPredictions$predictedTissue,1], legend = row.names(tissueCol)[rownames(tissueCol) %in% umapAndPredictions$predictedTissue], bty = "n", ncol = 2,cex = 0.7)


dev.off()




countTable <- table(umapAndPredictions$predictedTissue)
sum(countTable)
sum(countTable >= 500)
pdf("baplotTissues.pdf", width = 15, height = 10)
par(mar = c(25,5,2,0.1), xpd = NA)
b <- barplot(countTable, las =2, col = tissueCol[names(countTable),])
text(b, countTable + 280, countTable, font=1, srt = 90)
dev.off()
