#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55001, password = "laberkak", sync = T)


#remoter::client("localhost", port = 55001, password = "laberkak")

library(uwot)

setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Recount3\\")
setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

load("tissuePredictions/samplesWithPrediction_16_09_22_noOutliers.RData", verbose = T)

load(file = "DataForPredictions.RData")

pcsAndMeta <- pcsAndMeta[!pcsAndMeta$exclude,]


all(rownames(samplesWithPredictionNoOutliers) %in%  pcsAndMeta$Row.names[pcsAndMeta$selectedSamples])


table(cancerCelllineTissuePred$Cancer, useNA = "a")

pcsAndMeta$cancerTraining <- NA
cancerCelllineTissuePred <- merge(pcsAndMeta, samplesWithPredictionNoOutliers[,"predictedTissue", drop =F], all = T, by.x = 1, by.y = 0 )

cancerCelllineTissuePred$predictedCellineCancer <- cancerCelllineTissuePred$excludeBasedOnPredictionCellline2 | cancerCelllineTissuePred$excludeBasedOnPredictionCancer

table(cancerCelllineTissuePred$predictedCellineCancer)

cancerCelllineTissuePred$AnnotatedAsTissue <- (cancerCelllineTissuePred$Tissue != "" | cancerCelllineTissuePred$Tissue2 != "") & !cancerCelllineTissuePred$Cancer & !cancerCelllineTissuePred$Cellline
table(cancerCelllineTissuePred$AnnotatedAsTissue)

table(cancerCelllineTissuePred$predictedCellineCancer, cancerCelllineTissuePred$AnnotatedAsTissue)


table(cancerCelllineTissuePred$predictedCellineCancer, is.na(cancerCelllineTissuePred$predictedTissue), useNA = "a")

#only select samples predicted as cancer / celline and main tissue samples that passed the QC per tissue
cancerCelllineTissuePred <- cancerCelllineTissuePred[cancerCelllineTissuePred$predictedCellineCancer | !is.na(cancerCelllineTissuePred$predictedTissue), ]

dim(cancerCelllineTissuePred)



cancerCelllineTissuePred$predictedTissue <- as.factor(cancerCelllineTissuePred$predictedTissue)

bar



minSamplesTraining <- 50
maxFractionOfStudy <- 0.8

#Take only samples that have an annotation
cancerCelllineTissuePredClassified <- cancerCelllineTissuePred[!is.na(cancerCelllineTissuePred$predictedTissue),]
#First put all in test, algorithm will put some 
cancerCelllineTissuePredClassified$training <- FALSE

tissueClass <- levels(cancerCelllineTissuePredClassified$predictedTissue)[1]
#study <- "GTEx"

set.seed(42)
#for each tissue slecect samples for training
for(tissueClass in levels(cancerCelllineTissuePredClassified$predictedTissue)){
  thisTissueSamples <- cancerCelllineTissuePredClassified$predictedTissue==tissueClass
  studiesForThisTissue <- unique(cancerCelllineTissuePredClassified$study[thisTissueSamples])
  numberOfStudies <- length(studiesForThisTissue)
  numberOfSamplesPerStudy <- ceiling(minSamplesTraining / numberOfStudies)
  print(paste(tissueClass, length(studiesForThisTissue), numberOfSamplesPerStudy, sep = " - "))
  #for each study put samples to training or test
  for(study in studiesForThisTissue){
    
    thisTissueAndStudySamples <- thisTissueSamples & cancerCelllineTissuePredClassified$study == study
    thisTissueAndStudySamplesCount <- sum(thisTissueAndStudySamples)
    
    #Don't select more samples from study then the study has and also no more then set fraction. Do floor to put studies with single sample to testset
    potentialMax <- floor(thisTissueAndStudySamplesCount * maxFractionOfStudy)
    numberTrainingSamplesThisStudy <- if(potentialMax > numberOfSamplesPerStudy) numberOfSamplesPerStudy  else potentialMax
    if(numberTrainingSamplesThisStudy > 0){
      #The which will get all indices for the samples of this study-tissue combination. These are then samples for the samples used for training
      trainingSamplesThisStudy <- sample(which(thisTissueAndStudySamples), numberTrainingSamplesThisStudy)
      #Set selected to TRUE
      cancerCelllineTissuePredClassified$training[trainingSamplesThisStudy] <- TRUE
    }
    
    
    #print(paste0(thisTissueAndStudySamplesCount, " - ", numberTrainingSamplesThisStudy))
  }
  
}

sum(cancerCelllineTissuePredClassified$training)


cancerCelllineTissuePredClassifiedTraining <- cancerCelllineTissuePredClassified[cancerCelllineTissuePredClassified$training,]
table(cancerCelllineTissuePredClassifiedTraining$predictedTissue)
cancerCelllineTissuePredClassifiedTest <- cancerCelllineTissuePredClassified[!cancerCelllineTissuePredClassified$training,]
dim(cancerCelllineTissuePredClassifiedTest)



library(glmnet)
cfit <- cv.glmnet(x = as.matrix(cancerCelllineTissuePredClassifiedTraining[,paste0("PC_",1:compsToUse)]), y = cancerCelllineTissuePredClassifiedTraining$predictedTissue, family = "multinomial", type.measure = "class")
cfit

plot(cfit) 



assess.glmnet(cfit, newx = as.matrix(cancerCelllineTissuePredClassifiedTest[,paste0("PC_",1:compsToUse)]), newy = cancerCelllineTissuePredClassifiedTest$predictedTissue, family = "multinomial", type.measure = "class", keep = TRUE, alpha=1, lambda = "1se")



predictionsTest <- predict(cfit, s = "lambda.1se", newx = as.matrix(cancerCelllineTissuePredClassifiedTest[,paste0("PC_",1:compsToUse)]), type = "class")

str(predictionsTest)


sum(!cancerCelllineTissuePredClassifiedTest$predictedTissue == predictionsTest[,1])


predictions <- predict(cfit, s = "lambda.1se", newx = as.matrix(cancerCelllineTissuePred[,paste0("PC_",1:compsToUse)]), type = "class")

table(predictions)
cancerCelllineTissuePred$newPrediction <- predictions[,1]

predictionsScores <- predict(cfit, s = "lambda.1se", newx = as.matrix(cancerCelllineTissuePred[,paste0("PC_",1:compsToUse)]), type = "response")
predictionsScores <- predictionsScores[,,1]
cancerCelllineTissuePred$newPredictionScore <- apply(predictionsScores, 1, max)


hist(cancerCelllineTissuePred$newPredictionScore, breaks = 42)

layout(matrix(1:2, ncol = 2))
par(mar = c(5, 4, 4, 2) + 0.1)
hist(cancerCelllineTissuePred$newPredictionScore[!is.na(cancerCelllineTissuePred$predictedTissue)], breaks = 42, main = "Posterior probability of predicited tissue of normal", xlab = "Posterior probability")
hist(cancerCelllineTissuePred$newPredictionScore[cancerCelllineTissuePred$predictedCellineCancer], breaks = 42, main = "Posterior probability of predicited tissue of cancer/cell-line", xlab = "Posterior probability")

sum(cancerCelllineTissuePred$newPredictionScore >= 0.5)
#umapAndMeta$predictedTissue[umapAndMeta$predictedTissueScore <= 0.5] <- NA



#Subset to only the samples that are predited to be cancer or celllines
cancerCelllineTissuePred2 <- cancerCelllineTissuePred[cancerCelllineTissuePred$predictedCellineCancer,]
dim(cancerCelllineTissuePred2)

predictionCounts <- table(cancerCelllineTissuePred2$newPrediction[cancerCelllineTissuePred2$newPredictionScore >= 0.5])
layout(1)
par(mar  =c(3,30,1,1))
barplot(sort(predictionCounts, decreasing = F), horiz = T, las =1 )

table(cancerCelllineTissuePred2$Cohort)
aggregate(cancerCelllineTissuePred2$newPredictionScore, list(cancerCelllineTissuePred2$Cohort), mean)

hist(cancerCelllineTissuePred2$newPredictionScore[cancerCelllineTissuePred2$Cohort=="TCGA"])


sort(table(cancerCelllineTissuePred2$newPrediction[cancerCelllineTissuePred2$Cohort=="GTEx"]))

length(cancerCelllineTissuePred2$gtex.smtsd[cancerCelllineTissuePred2$Cohort=="GTEx" & (cancerCelllineTissuePred2$Tissue != "" | cancerCelllineTissuePred2$Tissue2 != "")])
barplot(sort(predictionCounts, decreasing = F), horiz = T, las =1 )



sort(table(paste0(cancerCelllineTissuePred2$newPrediction, "-", cancerCelllineTissuePred2$CelllineName)[cancerCelllineTissuePred2$CelllineName=="HeLa"]))


table(cancerCelllineTissuePred2$CelllineName[cancerCelllineTissuePred2$newPrediction=="hematopoietic progenitors" & cancerCelllineTissuePred2$newPredictionScore >= 0.5])
table(cancerCelllineTissuePred2$CelllineName[cancerCelllineTissuePred2$newPrediction=="fibroblasts_cell-lines_smooth-muscle-cell_mesenchymal-stem-cells" & cancerCelllineTissuePred2$newPredictionScore >= 0.5])
table(cancerCelllineTissuePred2$CelllineName[cancerCelllineTissuePred2$newPrediction=="T-Cells" & cancerCelllineTissuePred2$newPredictionScore >= 0.5])

"tcga.gdc_cases.project.primary_site" 
"Cohort"   
study
Tissue
Tissue2
CelllineName

cancerCelllineTissuePred2Tcga <- cancerCelllineTissuePred2[cancerCelllineTissuePred2$Cohort=="TCGA" & cancerCelllineTissuePred2$tcga.cgc_sample_sample_type=="Primary Tumor" & cancerCelllineTissuePred2$newPredictionScore >= 0.5,]
table(cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site, useNA = "a")

sort(table(paste0(cancerCelllineTissuePred2Tcga$newPrediction, " - ", cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site)[]))

table(cancerCelllineTissuePred2Tcga$newPrediction[cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site == "Liver"])
table(cancerCelllineTissuePred2Tcga$newPrediction[cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site == "Kidney"])
table(cancerCelllineTissuePred2Tcga$newPrediction[cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site == "Prostate"])
table(cancerCelllineTissuePred2Tcga$newPrediction[cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site == "Lung"])
table(cancerCelllineTissuePred2Tcga$newPrediction[cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site == "Brain"])
table(cancerCelllineTissuePred2Tcga$newPrediction[cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site == "Breast"])
table(cancerCelllineTissuePred2Tcga$newPrediction[cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site == "Colorectal"])
table(cancerCelllineTissuePred2Tcga$newPrediction[cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site == "Head and Neck"])
table(cancerCelllineTissuePred2Tcga$newPrediction[cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site == "Thyroid"])


sort(table(cancerCelllineTissuePred2Tcga$tcga.gdc_cases.project.primary_site))
