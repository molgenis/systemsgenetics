#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55001, password = "laberkak", sync = T)


remoter::client("localhost", port = 55001, password = "laberkak")
plot(2)
dev.off()

library(uwot)

setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs\\Recount3\\")
setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")

load("tissuePredictions/samplesWithPrediction_16_09_22_noOutliers.RData", verbose = T)

load(file = "DataForPredictions.RData")

colnames(pcsAndMeta)
pcsAndMeta <- pcsAndMeta[!pcsAndMeta$exclude,]


str(samplesWithPredictionNoOutliers)

dim(pcsAndMeta)
table(pcsAndMeta$selectedSamples, useNA = "always")
colnames(pcsAndMeta)
colnames(samplesWithPredictionNoOutliers)


head(rownames(samplesWithPredictionNoOutliers))
head(pcsAndMeta$Row.names)


str(pcsAndMeta$Row.names[pcsAndMeta$selectedSamples])


all(rownames(samplesWithPredictionNoOutliers) %in%  pcsAndMeta$Row.names[pcsAndMeta$selectedSamples])

rownames(samplesWithPredictionNoOutliers[,"predictedTissue", drop =F ])

pcsAndMeta$cancerTraining <- NA
cancerCelllineTissuePred <- merge(pcsAndMeta, samplesWithPredictionNoOutliers[,"predictedTissue", drop =F], all = T, by.x = 1, by.y = 0 )

dim(pcsAndMeta)
colnames(test)

table(cancerCelllineTissuePred$excludeBasedOnPredictionCellline2, !is.na(cancerCelllineTissuePred$predictedTissue), useNA = "a")
table(cancerCelllineTissuePred$excludeBasedOnPredictionCancer, !is.na(cancerCelllineTissuePred$predictedTissue), useNA = "a")
table(cancerCelllineTissuePred$excludeBasedOnPredictionCancer, cancerCelllineTissuePred$excludeBasedOnPredictionCellline2, useNA = "a")

table(cancerCelllineTissuePred$excludeBasedOnPredictionCellline2 | cancerCelllineTissuePred$excludeBasedOnPredictionCancer, !is.na(cancerCelllineTissuePred$predictedTissue), useNA = "a")

cancerCelllineTissuePred$predictedCellineCancer <- cancerCelllineTissuePred$excludeBasedOnPredictionCellline2 | cancerCelllineTissuePred$excludeBasedOnPredictionCancer

table(cancerCelllineTissuePred$predictedCellineCancer, !is.na(cancerCelllineTissuePred$predictedTissue), useNA = "a")

#only select samples predicted as cancer / celline and main tissue samples that passed the QC per tissue
cancerCelllineTissuePred <- cancerCelllineTissuePred[cancerCelllineTissuePred$predictedCellineCancer | !is.na(cancerCelllineTissuePred$predictedTissue), ]

dim(cancerCelllineTissuePred)



cancerCelllineTissuePred$predictedTissue <- as.factor(cancerCelllineTissuePred$predictedTissue)




minSamplesTraining <- 50
maxFractionOfStudy <- 0.8

#Take only samples that have an annotation
cancerCelllineTissuePredClassified <- cancerCelllineTissuePred[!is.na(cancerCelllineTissuePred$predictedTissue),]
#First put all in test, algorithm will put some 
cancerCelllineTissuePredClassified$training <- FALSE

tissueClass <- levels(cancerCelllineTissuePredClassified$predictedTissue)[1]
study <- "GTEx"

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

rpng()
plot(cfit) 
rpng.off()
#dev.off()



assess.glmnet(cfit, newx = as.matrix(cancerCelllineTissuePredClassifiedTest[,paste0("PC_",1:compsToUse)]), newy = cancerCelllineTissuePredClassifiedTest$umapFactor, family = "multinomial", type.measure = "class", keep = TRUE, alpha=1, lambda = "1se")



predictionsTest <- predict(cfit, s = "lambda.1se", newx = as.matrix(cancerCelllineTissuePredClassifiedTest[,paste0("PC_",1:compsToUse)]), type = "class")

