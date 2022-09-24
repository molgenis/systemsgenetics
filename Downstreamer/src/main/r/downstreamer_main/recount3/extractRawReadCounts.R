#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, password = "laberkak", sync = T)


remoter::client("localhost", port = 55501, password = "laberkak")



setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")



sraFiles <- list.files(path="rse-sra/SRA_Files/", pattern="sra*", full.names=TRUE, recursive=FALSE)
gtexFiles <- list.files(path="rse-gtex/rse_gtex", pattern="rse*", full.names=TRUE, recursive=FALSE)
allFiles <- c(sraFiles, gtexFiles, "rse-tcga/rseTCGA.rda")

load("tissuePredictions/samplesWithPrediction_16_09_22.RData")
selectedSamples <- rownames(samplesWithPrediction)
str(selectedSamples)


#file = sraFiles[10]



perChunkExp <- sapply(allFiles, function(file){
  
  loadedObject <- load(file)
  
  sreObjects <- get(loadedObject[1])
  
  #sreObject <- sreObject[[1]]
  
  perStudyExp <- lapply(sreObjects, function(sreObject){
    studyExp <- sreObject@assays@data@listData$raw_counts
    return(studyExp[,colnames(studyExp) %in% selectedSamples, drop = F])
  })
  
  return(do.call(cbind, perStudyExp))
  
})

selectedSamplesExp <- do.call(cbind, perChunkExp)

all(selectedSamples %in% colnames(selectedSamplesExp ))
table(selectedSamples %in% colnames(selectedSamplesExp ))

#Some samples are duplicated in the chunks, now make sure only one is in the matrix
uniqueSamplesIndex <- match(selectedSamples, colnames(selectedSamplesExp))
selectedSamplesExp <- selectedSamplesExp[,uniqueSamplesIndex]

save(selectedSamplesExp, file = "tissuePredictions/selectedSamplesRawExpression.RData")

tissueClasses <- levels(samplesWithPrediction$predictedTissue)

tissue <- tissueClasses[1]

rownames(samplesWithPrediction)[!rownames(samplesWithPrediction) %in% colnames(selectedSamplesExp)]
table(rownames(samplesWithPrediction) %in% colnames(selectedSamplesExp))
dim(samplesWithPrediction)
dim(selectedSamplesExp)

perTissueExp <- lapply(tissueClasses, function(tissue){
  tissueSamples <- rownames(samplesWithPrediction)[samplesWithPrediction$predictedTissue == tissue]
  tissueExp <- selectedSamplesExp[,tissueSamples]
})

