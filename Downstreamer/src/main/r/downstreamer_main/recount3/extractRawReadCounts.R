#srun --cpus-per-task=1 --mem=50gb --nodes=1 --qos=priority --time=168:00:00 --pty bash -i
#remoter::server(verbose = T, port = 55556, sync = T)


remoter::client("localhost", port = 55501)



setwd("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/")



sraFiles <- list.files(path="rse-sra/SRA_Files/", pattern="sra*", full.names=TRUE, recursive=FALSE)
gtexFiles <- list.files(path="rse-gtex/rse_gtex", pattern="rse*", full.names=TRUE, recursive=FALSE)
allFiles <- c(sraFiles, gtexFiles, "rse-tcga/rseTCGA.rda", "rse-tcga/rse_ESCA_TCGA.rda")

load("tissuePredictions/samplesWithPrediction_16_09_22.RData")
selectedSamples <- rownames(samplesWithPrediction)
str(selectedSamples)


#file = allFiles[10]

perChunkExp <- sapply(allFiles, function(file){
  
  loadedObject <- load(file)
  
  sreObjects <- get(loadedObject[1])
  
  #sometimes single RSE is not in list. Put in list of one to make code uniform
  if(!is.list(sreObjects)){
    sreObjects <- list(sreObjects)
  }
  
  #sreObject <- sreObjects[[1]]
  
  perStudyExp <- lapply(sreObjects, function(sreObject){
    studyExp <- sreObject@assays@data@listData$raw_counts
    return(studyExp[,colnames(studyExp) %in% selectedSamples, drop = F])
  })
  
  return(do.call(cbind, perStudyExp))
  
})

str(sreObject)

selectedSamplesExp <- do.call(cbind, perChunkExp)
str(selectedSamplesExp)
all(selectedSamples %in% colnames(selectedSamplesExp ))
table(selectedSamples %in% colnames(selectedSamplesExp ))



#Some samples are duplicated in the chunks, now make sure only one is in the matrix
uniqueSamplesIndex <- match(selectedSamples, colnames(selectedSamplesExp))
selectedSamplesExp <- selectedSamplesExp[,uniqueSamplesIndex]



#save(selectedSamplesExp, file = "perTissueNormalization/selectedSamplesRawExpression.RData")



