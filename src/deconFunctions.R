deconFormula <- function(countNames, minus=NULL){
  if(is.null(minus)){
    formulaConstruct <- paste("Exp~",  
                              paste0(countNames,"+", collapse=""),
                              paste0(countNames, ":Gt", "+" ,collapse=""), 
                              collaps="")
    formulaConstruct <- paste(substr(formulaConstruct , 1, nchar(formulaConstruct)-1), 
                              "-1", 
                              collapse="")
  } else {
    formulaConstruct <- paste("Exp~",  
                              paste0(countNames,"+", collapse=""),
                              paste0(countNames[-minus], ":Gt", "+" ,collapse=""),
                              collapse="")
    formulaConstruct <- paste(substr(formulaConstruct , 1, nchar(formulaConstruct)-1), 
                              "-1" , collapse="")

  } 
  return(formulaConstruct)
}

decon_eQTL_Test <- function(expVector, genotypeVector, counts){
  print('starting deconvolution')
  cellTypeNames <- colnames(counts)
  tData <- data.frame (Exp= expVector, Gt= genotypeVector, counts)
  print('start formula')
  formula <- as.formula(deconFormula(cellTypeNames))
  print('start full model')
  print(head(tData))
  fullModel <- lm(formula,data=tData)
  print('done full model')
  pVals <- sapply(1:ncol(counts), function(x){
              ctModel <- lm(as.formula(deconFormula(cellTypeNames, minus=x)),data=tData)
              z <- anova(fullModel,ctModel)[2,6]
              return(z)})
  names(pVals) <- cellTypeNames
  print('returning')
  return(pVals) 
}
#print(cellTypes)
expressionVector <- unlist(expressionVector)
names(expressionVector) <- geneNames
genotypeVector <- unlist(genotypeVector)
names(genotypeVector) <- geneNames
cellcountDf <-  as.data.frame(matrix(unlist(cellcountTable), nrow=length(unlist(cellcountTable[1]))))
rownames(cellcountDf) <- cellcountDf[,1]
cellcountDf[,1] <- NULL
cellcountDf <- cellcountDf[-1,]
colnames(cellcountDf) <- cellTypes
print(decon_eQTL_Test(expressionVector,genotypeVector,cellcountDf))
#write.table(expVector, file = "expVector.txt",
#      append = FALSE, sep = "\t", quote = FALSE,
#      col.names=FALSE)
#write.table(genotypeVector, file = "genotypeVector.txt",
#            append = FALSE, sep = "\t", quote = FALSE,
#            col.names=FALSE)
#write.table(counts1, file = "counts1.txt",
#            append = FALSE, sep = "\t", quote = FALSE,
#            col.names=NA)
