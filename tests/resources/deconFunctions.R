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

decon_eQTL_Test <- function(expVector, genotypeVector, counts, ...){
  cellTypeNames <- colnames(counts)
  tData <- data.frame (Exp= expVector, Gt= genotypeVector, counts)
  fullModel <- lm(as.formula(deconFormula(cellTypeNames)),data=tData)
  print(fullModel$coefficients)
  exit()
  pVals <- sapply(1:ncol(counts), function(x){
    ctModel <- lm(as.formula(deconFormula(cellTypeNames, minus=x)),data=tData)
    tmp <- anova(fullModel,ctModel)
    z <- tmp[2,6]
    #print(tmp)
    return(z)})
  names(pVals) <- cellTypeNames
  return(pVals) 
}

setwd("/Users/NPK/UMCG/projects/deconvolution/expData")
counts = read.delim("counts.txt",sep="\t", check.names=FALSE, header=TRUE,row.names=1) 
eqtl_dosage = read.delim("dsgTable_testing.txt",sep="\t", check.names=FALSE, header=TRUE) 
colnames(eqtl_dosage) <- gsub(".", "-", colnames(eqtl_dosage), fixed=TRUE)# changin "." for "-"
expTable_corr_addmean = read.delim("expTable_Corrected_addMean_snpname.txt",sep="\t", check.names=FALSE, header=TRUE) # check.names to avoid changin "-" for "."
rownames(expTable_corr_addmean) <- expTable_corr_addmean[,1]
expTable_corr_addmean[,1] <- NULL
df <- data.frame() 
for (i in 1:nrow(expTable_corr_addmean)){
  decon_eQTL_Test(unlist(expTable_corr_addmean[i,]), unlist(eqtl_dosage[i,]),counts)
  #write(decon_eQTL_Test(unlist(expTable_corr_addmean[i,]), unlist(eqtl_dosage[i,]),counts),append="TRUE",file="deconResults_expTable_corr_addmean_full_redone.txt")
}

shapiro <- function(x) {
  result = shapiro.test(as.numeric(x))
  pval = result$pvalue
}

readExp <- function(){
  expTable= read.delim("/Users/NPK/UMCG/projects/deconvolution/expData/expTable_Corrected_addMean_snpname.txt",sep="\t", check.names=FALSE, header=TRUE) # check.names to avoid changin "-" for "."
  rownames(expTable) <- expTable[,1]
expTable[,1] <- NULL
}




expTableMean = apply(expTable, 1, mean)
expTableStd = apply(expTable, 1, sd)
expTableNoMean = expTable - expTableMean
expTableNoMeanNoStd = expTableNoMean/expTableStd
expTableAddExponential = expTableNoMean**2
vector = c()
for (i in 1:nrow(expTable)){
 result = shapiro.test(as.numeric(expTableAddExponential[i,]))
 #vector <- c(vector, result$p.value)
 if(result$p.value > 0.05){
   print(result$p.value )
   hist(as.numeric(expTableAddExponential[i,]), title=result$p.value)
 }
}

