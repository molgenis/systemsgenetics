###--------
### r.a.aguirre.gamboa
### DeconCell
### Functions

#' Pre-process RNASeq based expression data.
#' \code{dCell.expProcessing} Performes a TMM normalization, log2(counts+1) and scaling.
#'
#' @param count.table A table of quantified gene expression, samples by columns and genes by row.
#' @param trim boolean, whather the expression data should be trimmed to include only genes present in DeconCell models.
#' @return A normalized and scaled expression matrix.
#' @export (do export this in NAMESPACE)
#' @examples
#' data(coun.table)
#' dCell.exp <- dCell.expProcessing(count.table)
#' dim(dCell.exp)
#' head(dCell.exp)
dCell.expProcessing <- function(count.table, trim=TRUE){
  D<-edgeR::DGEList(counts=count.table)
  d <- calcNormFactors(D)
  scalar <- d$samples$lib.size*d$samples$norm.factors/exp(mean(log(d$samples$lib.size*d$samples$norm.factors)))
  scal.mat <- outer(rep(1,nrow(d$counts)), scalar)
  normExp <- log2((d$counts/scal.mat)+1)
  sdExp <- apply(normExp, 1, sd)
  zeros <- which(sdExp==0)
  if(length(zeros) != 0){
    standarizedExp <- (normExp[-zeros,]-apply(normExp[-zeros,], 1, mean)) / sdExp[-zeros]
  } else{
    standarizedExp <- (normExp-apply(normExp, 1, mean)) / sdExp
  }
  if(trim){
    standarizedExp <- dCell.exp.trimming(standarizedExp)
  }
  return(standarizedExp)
}

#' Trims a dCell.exp matrix for only dCell genes
#' \code{dCell.predict} Loads a the gene list used in all dCell models and trims a dCell.exp matrix.
#'
#' @param dCell.exp A matrix of preprocessed gene expression data
#' @return Percentage of coeficient represented in the new data
#' @export (do export this in NAMESPACE)
#' @examples
#' dCell.exp.trimming(dCell.exp)
dCell.exp.trimming <- function(dCell.exp){
  if(!exists("dCell.geneList")){
    data(dCell.geneList)
  }
  gene.index <- intersect(dCell.geneList, rownames(dCell.exp))
  dCell.exp <- dCell.exp[gene.index,]
  propFound <- (length(gene.index)/length(dCell.geneList)) *100
  cat("[INFO]\t Total of", round(propFound,digits = 2), "% genes from dCell are found")
  return(dCell.exp)
}


#' Evaluate percentage of coeficient used in a model
#' \code{dCell.predict} Predicts cell proportions of 32 cell subpopulations usin decon Cell models.
#'
#' @param genes.ensembl Character vector of ensemble gene ids overlapping between the decon cell model and the new gene expression.
#' @param ct.model A single decon cell model. A numeric vector with ensembl IDs as names.
#' @return Percentage of coeficient represented in the new data
#' @export (do export this in NAMESPACE)
#' @examples
#' dCell.evaluate(genes.ensembl, ct.model)
dCell.evaluate <- function(genes.ensembl, ct.model){
  if(names(ct.model)[1] == "(Intercept)") {
    ct.model <- ct.model[-1]
  }
  return((sum(abs(ct.model[genes.ensembl]))/sum(abs(ct.model)))*100)
}


#' Predict a single cell subpopoulation using a DeconCell model
#' \code{dCell.predict.single} Predicts cell proportions of one out of the 32 cell subpopulations using the DeconCell models.
#'
#' @param dCell.exp A table of normalized and scaled expression, samples by columns and genes by row. Preferably comming from \code{dCell.expProcessing}.
#' @param dCell.models A list containing the DeconCell models.
#' @param res.type Whether the prediction per DeconCell models should be sumarized by either \code{median} or \code{mean}.
#' @return A vector of predicted cell proportions for a given set of models.
#' @return Mean of percentage of coeficients covered given
#' @export (do export this in NAMESPACE)
#' @examples
#' data(dCell.models)
#' x <- dCell.predict.single(dCell.exp, dCell.model=dCell.models[[1]])
#' x$Predicted
#' x$Prediction.eval
dCell.predict.single <- function(dCell.exp, dCell.model, res.type="median"){
  ## predicted matrix, samples by row
  pred.matrix <- matrix(NA,
                       nrow=ncol(dCell.exp),
                       ncol=length(dCell.model))
  pred.eval <- rep(NA , times=length(dCell.model))
  #for loop for possible inclussion of more samples or models
  for(i in 1:length(dCell.model)){
      i.model <- dCell.model[[i]]
      gene.index <- intersect(rownames(dCell.exp), names(i.model[-1]))
      # remove genes from the model that are not present in the dCell.exp
      pred.eval[i] <- dCell.evaluate(genes.ensembl = gene.index, ct.model = i.model)
      # remove genes from the model that are not present in the dCell.exp
      i.model <- c(i.model[1], i.model[gene.index])
      pred.matrix[,i]  <- as.numeric(i.model) %*% (as.matrix(rbind(1, dCell.exp[names(i.model[-1]),])))
  }

  if(res.type=="median"){
    pred.ct <- apply(pred.matrix, 1, median)
  } else if(res.type=="mean"){
    pred.ct <- apply(pred.matrix, 1, mean)
  }
  names(pred.ct) <- colnames(dCell.exp)
  return(list(Predicted= pred.ct, Prediction.eval= mean(pred.eval)))
}



#' Predict all available cell subpopoulation in DeconCell
#' \code{dCell.predict.single} Runs dCell.predict
#'
#' @param dCell.exp A table of normalized and scaled expression, samples by columns and genes by row. Preferably comming from \code{dCell.expProcessing}.
#' @param dCell.models A list containing the DeconCell models.
#' @param res.type Whether the prediction per DeconCell models should be sumarized by either \code{median} or \code{mean}.
#' @return A matrix of predicted cell proportions for all available cell subpopoulation in DeconCell.
#' @return A vector containing the performance evaluation of all available models.
#' @export (do export this in NAMESPACE)
#' @examples
#' data(dCell.models)
#' x <- dCell.predict(dCell.exp, dCell.models=dCell.models)
#' x$dCell.prediction
#' x$Evaluation
dCell.predict <- function(dCell.exp, dCell.models, res.type="median"){

  #check for dCell.names
  if(!exists("dCell.names")){
    data(dCell.names)
  }
  if(!exists("dCell.models")){
    data(dCell.models)
  }
  dCell.prediction <- matrix(NA,
                            ncol=length(dCell.models),
                            nrow=ncol(dCell.exp))
  dCell.eval <- rep(NA,
                    times=length(dCell.models))

  for(i in 1:length(dCell.models)){
    i.dCell <- dCell.predict.single(dCell.exp=dCell.exp , dCell.model = dCell.models[[i]], res.type =res.type )
    dCell.prediction[,i] <- i.dCell$Predicted
    dCell.eval[i] <- i.dCell$Prediction.eval
  }
  rownames(dCell.prediction) <- colnames(dCell.exp)
  colnames(dCell.prediction) <- names(dCell.eval) <- dCell.names[names(dCell.models),"finalName"]
  return(list(dCell.prediction=dCell.prediction, Evaluation= dCell.eval))
}



