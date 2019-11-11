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
#' @import edgeR
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' \code{dCell.predict} Runs dCell.predict
#'
#' @param dCell.exp A matrix or data.frame of normalized (by gene lenght and library size (RPKM, TMM, TPM...)) and scaled expression data, samples by columns and genes by row. Preferably comming from \code{dCell.expProcessing}.
#' @param dCell.models A list containing the DeconCell models.
#' @param res.type Whether the prediction per DeconCell models should be sumarized by either \code{median} or \code{mean}.
#' @return A matrix of predicted cell proportions for all available cell subpopoulation in DeconCell.
#' @return A vector containing the performance evaluation of all available models.
#' @export
#' @examples
#' data(dCell.models)
#' x <- dCell.predict(dCell.exp, dCell.models=dCell.models)
#' x$dCell.prediction
#' x$Evaluation
dCell.predict <- function(dCell.exp, dCell.models, res.type="median", custom=FALSE){

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

  if(custom){
    colnames(dCell.prediction) <- names(dCell.eval) <- names(dCell.models)
  }

  return(list(dCell.prediction=dCell.prediction, Evaluation= dCell.eval))
}


#' Run Decon-cell to generate new gene signature and characterize marker genes.
#' \code{dCell.run} Runs dCell.run to generate gene expression signatures to predict cell proportions of bulk tissues
#'
#' @param exp A matrix or data.frame of normalized (by gene lenght and library size (RPKM, TMM, TPM...)) and scaled expression data, samples by columns and genes by row. Preferably comming from \code{dCell.expProcessing}.
#' @param proportions A matrix or data.frame of cell proportions within each bulk RNA-Seq sample. Samples by row and cell types by column.
#' @param iterations Number of iterations (and models) which would be calculated
#' @param test.size Proportion of samples to be considered as test, 100 - \code{test.size} would be considered de proportion used for training.
#' @return a dCell.models-like object, a list (with a \code{length == ncol(proportion.mat))}, each element of the list contains the number of models defined by \code{iterations}
#' @export
#' @examples
#' dCell.run()
dCell.run <- function(exp,
                      proportions,
                      iterations= 10){

  n.cellTypes <- ncol(proportions)
  cat("INFO\t", "A total of ", n.cellTypes, "cell types are considered", "\n")

  ##Proportions and gene expression sample indexing
  if(all(colnames(exp) %in% rownames(proportions))  == FALSE ){
    stop("ERROR\t Same number of samples are needed in both expression and proportions matrix", "\n")
  } else {
    sample.index <- intersect(colnames(exp), rownames(proportions))
    exp <- exp[,sample.index]
    proportions <- proportions[sample.index,]
    cat("INFO\tSamples indexed", "\n")
  }


  ##Setting number of samples for training and tests
  deconCell.models.per.CT <- list()
  deconCell.marker.genes <- list()

  for(i.cellType in 1:n.cellTypes){

    i.ct.proportions <- proportions[,i.cellType]
    names(i.ct.proportions) <- rownames(proportions)

    i.dCell.run <- dCell.run.single(exp = exp, iterations = 5,
                     i.ct.proportions=i.ct.proportions)

    deconCell.models.per.CT[[i.cellType]] <- i.dCell.run$decon.cell.models
    deconCell.marker.genes[[i.cellType]] <- i.dCell.run$marker.genes.summary
  }
  names(deconCell.models.per.CT) <- colnames(proportions)
  names(deconCell.marker.genes) <- colnames(proportions)


  return(list(deconCell.models.per.CT= deconCell.models.per.CT,
              deconCell.marker.genes = deconCell.marker.genes))
}


#' Run DeconCell to generate new gene signature for a single cell type.
#' \code{dCell.run.single} Runs dCell.run to generate gene expression signatures to predict cell proportions of bulk tissues
#'
#' @param exp A matrix or data.frame of normalized (by gene lenght and library size (RPKM, TMM, TPM...)) and scaled expression data, samples by columns and genes by row. Preferably comming from \code{dCell.expProcessing}.
#' @param i.ct.proportions A matrix or data.frame of cell proportions within each bulk RNA-Seq sample. Samples by row and cell types by column.
#' @export
dCell.run.single <- function(exp, i.ct.proportions,
                             n.samples.train,
                             n.samples.test= NULL,
                             iterations= 10){

  decon.cell.models.list <- list()
  marker.gene.list <- list()
  cat(paste0("\nINFO\t Starting ",iterations," eNet iterations \n"))
  j <- 1
  while(j <= iterations) {
    cat(paste0("\nINFO\t Starting eNet iteration",j, "\n"))
    ## Define samples used of iteratrion
    eNetRes <- glmnet.wrapper(y=i.ct.proportions,
                              x= t(data.matrix(exp)),
                              parallel=F,
                              return.type="custom",
                              alpha.runs = c(0.01, 0.03 ,0.05, 0,075, 0.1, 0.15, 0.2))

    if(!length(eNetRes$coef) >= nrow(exp)) {
      marker.gene.list[[j]] <- names(eNetRes$coef)[-1]
      decon.cell.models.list[[j]] <- eNetRes$coef
      #cat("[INFO]\tDone iteration ", j, "\n")
      j <- j+1
    }
  }
  #return list of marker genes (genes (features) present in 70% of models)
  marker.genes.summary <- get.marker.genes(marker.gene.list = marker.gene.list, iterations = iterations)

  return(list(marker.genes.summary= marker.genes.summary,
              decon.cell.models=decon.cell.models.list))
}

#' Wrapper function for glmnet
#' \code{get.marker.genes}
#'
#' @param coef.list a list of character vectors
#' @export
#' @return A character vector with all marker genes for a dCell.run
#'
#'
get.marker.genes <- function(marker.gene.list, iterations, frequency.percentage=0.7){
  freqs <- table(unlist(marker.gene.list))
  marker.genes <- names(freqs)[which(freqs >= ceiling(iterations * frequency.percentage))]
  return(marker.genes)
}

#' Wrapper function for glmnet
#' \code{glmnet.wrapper}
#'
#' @param y numeric vector, in this case proportions which will be predicted
#' @param x numeric matrix, in this case gene expression data from bulk tissues
#' @param newx if not NULL, numeric matrix, in this case gene expression data from bulk tissues for which the proportions will be predicted for
#' @param newy if not NULL, numeric vector,
#' @param alpha.runs numeric vector of the alphas that will be testes
#' @param parallel boolean, if TRUE then parallel
#' @param return.type type of return: "custom", "model" or "all"
#' @param method if newy != NULL, which correlation method will be used to calculate correlation coeficient between predicted proportions and measured ones.
#' @param use usage of NAs within y, x, newx and newy
#' @param remove.outliers boolan, if TRUE all outliers (as defined in boxplots) will be removed from y.
#' @import glmnet
#' @return a summarized object of the glmnet function
#' @export
#' @examples
#' glmnet.wrapper()
glmnet.wrapper <- function(y, x, newx=NULL, newy=NULL, alpha.runs=NULL,
                           parallel=F, return.type="custom",
                           method="spearman", use="pairwise.complete.obs",
                           remove.outliers=FALSE) {
  y             <- y[!is.na(y)]
  if (remove.outliers) {
    # Remove outliers in the response
    # (improves prediction but reduces sample size)
    y <- y[!y %in% boxplot.stats(y)$out]
  }
  ol            <- intersect(rownames(x), names(y))
  y             <- y[ol]
  x             <- as.matrix(x[ol, , drop=F])
  models.alpha  <- list()
  alpha.vec     <- c()

  if(is.null(alpha.runs)){
    alpha.runs    <- c(0.0, 0.01, 0.025, 0.05, 0,075, 0.1, 0.2, 0.3, 0.4, 0.5 ,0.75, 0.85, 0.9, 1)
  }

  # for each alpha
  for (a in 1:length(alpha.runs)) {
    model             <- cv.glmnet(y=y, x=x, alpha=alpha.runs[a], parallel=parallel)
    #alpha.vec[a]      <- model$cvm[model$lambda == model$lambda.min]
    # Use 1se to prevent overfitting
    alpha.vec[a]      <- model$cvm[model$lambda == model$lambda.1se]
    models.alpha[[a]] <- model
  }

  # Arrange the results in a custom list format which contains predictions
  if (return.type=="custom") {
    result                <- list()
    best.alpha            <- which(alpha.vec == min(alpha.vec))
    best.model            <- models.alpha[[best.alpha]]
    best.lambda           <- best.model$lambda.1se
    coef                  <- as.matrix(coef(best.model, s=best.lambda))
    coef                  <- coef[coef != 0,]

    if (!is.null(newx)) {
      pred                <- predict(best.model$glmnet.fit, newx=as.matrix(newx),
                                     s=best.lambda, type="response")
      pred.names          <- rownames(pred)
      pred                <- as.numeric(pred)
      names(pred)         <- pred.names
      result$pred         <- pred

      if(!is.null(newy)){
        pred.cor            <- cor(pred, newy, method=method, use=use)

        result$pred.cor     <- pred.cor
        result$n.test       <- length(newy)
        result$mse          <- sum((pred - newy)^2)/length(pred)
      }

    }

    result$best.alpha      <- alpha.runs[best.alpha]
    result$best.lambda     <- best.lambda
    result$best.model      <- best.model
    result$coef            <- coef
    result$n.train         <- length(y)

  } else if (return.type == "model") {
    result <- models.alpha[[which(alpha.vec == min(alpha.vec))]]
  } else if (return.type == "all") {
    result <- models.alpha
  }

  return(result)
}


