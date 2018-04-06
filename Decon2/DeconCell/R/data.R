########
#### Documentation of data in the deconCell package
########


#' DeconCell gene list
#' Vector of gene names (ENSEMBL identifiers) used in all of the Decon cell models
#'
#' @format Lists of 100 elemented with coefficients, names are coded in ENSEMBL identifiers
"dCell.geneList"

#' DeconCell models
#' Lists fo 33 lists containing the coefficients used in the prediction of cell proportions using expression values
#'
#' @format Lists of 100 elemented with coefficients, names are coded in ENSEMBL identifiers
"dCell.models"

#' Cell subpopulation names
#' Code and names for the predicted cell subpopulations in the DeconCell package
#'
#' @format data.frame of 33 rows and 2 columns, rownames are the codes foe the dCell.models lists
"dCell.names"

#' Expression data
#' Example RNA-seq data from 5 healthy individuals
#'
#' @format A matrix with 42992 quantified genes and 5 samples
"count.table"

#' Cell Proportion data
#' Example cell proportions of circulating immune cell subpopulations from 5 healthy individuals
#'
#' @format A matrix with 33 quantified cell proportions 5 samples
"cell.proportions"
