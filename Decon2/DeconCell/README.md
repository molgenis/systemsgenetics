# NOTE: Development has moved to https://github.com/molgenis/systemsgenetics/tree/master/Decon2

---
title: "DeconCell"
author: "Raúl Aguirre-Gamboa and Niek de Klein"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction 
DeconCell is an r package containing models for predicting the proportions of circulating immune cell subpopulations using bulk gene expression data from  whole blood. Models were built using an elastic net and training in 95 healthy dutch volunteers from the [500FG cohort](http://www.humanfunctionalgenomics.org/site/?page_id=82) with FACS quantification of 73 circulating cell subpopulations as described in our previous [publication](http://www.cell.com/cell-reports/fulltext/S2211-1247(16)31473-5). 
For additional details on methods and results please go our [manuscript](link to be updated). 


## Install the package from github

```{r}
library(devtools)
install_github("raguirreg/DeconCell")
```


## Pre-processing example data
Let's load and pre-process our example data. These are 5 samples with > ~40k genes quantified. These are gene read counts, we need to approximate the example data to a normal-like distribution and account for library sizes. In order to do this, we use the `dCell.expProcessing` function. This function will perform a TMM normalization (as described in the [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html)package) a log2(counts+1) and scale (z-transformation) per gene.

```{r}
library(DeconCell)
library(edgeR)

data("count.table")
dCell.exp <- dCell.expProcessing(count.table, trim = TRUE)
```

## Prediction of cell propotions

```{r}
data("dCell.models")
prediction <- dCell.predict(dCell.exp, dCell.models, res.type = "median")
head(prediction$dCell.prediction)
head(prediction$Evaluation)
```

## Correlation coeficient between of predicted and measured values
```{r}
data("cell.proportions")
library(reshape2)
library(ggplot2)
data("dCell.names")
pData <- data.frame(PearsonCor= diag(cor(cell.proportions, prediction$dCell.prediction)), 
                    CTs = dCell.names[colnames(cell.proportions), "finalName"], 
                    Subpop = dCell.names[colnames(cell.proportions), "broadSubpopulations"])
ggplot(pData, aes(y=PearsonCor , x= CTs, fill=Subpop))+
  geom_bar(stat="identity", alpha=0.8)+
  geom_hline(yintercept = 0.5, alpha=0.5, color="red")+
  coord_flip()+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()
  
```

