---
title: "DeconCell"
author: "Raúl Aguirre-Gamboa and Niek de Klein"
date: "`r Sys.Date()`"
output: rmarkdown::github_document
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
install_github("molgenis/systemsgenetics/Decon2/DeconCell")
```

##Vignette 
For a quick start and tutorial please visit our  [vignette](http://htmlpreview.github.io/?https://github.com/molgenis/systemsgenetics/blob/master/Decon2/DeconCell/inst/doc/my-vignette.html)



