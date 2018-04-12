Authors: Raúl Aguirre-Gamboa and Niek de Klein

# Introduction 
Decon2 is a statistical framework for estimating cell counts using molecular profiling such as expression or methylation data from heterogeneous samples (Decon-cell) and consecutive 
deconvolution of expression quantitative trait loci (Decon-eQTL) into each cell subpopulation.  
See DeconCell/ and Decon-eQTL/ for their documentation.

# Installation

For DeconCell, in R do

```
library('devtools')
# NOTE: This can take a few minutes:
install_github("systemsgenetics/Decon2/DeconCell")
```

For Decon-eQTL, download ```Decon-eQTL-v*.*.*-jar-with-dependencies.jar``` from [here: Decon-eQTL/lastBuild](http://molgenis50.gcc.rug.nl/jenkins/job/systemsgenetics/Decon-eQTL$Decon-eQTL/lastBuild/).


# Example run

For this example run we will use faked expression and genotype data. We will go through both steps of deconvolution:  
  
1. Predict cellcount proportions  
2. Deconvolute cell type specific effects  
  
If you have measured cell types, you can skip to [Decon-eQTL](#decon-eqtl).

## Simulate count data

We will use the samples from the DeconCell vignette to simulate expression data. In R, do:

```r
library(matrixStats)
library(DeconCell)
library(truncnorm)

data("count.table")

# take the row means and row standard deviation of the vignette samples
rMeans <- rowMeans(count.table)
rSD <- rowSds(as.matrix(count.table))
names(rSD) <- names(rMeans)
means_and_SD <- data.frame(mean=rMeans, SD=rSD)

# number of samples to simulate
number_of_samples <- 100

# make sample names
sampleNames <- c()
for(i in 1:number_of_samples) {
  sampleNames <- c(sampleNames, paste0('sample_', i))
}

# For every gene sample from same mean and SD as vignette samples
# with expression >= 0
new_count_table <- data.frame(t(apply(
	means_and_SD, 1, function(x) rtruncnorm(number_of_samples, a=0,
	                                        mean=x[['mean']],
	                                        sd=x[['SD']]))
	))
colnames(new_count_table) <- sampleNames

# set NA's to 0
new_count_table[is.na(new_count_table)] <- 0
# Write the simulated data to a file
write.table(new_count_table, 
			  'example_data/count.table.simulated.txt', 
			  sep='\t',quote=F,col.names=NA)
```

## DeconCell

In R

```r
library(DeconCell)
# read the count data
count.table <- read.table("example_data/count.table.simulated.txt")
# load the model
data("dCell.models")

# Normalize gene expression (see DeconCell vignette for more details)
dCell.exp <- dCell.expProcessing(count.table, trim = TRUE)

# predict the cell counts
prediction <- dCell.predict(dCell.exp, dCell.models, res.type = "median")
```

The prediction object contains the predicted cellcounts and model evaluation. For more info on model information, see the DeconCell Vignette. For this example run, we will use only 6 main cell types: ```Granulocytes```, ```Monocytes```, ```CD4+```, ```CD8+```, ```NK```, and ```B cells```. For Decon-eQTL the cellcounts have to be scaled to sum to 100.

In R

```r
# select relevant cell types
predicted.cellcounts <- prediction$dCell.prediction[,c('Granulocytes',
                                                       'B cells (CD19+)',
                                                       'CD4+ T cells',
                                                       'CD8+ T cells',
                                                       'NK cells (CD3- CD56+)',
                                                       'Monocytes (CD14+)')]
# scale to sum to 100
predicted.cellcounts.scaled <- (predicted.cellcounts/rowSums(predicted.cellcounts))*100

# write the scaled cellcounts to outfile
write.table(predicted.cellcounts.scaled, 'example_data/predicted.cellcounts.scaled.txt', 
			sep='\t', quote=F, col.names=NA)
```

You now how scaled, predicted cell counts that can be used as input for Decon-eQTL.

## Decon-eQTL

With the predicted cell counts you can now run Deco-eQTL. In your terminal, run:

```bash
java -jar Decon-eQTL-v1.0.1-jar-with-dependencies.jar \
	--cellcount example_data/predicted.cellcounts.scaled.txt \
	--expression example_data/count.table.simulated.txt \
	--genotype example_data/fake_genotype_dosages.txt \
	--snpsToTest example_data/gene_snp_file.txt \
	--outfolder example_output/
```

A file called ```deconvolutionResults.csv``` will be written to ```example_output/```. Below is one line from the deconvolution result file. The first 6 columns are the p-values for each of the celltypes, the middle 6 columns are the beta (effect size) of the cell type - genotype interaction effect, and the last 6 columns are the direction of the effect relative to the allele that is coded as 2 in the genotype file.


|	|Granulocytes_pvalue	|B cells (CD19+)_pvalue	|CD4+ T cells_pvalue	|CD8+ T cells_pvalue	|NK cells (CD3- CD56+)_pvalue	|Monocytes (CD14+)_pvalue	|Beta1_Granulocytes	|Beta2_B cells (CD19+)	|Beta3_CD4+ T cells	|Beta4_CD8+ T cells	|Beta5_NK cells (CD3- CD56+)	|Beta6_Monocytes (CD14+)	|Beta7_Granulocytes:GT	|Beta8_B cells (CD19+):GT	|Beta9_CD4+ T cells:GT	|Beta10_CD8+ T cells:GT	|Beta11_NK cells (CD3- CD56+):GT	|Beta12_Monocytes (CD14+):GT	|effectDirectionDosage2_Granulocytes	|effectDirectionDosage2_B cells (CD19+)	|effectDirectionDosage2\_CD4+ T cells	|effectDirectionDosage2_CD8+ T cells	|effectDirectionDosage2_NK cells (CD3- CD56+)	|effectDirectionDosage2_Monocytes (CD14+)
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---
ENSG00000103855\_snp_244	|0.8754191650343855	|1.0	|0.9855839296645753	|1.0	|0.8801056746139749	|1.0	|0.0	|0.0	|0.0	|0.0	|0.19545183472374217	|0.0	|0.008312507740549976	|0.0	|0.0014092158690805837	|0.0	|0.05818753733261271	|0.0	|-	|-	|-	|-	|+	|-

