
Authors: Raúl Aguirre-Gamboa and Niek de Klein

## Introduction 
Decon-eQTL_simple is as simple y ~ snp + snp : gt interaction  model to use as comparison to our full model. We suggest using 
[(Decon-eQTL)](https://github.com/molgenis/systemsgenetics/tree/master/Decon2/Decon-eQTL) for permorming deconvolution.


## Requirements
`Java version 1.8 or higher` (Tested using `1.8.0_121`)  
Decon-eQTL-v*.*.*.jar


### For compiling
commons-cli v1.3.1  
commons-csv v1.2  
commons-io v2.4  
commons-lang3 v3.4  
commons-math3 v3.6  
jsci-core  
jsci-mathmlimpl  
jsci-sci

## minimal usage example
    
    java -jar deconvolution -c <file containing cellcounts> \
                            -e <file containing expression data> \
                            -g <file containing genotypes>  \
                            -o <output directory> \
                            -sn <file with SNP and gene combination to test>
    
For a full guide see [Decon2 manual](https://github.com/molgenis/systemsgenetics/tree/master/Decon2)

## Options overview

    -c,--cellcount <file>                     Cellcount file name
    -e,--expression <file>                    Expression file name
    -g,--genotype <file>                      Genotype file name
    -help                                     print this message
    -o,--outfolder <path>                     Path to folder to write output to
    -oe,--outputPredictedExpression           Write output file with predicted expression
    -of,--outfile <file>                      Outfile name of deconvolution results (will be written in outfolder)
    -sn,--snpsToTest <file>                   Tab delimited file with first column gene name, second column SNP name. Need to match with names from genotype and expression files.
