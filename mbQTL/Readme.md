# MbQTL - QTL analysis with beta-distributed null
This repository contains a version of MetaQTL (https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline)
that uses a different multiple testing strategy than the original MetaQTL that was used previously used in the eQTLgen and sc-eQTLgen studies.
The main differences between MetaQTL and this version are that this version no longer makes use of the TriTyper genotype format, but rather uses
tabix indexed VCF files as input, and that this version implements a different multiple testing correction strategy. 
While the previous version of MetaQTL determined a null distribution over the top associated variants of all genes, the method
implemented here is more akin to FastQTL (http://fastqtl.sourceforge.net/) or QTLTools (https://qtltools.github.io/qtltools/).
The main difference between FastQTL and this QTL software, is that while FastQTL assesses associations over all included samples,
this software allows on the fly meta-analysis of different datasets. Additionally, while FastQTL generally uses Pearson correlations to
determine relationships between phenotype and genotype, this software uses Spearman correlations by default.

# Compilation
This software depends on the 'genetica-libaries' repository (https://github.com/molgenis/systemsgenetics/tree/master/genetica-libraries) for various functions.
Download the full https://github.com/molgenis/systemsgenetics repository, and import as a maven project in e.g. IntelliJ Idea. Then, build the genetica-libraries module first,
and then this module.

# Download
A binary version of this software will be made available soon.

#Usage
A usage description will be added soon
