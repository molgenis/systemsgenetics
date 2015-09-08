#Cell type specific allele specific expression

##Introduction
This module can be used for the detection of allele specific expression, 
with the novel option of including cell type specific effects when dealing 
with heterogenous tissues.


### Internal workings

The cell type specific allele specific expression module is separated into three sub-modules:

1. ASreads
2. ASEperSNP
3. ASEperRegion

#####1. ASreads

This sub module is able to quantify the number sequencing reads that overlap an
allele at heterozygous SNP sites, by integrating any quantitative sequencing 
experiment (like RNA-seq or proseq) with genotype information.

The data created in this module can be used in the later modules for the 
quantification of allele specific expression.    


#####2. ASEperSNP

This module combines multiple ASfiles (Produced previously in the ASreads submodule)
of samples, and determines the probability of allelic imbalance per SNP.
This module uses multiple methods of finding the allelic imbalance using a likelihood ratio test (LRT):

1. Binomial LRT
2. Beta binomial LRT
When dealing with heterogeneous tissue, and cell proportions are available:
3. Cell type specific Binomial LRT
4. Cell type specicic beta binomial LRT