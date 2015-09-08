#Cell type specific allele specific expression

##Introduction
This module can be used for the detection of allele specific expression, 
with the novel option of including cell type specific effects when dealing 
with heterogenous tissues.


##Summary of internal workings

The cell type specific allele specific expression module is separated into three sub-modules:

1. ASreads
2. ASEperSNP
3. ASEperRegion

###1. ASreads

This sub module is able to quantify the number sequencing reads that overlap an
allele at heterozygous SNP sites, by integrating any quantitative sequencing 
experiment (like RNA-seq or proseq) with with genotype information.

The data created in this module can be used in the later modules for the 
quantification of allele specific expression.    


###2. ASEperSNP

This module combines 