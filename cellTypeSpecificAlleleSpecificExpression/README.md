#Cell type specific allele specific expression

##Introduction
This module can be used for the detection of allele specific expression, 
with the novel option of including cell type specific effects when dealing 
with heterogenous tissues.


### Internal workings

The cell type specific (CTS) allele specific expression (ASE) module is separated into three sub-modules:

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

This sub module combines multiple ASfiles (Produced previously in the ASreads submodule)
of samples, and determines the probability of allelic imbalance per **SNP**.
This module uses multiple methods of finding the allelic imbalance using a likelihood ratio test (LRT):

- Binomial LRT
- Beta binomial LRT

When dealing with heterogeneous tissue, and cell proportions are available:

- CTS Binomial LRT
- CTS beta binomial LRT

Per SNP, a P value is determined. The P value is the probability of finding **no** ASE.
Lower P values mean more evidence for ASE at this SNP.
Usually the Beta Binomial (BB) test is better at correction for sample specific biasses,
thus, the Results of the Beta Binomial should be given more weight.  


#####2. ASEperRegion

This sub module combines multiple ASfiles (Produced previously in the ASreads submodule)
of samples, and determines the probability of allelic imbalance per genomic **region**.
This module uses multiple methods of finding the allelic imbalance using a likelihood ratio test (LRT):

- Binomial LRT
- Beta binomial LRT

