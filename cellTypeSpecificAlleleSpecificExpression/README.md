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

This sub-module combines multiple ASfiles (Produced previously in the ASreads submodule)
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
thus, the results of the Beta Binomial should usually be considered to be more accurate. 


#####3. ASEperRegion

This sub-module combines multiple ASfiles (Produced previously in the ASreads submodule)
of samples, and determines the probability of allelic imbalance per genomic **region**.
This module uses multiple methods of finding the allelic imbalance using a likelihood ratio test (LRT):

- Binomial LRT
- Beta binomial LRT

Per Region, and test SNP, a P-value is determined, in the same was as the ASEperSNP module.
Currently this module is functional, but not yet sufficiently tested to ensure  
a user obtains correct results. 


##Start: Binomial and Beta Binomial test

This section will describe how to do the most simple analysis: non-CTS ASEperSNP, 
in later sections, more options will be added to the analysis described here, 
to make the analysis better. 


**Minimum of data required:**

- A file with aligned sequence reads for multiple samples, and an index.
- A genotype file containing at least the genotype information of the bam file.
- A coupling file that indicates which bam files belong to which genotypes.


#####STEP 0: Preprocessing of data

**preprocessing your sequence alignment files**

This file will be used for the ASreads sub-module. It requires files to be in the 
bam format, a binary form of the [sam](http://genome.sph.umich.edu/wiki/SAM) format. You may be able to convert your file 
into the bam format using the program [samtools](https://github.com/samtools/samtools).


When `Suzie-RNAseq1.bam` is your bam file, you can provide it with an index by using samtools.
you can then index the bam using:
```
samtools index Suzie-RNAseq1.bam
```
A file named: `Suzie-RNAseq1.bam.idx` will be add in the working directory, 
this contains the index for the specific bam.
Now the bam file can be used for reading by the ASreads sub-module.
Please index all the bam files that you want to use for analysis.


**Preprocessing genotype files**

The sub-module ASreads accepts two types of file formats:

1. TriTyper
2. VCF

The TriTyper format is considerably faster to read than the VCF, therefore we will continue with the TriTyper format.
Conversion of the genotype format into TriTyper can be done using the [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer)
This guide refers to their wiki for the conversion of your genotype format into TriTyper format.

**Creating a coupling file**
The coupling file is used to match the sample data (bam file name) to the individual genotype data (name in the genotype file).
A coupling file is a tab delimited file, with one individual - sample pair per line. 
Individuals are in the first column (same name as in the genotype file), samples are in the second column (same name as the bam filenames, without the .bam extension).


In the most simple case you would have one individual in the genotype file: _Suzie_, with the subsequent sample file _Suzie-RNAseq1.bam_
```
Suzie\tSuzie-RNAseq1
```

Adding an individual is as simple as adding a similar line to the file, for instance _Peter_ with the sample Peter-RNAseq1.bam
Making the file the following:
```
Suzie\tSuzie-RNAseq1
Peter\tPeter-RNAseq1
```

For purposes of examples the guide will keep using the names Suzie and Peter for clarity.

#####STEP 1: ASREADS

After preprocessing the data into your 



 






 