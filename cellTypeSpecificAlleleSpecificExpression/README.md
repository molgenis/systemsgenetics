#Cell type specific allele specific expression

##Introduction
This module can be used for the detection of allele specific expression, 
with the novel option of including cell type specific effects when dealing 
with heterogenous tissues.

**Requirements:**

This program requires java version 8, and was tested on ubuntu 14.04 LTS.
other Linux- and OS X based systems should also work.


### Internal organisation

The cell type specific (CTS) allele specific expression (ASE) module is separated into three sub-modules:

1. ASreads
2. ASEperSNP
3. ASEperRegion

#####1. ASreads

This sub-module is able to quantify the number sequencing reads that overlap an allele at heterozygous SNP sites, by integrating any quantitative sequencing experiment (like RNA-seq or proseq) with genotype information.

The data created in this module can be used in the later modules for the 
quantification of allele specific expression.    


#####2. ASEperSNP

This sub-module combines multiple ASfiles (Produced previously in the ASreads sub-module) of samples, and determines the probability of allelic imbalance per **SNP**.
This module uses multiple methods of finding the allelic imbalance using a likelihood ratio test (LRT):

- Binomial LRT
- Beta binomial LRT

When dealing with heterogeneous tissue, and cell proportions are available:

- CTS Binomial LRT
- CTS beta binomial LRT

Per SNP, a P value is determined. The P value is the probability of finding **no** ASE.
Lower P values mean more evidence for ASE at this SNP.
Usually the Beta Binomial (BB) test is better at correction for sample specific biases, thus, the results of the Beta Binomial should usually be considered to be more accurate. 


#####3. ASEperRegion

This sub-module combines multiple ASfiles (Produced previously in the ASreads submodule) of samples, and determines the probability of allelic imbalance per genomic **region**.
This module uses multiple methods of finding the allelic imbalance using a likelihood ratio test (LRT):

- Binomial LRT
- Beta binomial LRT

Per Region, and test SNP, a P-value is determined, in the same was as the ASEperSNP module.

**WARNING**

Currently this module is functional, but not yet sufficiently tested to ensure  
a user obtains correct results. 


## Basic operation: Binomial and Beta Binomial test

This section will describe how to do the most simple analysis: non-CTS ASEperSNP, 
in later sections, more options will be added to the analysis described here, to 
be able to add extra functionality


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
you can then index the bam using the following command in your terminal:
```
samtools index Suzie-RNAseq1.bam
```
A file named: `Suzie-RNAseq1.bam.bai` will be add in the working directory, 
this contains the index for the specific bam.
Now the bam file can be used for reading by the ASreads sub-module.
Please index all the bam files that you want to use for analysis.


**Preprocessing genotype files**

The sub-module ASreads accepts two types of file formats:

1. TriTyper
2. VCF

The TriTyper format is considerably faster to read than the VCF at this moment, therefore this guide will continue with the TriTyper format.
Conversion of the genotype format into TriTyper can be done using the [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer). 
This guide refers to their wiki for the conversion of your genotype format into TriTyper format.
For example purposes, the directory containing TriTyper information will be set to the following: `Suzie-Peter_Genotype/`

```
bgzip -c Suzie-RNAseq1.vcf > Suzie-RNAseq1.vcf.gz
tabix -p vcf Suzie-RNAseq1.vcf.gz 
sh GenotypeHarmoinzer.sh --input Suzie-RNAseq1.vcf.gz  --outputType TRITYPER --output Suzie-Peter_Genotype/
```


**Creating a coupling file**

The coupling file is used to match the sample data (bam file name) to the individual genotype data (name in the genotype file).
A coupling file is a tab delimited file, with one individual - sample pair per line. 
Individuals are in the first column (same name as in the genotype file), samples are in the second column (same name as the bam filenames, without the .bam extension).


In the most simple case you would have one individual in the genotype file: _Suzie_, with the subsequent sample file _Suzie-RNAseq1.bam_
```
Suzie    Suzie-RNAseq1
```

Adding an individual is as simple as adding a similar line to the file, for instance _Peter_ with the sample Peter-RNAseq1.bam
Making the file the following:
```
Suzie    Suzie-RNAseq1
Peter    Peter-RNAseq1
```
Please note that the spacing between individual and sample needs to be a tab character (\t), not spaces.
This file can be saved as: `suzie-peter_coupling.txt`
For example purposes, this guide will keep using the names Suzie and Peter for clarity.

#####STEP 1: ASreads

After preprocessing, ASread files can be created.
ASread files provide information per Bi-allelic SNP (other variants are not supported) on where a read maps.
Formatting of ASfiles is in tab delimited fashion with the following data per column:

```
1. Chromosome name
2. Position on the chromosome
3. Variant name (rs name)
4. Reference base (based on the genotype file)
5. Alternative base (based on the genotype file)
6. Number of reads mapping to the reference base
7. Number of reads mapping to the alternative base
8. Number of reads mapping to neither the reference not alternative base.
9. Genotype of the individual in a "[C, G]" formatting.s 
```

To create this file for individual _Suzie_, 

This is the first time we're going to use the program, please download or build it, currently 
it is named as the following but this will probably change (TODO):
`cellTypeSpecificAlleleSpecificExpression-1.0.2-SNAPSHOT-jar-with-dependencies.jar`

To create an AS file for suzie named: `Suzie-ASreads.txt`,
 we can run the following command in bash:
```
java -jar cellTypeSpecificAlleleSpecificExpression-1.0.2-SNAPSHOT-jar-with-dependencies.jar \
    --action ASreads \
    --output Suzie-ASreads.txt \
    --coupling_file suzie-peter_coupling.txt \
    --genotype_location Suzie-Peter_Genotype/ \
    --bam_file Suzie-RNAseq1.bam\
    > Suzie_ASreads_output.txt
```

When this run has finished (takes quite some time, depending on the number of SNPs.), 
We need to do the other individuals as well (Peter in the example case), 
by changing the file which to output and the bam file for input.

Now all individuals in this example have an AS file, we can start the ASE testing.

**Check your output**
Based on the AS files, you can determine if your genotypes and coupling are correct.
When you open an ASfile and find high numbers that map to no allele and or hardly 
any balance in the alleles with heterozygotes, you could have some problems with 
incorrect coupling, or incorrect genotyping.


#####STEP 2: Testing for ASE at single SNPs.

Before testing for ASE, one needs to specify the parameters for each test:
These are:

-Which ASfiles do you want to include
-How many heterozygotes does a variant need to have before starting to test. default:1
-How many reads must overlap a variant before starting to test. default:10

Which ASfiles you need to include must be stored in a file with a file per line as such(in case of testing Suzie and Peter):
```
Suzie_ASreads_output.txt
Peter_ASreads_output.txt
```
We will save this as `Suzie-Peter-ASFiles.txt`

**Warning**

Make sure you the ordering of the AS files is the same in all files specified in the above file. 
Otherwise, the testing will not be done.

**/Warning**

Now head to your terminal and run the following command:

```
java -jar cellTypeSpecificAlleleSpecificExpression-1.0.2-SNAPSHOT-jar-with-dependencies.jar \
    --action ASEperSNP \
    --output PeterSuzieTests
    --as_locations Suzie-Peter-ASFiles.txt \
    --minimum_hets 1 \
    --minimum_reads 10\
    > PeterSuzieBinomialBetaBinomialOut.txt
```

This command will run for some time depending on the number of SNPs and the number of individuals used for testing.

The output in PeterSuzieBinomialBetaBinomialOut.txt will show, in a human readable 
form, what tests are performed, their input and their output.

The Specific testing output will be saved in a number of files with the test 
type appended to the original output parameter, namely:
`PeterSuzieTests_dispersionFile.txt, PeterSuzieTests_BinomialResults.txt and PeterSuzieTests_BetaBinomialResults.txt`
These files show the resulting statistics for specific tests.

When running Cell type specific tests (not shown in this example), 
the output file will be saved as <output>_CTSBinomialResults.txt and 
<output>_CTSBetaBinomialResults.txt, for the binomial and beta binomial cell type specific tests respectively.

**Check your output**
To get a feel for the data I suggest taking a look at what is piped through standard out (saved in PeterSuzieBinomialBetaBinomialOut.txt in this example case)
and you may be able to sort the Results files based on P-value or chi squared value.

Currently, the overall structure of results files is not yet the same everywhere, 
but the first 6 columns are the same throughout the ASEperSNP test:

```
1. Chromosome 
2. Position on chromosome
3. SNP name (rs number)
4. Number of samples heterozygote for this SNP
5. P value
6. Chi squared value
```

In our example case, if you would like to see the 5 SNPs that are most likely 
ASE SNPs in the beta binomial tests, you could sort based on the chi squared value, like 
this in the terminal:

```
sort -n -k 6,6 PeterSuzieTests_BetaBinomialResults.txt | tail -n 5
```

This concludes the basic usage section, preprocessing, isolating ASreads and 
testing for ASE using a binomial and beta binomial were discussed.

**Warning**

Please note that having two individuals in your analysis is usually not enough to provide you with robust estimates of ASE.

**\Warning**


## Additional features (1): Cell type specific tests:

When your sequenced sample is a mixture of tissues (for instance when sequencing 
blood samples), you may be interested in how allele specific expression is a 
factor in a specific tissue. For this purpose, a cell type specific testing 
feature is available in this module that  will automatically produce results 
when a phenotype file is available.


A phenotype file is a file with some ratio per individual:
```
0.25
0.43
```
Where the ordering of the file is the same ordering as the ordering of individuals in the `--as_locations` option.
When integrating this phenotype file with the Suzie Peter example (see Step 2 in the basic usage), this will mean 
that for some cell type, the proportion of the celltype in the tissue is 0.25  (25%) in Suzie and 0.43 (43%) in Peter

Saving this file as `Suzie-Peter-Phenotype.txt`, the command run will look like this:

```
java -jar cellTypeSpecificAlleleSpecificExpression-1.0.2-SNAPSHOT-jar-with-dependencies.jar \
    --action ASEperSNP \
    --output PeterSuzieTests
    --as_locations Suzie-Peter-ASFiles.txt \
    --minimum_hets 1 \
    --minimum_reads 10\
    --pheno_file Suzie-Peter-Phenotype.txt\
    > PeterSuzie_CTS_Out.txt
```

Results of the Cell type specific effects can be found in the results from the standard out and from the output.


**Warning**

Only 2 samples is very little information for the detection of CTS
effects, even more so than looking at non CTS ASE for itself.

**\Warning**


##Additional features(2): ASE per Region

**Warning**

This module is still under development and only sparsely tested, only non CTS effects are currently supported.

**\Warning**

###### Basic mechanism

The ASEperRegion sub-module can be used to determine the ASE effects per region to 
increase statistical power by combining multiple SNPs in the same transcript or 
_gene region_ and comparing them to some _test SNP_, as shown in the following scheme below:

```
## Fictitious example of ASE per region:


           -----------------------Some test region--------------------------

            Test SNP        |///////////gene region////////////|
               X
                                       region SNP #
                                1        2        3        4

Suzie:        Het              Het      Het      Het      Hom
Allele 1:   ~~~T~~~~~~~~~~~~|___A________T________G________C___|~~~~~~~~~~~~

Allele 2:   ~~~C~~~~~~~~~~~~|___G________A________A________C___|~~~~~~~~~~~~
                                
                                
Peter:        Het              Het      Hom      Hom      Hom
Allele 1:   ~~~C~~~~~~~~~~~~|___A________A________G________A___|~~~~~~~~~~~~

Allele 2:   ~~~T~~~~~~~~~~~~|___G________A________G________A___|~~~~~~~~~~~~


Walt:         Hom              Het      Hom      Het      Het
Allele 1:   ~~~C~~~~~~~~~~~~|___A________A________A________C___|~~~~~~~~~~~~

Allele 2:   ~~~C~~~~~~~~~~~~|___G________A________G________A___|~~~~~~~~~~~~


-------------------------------Legend:--------------------------------------

Het: Hetererozygote                     Hom:        Homozygote                         
"X": Position of the test SNP           [A,G,T,C]:  Base at SNP position
"~": Non SNP in  region                 "_":        Non SNP in gene region
```

In this example, only individuals with heterozygous _test SNPs_ are used in
the analysis. 
Which means that only Suzie and Peter will have data that will be used, data from Walt will not be considered. 
In this example the test SNP is a C/T polymorphism.
All reads that overlap the _gene region_ on the allele _test SNP_ with genotype C are counted, and compared to
all the reads from the _gene region_ on the allele with _test SNP_ with genotype T.

For indicative purposes, reads overlapping _region SNP_ 1, 2 and 3 can be used 
in Suzie, while _region SNP_ 1 is the only one that can be used in Peter, due to it being the only heterozygote in Peters _gene region_.
It is very important to make sure the region is well defined as this can be a confounding factor, 
where some gene may be a different gene or splice variant.   


##### Required data:

Required data for this sub-module is as follows (not yet final)


- A file containing gene regions and test regions (at time of writing, the same)
- Locations of AS files
- Genotype information with phasing information per SNP. 
- Coupling file for the phased genotypes and AS files.

The genotype information should be phased information of at least some of the 
SNPs in the test region.

The coupling file is the same as in the Basic usage example.

The AS reads file is also the same as in the Basic usage example.

The regions file will specify the test regions.

Currently only binomial and beta binomial work.






