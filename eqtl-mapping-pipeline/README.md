#EQTL mapping pipeline manual


This software allows for (e)QTL mapping using linear models and direct meta-analysis of such data.

##Downloading the software
You can download the latest version of the software here: [Latest version](http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$eqtl-mapping-pipeline/lastBuild/).
Make sure to download the stand alone jar: `eqtl-mapping-pipeline-******-jar-with-dependencies.jar`

Please note that the manual refers to eQTLMappingPipeline.jar, while the name of the package described above may be different (because of different version numbers etc)

##Before you start


###Java Virtual Machine
Our software is written in Java, which makes the software both fast and portable across multiple operating systems. Executing a java program is similar to executing a normal program or app, although there are some considerations:

* eQTL mapping heavily relies on available memory. Make sure your machine is 64-bit and has lots of memory installed (at least 4Gb). 

* Please make sure your version of java is up-to-date. Use at least the 64-bit version of Java version 6 (also called version 1.6). If you are running on Windows, you can download the so-called Java Runtime Environment (JRE) from: http://www.java.com. If you are running Linux, the virtual machine may be present in the proprietary section of your package manager, or may be available via http://www.java.com.

* Java executables are called jar-files. The eQTL mapping pipeline is such a jar-file. You can execute it by using the following command from a terminal / console:
```    
    java –jar eQTLMappingPipeline.jar 
```

* You need to specify the amount of memory available to the program using the command-line switch –Xmx. The amount of memory can be specified in megabytes (using an m suffix) or in gigabytes (using a g suffix). To be sure your computer is running java at 64-bit, please add the switch –d64. Both switches (–Xmx and –d64) should be called prior to the –jar switch. An example where the program is allowed to use 4gb of memory:
```	
    java –d64 –Xmx4g –jar eQTLMappingPipeline.jar 
```
* Try to increase the –Xmx amount when you get Out-Of-Memory-Errors or errors involving ‘heap space’

**IMPORTANT NOTE: In this manual, we assume you understand the principle that you need to allocate sufficient amounts of RAM and therefore we excluded the –Xmx switch from the example commands. Please be aware that you should use it, as most of the commands require a substantial amount of memory!**

**The eQTL mapping pipeline is a command line program, which makes the user interface not very intuitive. In order to help you a bit, an overview of available switch options is displayed when a command is incomplete or incorrect. Furthermore, each mode of the program also has its own overview of available switches with a small description of its functionality. For example: “java –jar eQTLMappingPipeline.jar” produces a list of available modes, while “java –jar eQTLMappingPipeline.jar  --mode metaqtl” produces a list of all available options for metaqtl. You can also find all possible switches in the manual section “command line options”.**

#Step by step eQTL analysis
This is a step by step guide which will guide you through the QTL mapping process using our software. Please note that this step by step guide illustrates only a part of the capabilities of our software. Our general method consists of six steps described below:

<ol>
<li>Step 1 - Preparation phenotype data</li>
<li>Step 2 - Preparation genotype data</li>
<li>Step 3 - Phenotype data normalization</li>
<li>Step 4 - Mapping mix-ups</li>
<li>Step 5 - Determining the optimum number of PCs to remove</li>
<li>Step 6 - Perform the final QTL analysis</li>
</ol>

##Step 1 - Preparation phenotype data
Because our software uses a nonparametric test by default, you can use virtually any continuous data as trait values to map a variety of QTL effects. However, currently the normalization tools provided with this package are focused on array based methylation data, array based (Illumina) expression data, and preprocessed RNA-seq data (e.g. transcript level quantified data) or (GC)-RMA processed Affymetrix data. 
Please format your phenotype data in a simple tab separated text file. The format of this file is described here: [Data Formats - Phenotype data](https://github.com/harmjanwestra/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file).

**Throughout the manual, we will refer to your phenotype file as `traitfile`**

###Illumina array based expression data

Because Illumina has several ways to annotate and preprocess their gene expression arrays, we provide a detailed instruction below, that will yield raw, untouched gene expression levels.

Use the GenomeStudio files of your expression arrays. **Please do not use any normalization, transformation, imputation, or background correction method that is offered by Illumina’s GenomeStudio! Just load your raw data, without doing any normalization, transformation, or background correction.**

-	Export a so-called ‘Final Report’ using GenomeStudio, including array_address_id (as the probe identifier) and the average expression signals per sample for all probes (click 'Analysis' in the top menu bar, then 'Reports').
-	If you create more than one Final Report, merge the Final Reports, so that all data will be combined in one Expression Matrix File
-	Rows should contain the different probes, and columns should contain the different sample IDs
-	Remove any header information that Genome studio might produce: the header of the matrix is expected at the first row
-	Create a new header (expected at the first row and column) for both rows and columns, describing the sample names (columns) and the probe IDs (rows). Note that you use the probe array address for your platform as the probe name.
-   The final file should look like this format: [Data Formats - Phenotype data](https://github.com/harmjanwestra/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file).


##Step 2 - Preparation of genotype data
Our software is able to use both unimputed called genotypes, as well as imputed genotypes and their dosage values. However, currently our software can only interpret files that are in the [TriTyper](http://genenetwork.nl/wordpress/trityper/) format. We are working on a generic method of genotype input through [Genotype IO](https://github.com/harmjanwestra/systemsgenetics/tree/master/Genotype-IO). In the mean time, users of the eQTL mapping pipeline should convert their data to TriTyper using ImputationTool. This tool is also integrated in the eQTL mapping pipeline, and can be called using the following command (no separate download required):

````
java –d64 –Xmx4g –jar eQTLMappingPipeline.jar --imputationtool
````

The documentation for the ImputationTool can be found at the [ImputationTool repository page](https://github.com/molgenis/systemsgenetics/tree/master/imputation-tool).

###Check your data
After converting your genotype data to TriTyper format, a number of files should have been created in your output directory. For a description of the different files, please refer to [Data formats - TriTyper](link)

##Step 3 - Phenotype data normalization
Generally, continuous trait data needs to be normalized prior to applying statistical testing. Our software package provides different ways of normalizing your data contained in the `traitfile`. 

Our general normalization strategy for Illumina based array data consists of the following steps:
1. Quantile normalization
2. Log<sub>2</sub> transformation
3. Probe centering and scaling (Z-transform): (Expression<sub>Probe,Sample</sub> = Expression<sub>Probe,Sample</sub> – MeanProbe) / Std.Dev.<sub>Probe</sub>
4. (Optionally) Removal of covariates - Correct gene expression data for first 4 PCs of the GWAS data – to remove possible population stratification
    - *We run a Generalized Linear Model with the probes as dependent variables, and the GWAS PCs as orthogonal covariates. For the remainder of the analysis, we use the residuals of this model.*
5. Principal component adjustment
    - *In the Principal Component Analysis, the software tries to convert the normalized and standardized expression results into a set of values of uncorrelated variables, called Principal Components (PCs). The number of PCs is equal to the number of samples, and the first PCs explain a higher variance than the last PCs. By adjusting the normalized and standardized expression data for a set of PCs (like we did with the removal of the covariates – use of a Generalized Linear Model), we try to remove batch effects in the data. By removing them in an incremental way, we try to find the optimal number of PCs to remove.*

After probe centering, sample z-transformation, and removal of covariates, a tab-separated gzipped plaintext file is created with an identical number of rows and columns as the input file. This means these files can subsequently be used during eQTL mapping.

Apart from the methods described above, our software can also perform Mtransformation for methylation beta values obtained from Illumina methylation arrays.

###Preparations
Note down the full path to your `traitfile`. The output of the normalization tool will be written to the same directory by default. Optionally, note down the full path to the file containing covariates. We will refer to the covariate file as `covariatefile`. The format of this file is identical to the `traitfile`: [File formats - Phenotype data](link). 

###Commands to be issued
To run the general normalization strategy described above, you can run the following command:
```
java –d64 –Xmx4g –jar eQTLMappingPipeline.jar --mode normalize --in traitfile
```

You can specify an output directory with the following command (specifying an `outdir`):
```
java –d64 –Xmx4g –jar eQTLMappingPipeline.jar --mode normalize --in traitfile --out outdir
```

To run the general normalization strategy described above, and correct for covariates:
```
java –d64 –Xmx4g –jar eQTLMappingPipeline.jar --mode normalize --in traitfile --adjustcovariates --cov covariatefile
```

Individual elements of the normalization strategy can also be separately executed. For example to only run Quantile normalization, and Log<sub>2</sub> transformation run the following command:
```
java –d64 –Xmx4g –jar eQTLMappingPipeline.jar --mode normalize --in traitfile --qqnorm --logtransform
```

**Note:**
Several other parameters can be set to customize your normalization strategy (e.g. number of PCs to remove, step size for PC removal, handling of missing values, etc). However, the order of procedures is fixed (Quantile Normalize > Log<sub>2</sub> transform > covariate adjustment > centering and scaling > PCA adjustment), irregardless of the order of each command line switch. To review the available options for normalization, issue the following command:
```
java –d64 –Xmx4g –jar eQTLMappingPipeline.jar --mode normalize
```

###Check your data
Running the general normalization procedure yields a number of files in the directory of your `traitfile`, or in the `outdir` if you specified one. 

|File|Description|
|----|-----------|
|**ExpressionData.txt.QuantileNormalized.txt.gz**|Quantile Normalized Expression Data|
|**ExpressionData.txt.QuantileNormalized.Log2Transformed.txt.gz**|Quantile Normalized Expression Data which is also Log<sub>2</sub> Transformed.|
|**ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.txt.gz**|Quantile Normalized Expression Data which is also Log2 Transformed. Probes were centered.|
|**ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz**|Quantile Normalized Expression Data which is also Log2 Transformed, probes were centered and samples were Z-transformed.|
|**ExpressionData.txt.PCAOverSamplesEigenvalues.txt.gz**|Eigenvalues created during eigenvalue decomposition of the gene expression sample correlation matrix (created from the Quantile Normalized, Log2 Transformed, Z-transformed data)|
|**ExpressionData.txt.PCAOverSamplesEigenvectors.txt.gz**|Eigenvectors created during eigenvalue decomposition of the gene expression sample correlation matrix (created from the Quantile Normalized, Log2 Transformed, Z-transformed data)|
|**ExpressionData.txt.PCAOverSamplesEigenvectorsTransposed.txt.gz**|Eigenvectors transposed|
|**ExpressionData.txt.PCAOverSamplesPrincipalComponents.txt.gz**|Principal Components describing the sample correlation matrix (created from the Quantile Normalized, Log2 Transformed, Z-transformed data)|
|**ExpressionData.txt.nPCAsOverSamplesRemoved.txt.gz**|Expression data, Quantile Normalized, Log2 Transformed, Z-transformed, with n Principal Components regressed out.|

##Step 4 - Mapping mix-ups
We have shown in a paper published in Bioinformatics (Westra et al.: [MixupMapper: correcting sample mix-ups in genome-wide datasets increases power to detect small genetic effects](http://bioinformatics.oxfordjournals.org/content/27/15/2104)), that sample mix-ups often occur in genetical genomics datasets (i.e. datasets with both genotype and gene expression data). Therefore, we developed a method called MixupMapper, which is implemented in the eQTL Mapping Pipeline. This program performs the following steps:
1.    At first a cis-eQTL analysis is conducted on the dataset:
    -	using a 250 kb window between the SNP and the mid-probe position
    -	performing 10 permutations to control the false discovery rate (FDR) at 0.05
2.	Calculate how well a gene expression array matches a genotype array. 
For details how this exactly works, please have a look at the paper or read the attachment. 

###Preparations

##Step 5 - Determining the optimum number of PCs to remove
##Step 6 - Perform the final QTL analysis

#File formats
This section lists the different file formats used by the software package.
##Probe annotation file
A probe annotation file is required when running a *cis*-eQTL analysis. This file describes where each probe/trait/gene is located on the genome. This file is a simple tab-separated text file, with a line for each probe, and a header.

**Please note: for reliable meta-analysis, probe annotations should be identical between datasets. Therefore, we recommend use of the array address for your platform as a probe identifier if you are using Illumina array based expression data.**

###File example
<pre>
Platform    HT12v4-ArrayAddress Symbol	Chr	ChrStart	ChrEnd	Probe     Seq
HT12v4      00001               GeneX	1       1504        1554        0       CGCTCCCCTTATAACTT-etc.
HT12v4      00002               GeneY	11      19900       19950       1       GGATCCCAGATTCCCT-etc.
HT12v4      00003               GeneZ	23      101         151         2       TTCTCCAGAGTCGAGC-etc.
</pre>


##Phenotype file, covariate file
Phenotype and covariate files have the same basic format. We use a tab separated table, with individuals on columns and probes or covariates on the rows. 

###File example
<pre>
ProbeArrayAddress    Sample1	Sample2		Sample3
00001		    	1.640		1.553		1.441
00002		    	1.671		0.201		5.321
00003		    	1.126		1.710		2.569
00004		    	1.129		1.002		1.313
</pre>

**Please note: if you are exporting your data from R, note that the first column often does not have a column identifier. Our software expects n+1 columns in both the header as well as the data (where n equals the number of individuals in your dataset).**

##TriTyper genotype data
The TriTyper format consists of several files, each describing an aspect of the genotype data. 

|File name | Required | Description |
|----------|----------------|
| **GenotypeMatrix.dat** | X | Binary file containing genotype data. GenotypeMatrix.dat has the following file size: (number of SNPs * 2) * number of individuals. The ImputedDosageMatrix.dat should be half this size.|
| **ImputedDosageMatrix.dat** | - | Binary file containing imputed genotype dosage values. The ImputedDosageMatrix.dat should be half this size in bytes compared to the GenotypeMatrix.dat.|
| **SNPs.txt** | X | The list of SNPs that are encoded within the GenotypeMatrix.dat file. One line per SNP.|
| **SNPMappings.txt** | X | The list of SNPs that are encoded within the GenotypeMatrix.dat file. One line per SNP, tab-separated: first column contains the chromosome number, second column contains the SNP position, and third column contains the SNPID (rs ID).|
| **Individuals.txt** | X | The list of individuals that are encoded within the GenotypeMatrix.dat file. One line per individual. **Do not change the order of the individuals in this file, or the number of individuals in this file. You can change the individual identifiers, although duplicates are not allowed.**|
| **PhenotypeInformation.txt** | X | This file describes the phenotypes of the individuals. One line per individual, 4 columns per individual: individual ID, case/control status, include/exclude a certain individual, gender (female/male). **This file does not have to contain all individuals contained in Individuals.txt and can be used to exclude certain individuals from the analysis**|

Please note that you can update the PhenotypeInformation.txt file. For a population based approach, you should designate all participants “control”. “Include” or “Exclude” determines whether you include or exclude a participant into the analysis. Finally, you need to add gender information. For individuals of unknown gender, you can use a random string, as long as this string does not match 'female' or 'male'. Individuals that are in the Individuals.txt file, but are not present in the PhenotypeInformation.txt file, will be excluded from analysis.

The SNP mappings will depend upon the genome build of the reference dataset used during imputation. To reliably perform a meta-analysis, SNP mappings should be identical across datasets. If you want to update the SNP mappings (for example to a newer build), you can either rename or remove the original SNPMappings file. 

###File examples

SNPs.txt: example of 3 SNPs
<pre>
rs11511647
rs12218882
rs10904045
</pre>

SNPMappings.txt: example of 3 SNPs
<pre>
10	62765	rs11511647	
10	84172	rs12218882	
10	84426	rs10904045	
</pre>


SNPMappings.txt: example of 3 SNPs
<pre>
10	62765	rs11511647
10	84172	rs12218882
10	84426	rs10904045
</pre>

Individuals.txt: example of 3 samples
<pre>
Sample1
Sample2
Sample3
</pre>

PhenotypeInformation.txt: example of 3 samples
<pre>
Sample1		control	include	female
Sample2		control	include 	male
Sample3		control	include	female
</pre>