#EQTL mapping pipeline manual
This software allows for (e)QTL mapping using linear models and direct meta-analysis of such data through weighted Z-score analysis.

##Questions or suggestions?
You can contact the authors of this software at westra.harmjan@gmail.com, or lude@ludesign.nl for questions or suggestions regarding the software or the manual. This manual was written in cooperation with Marjolein Peters and Joyce van Meurs.

##Manual contents

1. Downloading the software
2. Before you start
3. Step by step QTL analysis
    * Step 1 - Preparation phenotype data
    * Step 2 - Preparation genotype data
    * File Checklist
    * Step 3 - Phenotype data normalization
    * Step 4 - MixupMapper
    * Step 5 - Determining the optimum number of PCs to remove
    * Step 6 - Perform the final QTL analysis
4. Additional analyses and advanced settings
    * Meta-analysis and settings file
    * Multiple linear regression using interaction model
    * Conditional analysis
5. File formats
    * Probe annotation file
    * Phenotype file, covariate file
    * Genotype - phenotype coupling
    * TriTyper genotype data
6. Frequently asked questions

##Downloading the software
You can download the latest version of the software here: [Latest version](http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$eqtl-mapping-pipeline/lastBuild/).

Make sure to download the stand alone jar: `eqtl-mapping-pipeline-******-jar-with-dependencies.jar`

Please note that the manual refers to eQTLMappingPipeline.jar, while the name of the package described above may be different (because of different version numbers etc).

##Before you start

###Path definitions and commands
Please note that our software expects full paths, although shorter paths wil also work in most cases. So, if you are on Windows, a full path to a genotype directory would be similar to `c:\path\to\genotype\dir\` and a full path to a file would be `c:\path\to\genotype\directory\file.txt`. Linux and Mac OS use different path separators. On these systems, these paths would be similar to the following `/path/to/genotype/dir/` and `/path/to/genotype/dir/file.txt`. Our main point here is that when pointing to a directory, use a 'trailing slash'

This manual will combine references to paths with commands that need to be issued for a certain task. For example, at some point in this manual we refer to your phenotype data as `traitfile`, which will be printed in a grey box. Commands will also be in grey boxes, and can make references to paths defined earlier (to keep the manual readable), as follows:

````
java -jar eQTLMappingPipeline.jar --mode metaqtl --inexp traitfile
````

To run the command in the example above, you have to replace the string `traitfile` with the full path of your `traitfile` after the command line switch `--inexp`. So if the full path to your `traitfile` would be `/path/to/traitfile.txt`, the final command would be:

````
java -jar eQTLMappingPipeline.jar --mode metaqtl --inexp /path/to/traitfile.txt
````

###Java Virtual Machine
Our software is written in Java, which makes the software both fast and portable across multiple operating systems. Executing a java program is similar to executing a normal program or app, although there are some considerations:

* QTL mapping heavily relies on available memory. Make sure your machine is 64-bit and has lots of memory installed (at least 4Gb). 

* Please make sure your version of java is up-to-date. Use at least the 64-bit version of Java version 6 (also called version 1.6). The software should also work on version 7 and higher. If you are running on Windows, you can download the so-called Java Runtime Environment (JRE) from: http://www.java.com. If you are running Linux, the virtual machine may be present in the proprietary section of your package manager, or may be available via http://www.java.com.

* Java executables are called jar-files. The eQTL mapping pipeline is such a jar-file. You can execute it by using the following command from a terminal / console:
```    
    java –jar eQTLMappingPipeline.jar 
```

* You need to specify the amount of memory available to the program using the command-line switch –Xmx (case-sensitive!). The amount of memory can be specified in megabytes (using an m suffix) or in gigabytes (using a g suffix). To be sure your computer is running java at 64-bit, please add the switch –d64. Both switches (–Xmx and –d64) should be called prior to the –jar switch. An example where the program is allowed to use 4gb of memory:
* 
```    
    java –d64 –Xmx4g –jar eQTLMappingPipeline.jar 
```

* Try to increase the –Xmx amount when you get Out-Of-Memory-Errors or errors involving ‘heap space’

**IMPORTANT NOTE:** In this manual, we assume you understand the principle that you need to allocate sufficient amounts of RAM and therefore we excluded the –Xmx switch from the example commands. Please be aware that you should use it, as some of the commands may require a substantial amount of memory (depending on your dataset and settings)!

###General information about software
* The eQTL mapping pipeline is a command line program, which makes the user interface not very intuitive. In order to help you a bit, an overview of available switch options is displayed when a command is incomplete or incorrect. Furthermore, each mode of the program also has its own overview of available switches with a small description of its functionality. For example: ```java –jar eQTLMappingPipeline.jar``` produces a list of available modes, while ```java –jar eQTLMappingPipeline.jar  --mode metaqtl``` produces a list of all available options for metaqtl.
* The software is able to process GZipped text files for most of the input files (not files in .tar archives however), which allows you to save some space on your hard drive.


#Step by step eQTL analysis
This is a step by step guide which will guide you through the QTL mapping process using our software. Please note that this step by step guide illustrates only a part of the capabilities of our software. However, this guide does explain the different command line switches. [Additional types of analyses can be found in this section.](link). Our general method consists of six steps described below:

* Step 1 - Preparation phenotype data
* Step 2 - Preparation genotype data
* File Checklist
* Step 3 - Phenotype data normalization
* Step 4 - MixupMapper
* Step 5 - Determining the optimum number of PCs to remove
* Step 6 - Perform the final QTL analysis

##Definitions
Througout the manual, references to different full paths will be made. Here is an overview of these paths:

* The phenotype file will be referred to as `traitfile`
* The probe/trait/gene annotation file will be referred to as `annotationfile`
* The full path of your genotype data will be referred to as `genotypedir`
* The file linking phenotype individuals to genotype individuals will be referred to as `genotypephenotypecoupling`
* The file containing covariates will be defined as `covariatefile`

##Step 1 - Preparation phenotype data
Because our software uses a nonparametric test by default, you can use virtually any continuous data as trait values to map a variety of QTL effects. However, currently the normalization tools provided with this package are focused on array based methylation data, array based (Illumina) expression data, and preprocessed RNA-seq data (e.g. transcript level quantified data) or (GC)-RMA processed Affymetrix data. 

Please format your phenotype data in a simple tab separated text file. The format of this file is described here: [Data Formats - Phenotype data](https://github.com/harmjanwestra/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file).

Some analyses will require an annotation of the probes/traits/genes in your `traitfile`. We store this annotation in a separate file. The format of this file is described here: [File Formats - Annotation File](https://github.com/harmjanwestra/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file). 

###Illumina array based expression data

Because Illumina has several ways to annotate and preprocess their gene expression arrays, we provide a detailed instruction below, that will yield raw, untouched gene expression levels.

Use the GenomeStudio files of your expression arrays. **Note:** Please do not use any normalization, transformation, imputation, or background correction method that is offered by Illumina’s GenomeStudio! Just load your raw data, without doing any normalization, transformation, or background correction.

* Export a so-called ‘Final Report’ using GenomeStudio, including array_address_id (as the probe identifier) and the average expression signals per sample for all probes (click 'Analysis' in the top menu bar, then 'Reports').
* If you create more than one Final Report, merge the Final Reports, so that all data will be combined in one Expression Matrix File
* Rows should contain the different probes, and columns should contain the different sample IDs
* Remove any header information that Genome studio might produce: the header of the matrix is expected at the first row
* Create a new header (expected at the first row and column) for both rows and columns, describing the sample names (columns) and the probe IDs (rows). Note that you use the probe array address for your platform as the probe name.
* The final file should look like this format: [Data Formats - Phenotype data](https://github.com/harmjanwestra/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file).


##Step 2 - Preparation of genotype data
Our software is able to use both unimputed called genotypes, as well as imputed genotypes and their dosage values. However, currently our software can only interpret files that are in the [TriTyper](http://genenetwork.nl/wordpress/trityper/) format. We are working on a generic method of genotype input through [Genotype IO](https://github.com/harmjanwestra/systemsgenetics/tree/master/Genotype-IO). In the mean time, users of the eQTL mapping pipeline should convert their data to TriTyper using ImputationTool. This tool is also integrated in the eQTL mapping pipeline, and can be called using the following command (no separate download required):

```
java –jar eQTLMappingPipeline.jar --imputationtool
```

The documentation for the ImputationTool can be found at the [ImputationTool repository page](https://github.com/molgenis/systemsgenetics/tree/master/imputation-tool).

###Check your data
After converting your genotype data to TriTyper format, a number of files should have been created in your output directory. For a description of the different files, please refer to [Data formats - TriTyper](link)

#File Checklist
Before you continue in this manual, this is a good time check whether your files are in the correct format and whether you have all the required files ready.

* Check whether all required files are in your `genotypedir`.
* Check whether the `traitfile` is properly normalized.
* Check whether you have an `annotationfile`. The format is described here [File formats - Probe annotation file](link). 
    - Some QTL data may not have annotation (for example if you want to do a GWAS on lipid levels). In such cases an `annotationfile` is not required for QTL mapping. Please see [Step 6 - Running the final QTL mapping](link) for more details. 
    - For the other steps, the annotation file is required, although you can fool the software by setting your genotype genomic locations to 1 (both chromosome and chromosome position) and doing the same for your phenotype data using the `annotationfile`.
* If the sample identifiers differ between genotype and phenotype data, you have to create a file that links these identifiers together. This format is described here: [File formats - Genotype Phenotype coupling](link). This file is optional, although we will mention this file in this manual as `genotypephenotypecoupling`

##Step 3 - Phenotype data normalization
Generally, continuous trait data needs to be normalized prior to applying statistical testing. Our software package provides different ways of normalizing your data contained in the `traitfile`. 

Our general normalization strategy for Illumina based array data consists of the following steps:

1. Quantile normalization
2. Log<sub>2</sub> transformation
3. Probe centering and scaling (Z-transform): (Expression<sub>Probe,Sample</sub> = Expression<sub>Probe,Sample</sub> – MeanProbe) / Std.Dev.<sub>Probe</sub>
4. (Optionally) Removal of covariates - Correct gene expression data for first 4 PCs of the GWAS data – to remove possible population stratification
    * We run a Generalized Linear Model with the probes as dependent variables, and the GWAS PCs as orthogonal covariates. For the remainder of the analysis, we use the residuals of this model.*
5. Principal component adjustment
    * In the Principal Component Analysis, the software tries to convert the normalized and standardized expression results into a set of values of uncorrelated variables, called Principal Components (PCs). The number of PCs is equal to the number of samples, and the first PCs explain a higher variance than the last PCs. By adjusting the normalized and standardized expression data for a set of PCs (like we did with the removal of the covariates – use of a Generalized Linear Model), we try to remove batch effects in the data. By removing them in an incremental way, we try to find the optimal number of PCs to remove.*

After probe centering, sample z-transformation, and removal of covariates, a tab-separated gzipped plaintext file is created with an identical number of rows and columns as the input file. This means these files can subsequently be used during eQTL mapping.

Apart from the methods described above, our software can also perform Mtransformation for methylation beta values obtained from Illumina methylation arrays.

###Preparations
Note down the full path to your `traitfile`. The output of the normalization tool will be written to the same directory by default. Optionally, note down the full path to the file containing covariates. We will refer to the covariate file as `covariatefile`. The format of this file is identical to the `traitfile`: [File formats - Phenotype data](link). 

###Commands to be issued
To run the general normalization strategy described above, you can run the following command:

```
java –jar eQTLMappingPipeline.jar --mode normalize --in traitfile
```

You can specify an output directory with the following command (specifying an `outdir`):

```
java –jar eQTLMappingPipeline.jar --mode normalize --in traitfile --out outdir
```

To run the general normalization strategy described above, and correct for covariates:

```
java –jar eQTLMappingPipeline.jar --mode normalize --in traitfile --adjustcovariates --cov covariatefile
```

Individual elements of the normalization strategy can also be separately executed. For example to only run Quantile normalization, and Log<sub>2</sub> transformation run the following command:

```
java –jar eQTLMappingPipeline.jar --mode normalize --in traitfile --qqnorm --logtransform
```

**Note:**
Several other parameters can be set to customize your normalization strategy (e.g. number of PCs to remove, step size for PC removal, handling of missing values, etc). However, the order of procedures is fixed (Quantile Normalize > Log<sub>2</sub> transform > covariate adjustment > centering and scaling > PCA adjustment), irregardless of the order of each command line switch. To review the available options for normalization, issue the following command:

```
java –jar eQTLMappingPipeline.jar --mode normalize
```

###Check your data
Running the general normalization procedure yields a number of files in the directory of your `traitfile`, or in the `outdir` if you specified one. 

|File|Description|
|----|-----------|
|**ExpressionData.txt.&shy;QuantileNormalized.&shy;txt.gz**|Quantile Normalized Expression Data|
|**ExpressionData.txt.&shy;QuantileNormalized.&shy;Log2Transformed.txt.gz**|Quantile Normalized Expression Data which is also Log<sub>2</sub> Transformed.|
|**ExpressionData.txt.&shy;QuantileNormalized.&shy;Log2Transformed.&shy;ProbesCentered.txt.gz**|Quantile Normalized Expression Data which is also Log2 Transformed. Probes were centered.|
|**ExpressionData.txt.&shy;QuantileNormalized.&shy;Log2Transformed.&shy;ProbesCentered.&shy;SamplesZTransformed.txt.gz**|Quantile Normalized Expression Data which is also Log2 Transformed, probes were centered and samples were Z-transformed.|
|**ExpressionData.txt.&shy;PCAOverSamplesEigenvalues.&shy;txt.gz**|Eigenvalues created during eigenvalue decomposition of the gene expression sample correlation matrix (created from the Quantile Normalized, Log2 Transformed, Z-transformed data)|
|**ExpressionData.txt.&shy;PCAOverSamplesEigenvectors.txt.gz**|Eigenvectors created during eigenvalue decomposition of the gene expression sample correlation matrix (created from the Quantile Normalized, Log2 Transformed, Z-transformed data)|
|**ExpressionData.txt.&shy;PCAOverSamplesEigenvectorsTransposed.txt.gz**|Eigenvectors transposed|
|**ExpressionData.txt.&shy;PCAOverSamplesPrincipalComponents.txt.gz**|Principal Components describing the sample correlation matrix (created from the Quantile Normalized, Log2 Transformed, Z-transformed data)|
|**ExpressionData.txt.&shy;nPCAsOverSamplesRemoved.txt.gz**|Expression data, Quantile Normalized, Log2 Transformed, Z-transformed, with n Principal Components regressed out.|

##Step 4 - MixupMapper
We have shown in a paper published in Bioinformatics ([Westra et al.: *MixupMapper: correcting sample mix-ups in genome-wide datasets increases power to detect small genetic effects*](http://bioinformatics.oxfordjournals.org/content/27/15/2104)), that sample mix-ups often occur in genetical genomics datasets (i.e. datasets with both genotype and gene expression data). Therefore, we developed a method called *MixupMapper*, which is implemented in the eQTL Mapping Pipeline. This program performs the following steps:
1. At first a *cis*-eQTL analysis is conducted on the dataset:
    -	using a 250 kb window between the SNP and the mid-probe position
    -	performing 10 permutations to control the false discovery rate (FDR) at 0.05
2. Calculate how well a gene expression array matches a genotype array. 
For details how this exactly works, please have a look at the paper or read the attachment. 

###Preparations/variables
1. Note down the full path to your TriTyper genotype data directory. We will refer to this directory `genotypedir`. 
2. Determine the full path of the trait data you want to use, and make sure this data is normalized (e.g. use Quantile Normalized, Log<sub>2</sub> transformed Illumina gene expression data). We will refer to this path as `traitfile`. 
3. Locate your phenotype annotation file: `annotationfile`. Also note down the platform identifier `platformidentifier`.
4. Locate your `genotypephenotypecoupling` if you have such a file. You can also use this file to test specific combinations of genotype and phenotype individuals. 
5. Find a location on your hard drive to store the output. We will refer to this directory as `outputdir`.

###Commands to be issued
The *MixupMapper* analysis can be run using the following command:

```
java –jar eQTLMappingPipeline.jar --mode mixupmapper --in genotypedir --out outdir --inexp traitfile --inexpplatform platformidentifier --inexpannot annotationfile --gte genotypephenotypecoupling
```

By default, the software tests all SNPs in your genotype data, having a minor allele frequency of > 0.05, a Hardy-Weinberg P-value > 0.001 and a call-rate > 0.95. If you want to test a subset of SNPs, create a text file (`snpfile`) with one column, one SNP identifier per row, and append the command above with the following command line switch `--snps snpfile` (remember to use the full path).

By default, the software uses 10 permutations to determine the False Discovery Rate (FDR) p-value threshold during the *cis*-eQTL mapping step. If you want to change the number of permutations (`nrperm`), you can append the command above with the following command line switch `--perm nrperm` (nrperm should be an integer).

If you are running the software in a cluster environment, you can specificy the number of threads to use (`nrthreads`) by appending the command above with the following command line switch `--threads nrthreads` (nrthreads should be an integer).

If you want to test all possible combinations in your dataset, you can append the command above using the following command line switch `--testall`.

If you want to use a set of QTLs that you have previously calculated (`eqtlfile`, for example from another dataset, which we do not recommend because of technical and biological differences between datasets), you can append the command above using the following command line switch `--eqtls eqtlfile` (remember to use the full path). The format of such file is described here: [File Formatas - eQTL file](link).

###Check your data
*MixupMapper* is a two stage approach. As such, the default procedure creates two directories in the `outdir` you specified: cis-eQTLs and MixupMapping. Both folders contain a different set of output files, described below.

####*cis*-eQTLs directory
This directory contains output from a default *cis*-eQTL mapping approach. The contents of this directory are detailed here: [QTL mapping output](link).

####MixupMapping directory

|File|Description|
|----------------|
|Distribution-Frequency.pdf|Frequency distribution of overall Z-scores. The Z-scores are plotted on the X-axis and the frequencies on the Y-axis.|
|Distribution-Frequency-BestMatches.pdf|This file is identical to Distribution-Frequency.pdf, except for the samples which have been removed after sample mix-up correction.|
|Heatmap.pdf|Visualization of overall Z-scores per assessed pair of samples. The genotyped samples are plotted on the X-axis, and the gene expression samples are plotted on the Y-axis. The brightness of each box corresponds to the height of the overall Z-score, with lower values having brighter colours. The grey bars next to the sample names indicate the correlation of the samples with the first Principal Component, which is an indicator for sample quality. Samples are sorted alphabetically on both axes.|
|MixupMapper-Dataset-&shy;ExpressionCorrelationMatrix.txt|Correlation matrix on the gene expression data that is used for calculation of the Principal Components (PCA). The PCA scores are visualizes in the Heatmap.pdf.|
|MixupMapper-Dataset-&shy;GenotypeCorrelationMatrix.txt|Correlation matrix of the genotype data that is used for calculation of the Principal Components (PCA). The correlation of the PCA scores for the genotype data are visualized in the Heatmap.pdf.|
|ROC.pdf|Receiver Operator Curve (ROC), describing the distance between the overall Z-scores of the diagonal in the heatmap, and all other assessed pairs. Several lines are plotted, each corresponding to an increasing a priori chance of mixed-up samples.|
|ROC.txt|A tab-separated text file representing the ROC curve.|
|SampleMixups.txt|A tab-separated text file describing the best matched (gene expression) sample for each genotyped sample.|
|ScoresMatrix.txt|A tab-separated text matrix of overall Z-scores as calculated by the method for each of the possible pairs of genotyped and gene expression samples.|
|SuggestedCouplings.txt|This file is identical to SampleMixups.txt, except for the lack of the header and the samples which have been removed after sample mix-up correction.|


In the SampleMixups.txt file, you can find the best matching expression sample for each genotyped sample:

-   1st column = genotyped sample ID
-	2nd column = best matching expression set 
-	3rd column = original genotype of the best matching expression set
-	4th column = concordant (true or false)
-	5th column = best match for expression set (true or false)
-	6th column = z-score (the lower, the better)
-	7th column = include sample? (true or false)

**Example of SampleMixups.txt**
<pre>
Genotype-ID    BM-Exp-ID	Ori-Gen-BM-Exp-ID	Concordant	BM-Exp	Z-Score	Include?
Sample1		Sample1		Sample1			true		true	-11.357	true
Sample2		Sample2		Sample2			true		true	-15.232	true
Sample3		Sample6		Sample6			false		true	-3.774	false
Sample4		Sample5		Sample4			false		false	-3.892	true
</pre>

###Resolving Sample mix-ups
As described above, the SampleMixups.txt describes possible sample mix-ups in your data. Resolving sample mix-ups is however a bit of a puzzle: many things could have happened during hybridization to the genotype and gene expression chips. For example, samples could have been duplicated (hybridized to the array twice), contaminated, swapped or not hybridized at all. The MixupMapper tries to resolve these issues automatically, although you should always check in your logbook whether the proposed sample mix-ups actually make sense (for example, are the mixed-up samples located on the same row or column of the chip, was a complete row inverted, or was the RNA quality poor). For each match, a z-score is presented which describes how much the samples are alike. You should interpret the z-score as a distance measure, so the lower the z-score, the better the match.

You should read the output of the program as follows: for every genotype sample, a single gene expression sample is matched (it is the gene expression sample with the lowest z-score for the genotype sample). The assessed genotype sample ID is located in column 1, the matching gene expression sample is in column 2. The third column describes what assignment you gave to the gene expression sample (eg: to which genotype sample you think the gene expression sample belongs). Now, if the original assigned genotype sample in column 3 is identical to the genotype sample in column 1, there is no indication of a sample mix-up, and column 4 will be TRUE. However, if the sample names in column 1 and 3 are different, column 4 is FALSE and something might be going on with the samples.

Consider the following example:
<pre>
GT-1    Ex-2	GT-2	FALSE	TRUE	-10.4	TRUE
GT-2	Ex-1	GT-1	FALSE	TRUE	-9.60	TRUE
</pre>	

The example described here is to be considered a classical sample swap: you observe that Ex-2 matches GT-1 best, and that gene Ex-1 matches GT-2 best. In this case, you can see that column 5 is also TRUE for both samples: for GT-1, this means that not only Ex-2 is the best match for GT-1, but also GT-1 is the best match for Ex1 (eg: the relationship is bidirectional). If we now observe that for GT-2 this relationship is also bidirectional, and the z-scores in column 6 are very low (eg: below -4, although this depends on the dataset), we get a strong indication that these samples are swapped.  

Now consider the following example:
<pre>
GT-1	Ex-2	GT-2	FALSE	TRUE	-10.4	FALSE
GT-2	Ex-2	GT-2	TRUE	FALSE	-9.60	TRUE
</pre>

In this case, Ex-2 is matched to two genotype samples. This means that either GT-1 and GT-2 are identical, or Ex-2 is contaminated (eg: a mix of RNA of both GT-1 and GT-2). The problem here is to decide which sample to include and which sample to exclude. The program decides here to stick with the original assignment provided by you, and exclude GT-1 (see column 7) even though the z-score is lower for the GT-1-Ex-2 match. The choice to exclude a sample should however not be made by a program: like described above, you should check whether things makes sense from the lab. The z-score can however give you an indication of how  much a gene expression sample resembles the genotype.

Finally, consider the following example:
<pre>
GT-1	Ex-7	GT-7	FALSE	TRUE	-10.4	TRUE
GT-2	Ex-6	GT-6	FALSE	TRUE	-9.60	TRUE
GT-3	Ex-5	GT-5	FALSE	TRUE	-10.4	TRUE
GT-4	Ex-4	GT-4	FALSE	TRUE	-9.60	TRUE
GT-5	Ex-3	GT-3	FALSE	TRUE	-10.4	TRUE
GT-6	Ex-2	GT-2	FALSE	TRUE	-9.60	TRUE
</pre>

The above example shows you an example of what a row inversion on a chip would look like: GT-1 matches Ex-7, GT-2 matches Ex-6, etcetera.

After checking the SampleMixups.txt file, some samples can be identified as sample mix-ups. You can easily replace the mixed IDs in the `genotypephenotypecoupling` file, without the need to change the ExpressionData files (adjusted or not) themselves. Please check in your plate- and array layouts whether it was possible to make these sample mix-ups. If your layouts don’t give you any information on why a sample could have been mixed-up we have chosen to exclude samples for which column 4 indicates FALSE (eg: where the match is not concordant to what was initially defined).

If you want to remove a sample after for example the mix-up step, remove the sample by either deleting it from your  `genotypephenotypecoupling`, by setting the sample to ‘exclude’ in the ‘PhenotypeInformation.txt’, or by removing the sample from the gene expression data.
**Make sure that you never remove lines from the Individuals.txt as this will result in erroneous genotypes for the remainder of the samples.**

Always rerun the Sample Mix-up Mapper after removing or changing sample IDs and check whether the results become better. 

##Step 5 - The optimum number of PCs to remove
Prior to eQTL mapping, we would like to determine whether removing PCs increases power to detect *cis*- and *trans*-QTLs. For each PC removing step (during normalization), this method runs both *cis*- and *trans*-eQTLs mapping(to refresh your memory: [see the section on normalization](link)). 

1. To be able to determine the optimum number of PCs to remove, and to reduce calculation time, we ordinarily run the *cis*-analysis on a selection of about 300.000 SNPs (Illumina HumanHap 300K content) and the *trans*-analysis on a selection of about 5.000 SNPs (content of the GWAS database: Genome.gov/GWAS/). However, other selections of SNPs can be used as well as detailed below. At the end of the analyses, the program will produce a list with the number of significant *cis*- and *trans*-QTLs, for each increment of PC removal. Based on this list, a table will be created showing the optimum number of PCs to remove. **This analysis is optional**.
2. Because we have seen that some PCs explain genetic variation, we perform an additional analysis to only adjust for PCs that have no genetic association. To identify the PCs that have no genetic association, we perform a QTL mapping on the principal component eigenvector matrix (PCAOverSamplesEigenvectorsTransposed.txt.gz). We call PCs genetically associated when they have a FDR of 0 (thus selecting truly significantly affected components only). Subsequently, we repeat the *cis*- and *trans*-QTL analyses, although this time we do not remove PCs that are genetically associated. Additionally, this second step also writes new PC Corrected phenotype files, although this time, excluding those compontents with a genetic association.

###Preparations
1. Note down the full path to your TriTyper genotype data directory. We will refer to this directory `genotypedir`. 
2. Determine the full path of the trait data you want to use, and make sure this data is normalized (e.g. use Quantile Normalized, Log<sub>2</sub> transformed Illumina gene expression data). We will refer to this path as `traitfile`. 
3. Locate your phenotype annotation file: `annotationfile`. Also note down the platform identifier `platformidentifier`.
4. Locate your `genotypephenotypecoupling` if you have such a file. You can also use this file to test specific combinations of genotype and phenotype individuals. 
5. Find a location on your hard drive to store the output. We will refer to this directory as `outputdir`.
6. Create (or download) a list of SNPs to use for the *cis*- and *trans*-QTL analyses, and save this in (a) text-file(s). Use a single column, and one line per SNP identifier. We will refer to these files as `cissnpfile` and `transsnpfile`. As with *MixupMapper*, SNPs will only be tested with a minor allele frequency of > 0.05, a Hardy-Weinberg P-value > 0.001 and a call-rate > 0.95.  

###Commands to be issued
To run the analysis, not taking into account the genetic association of PCs with SNPs, use the following command:


```
java –jar eQTLMappingPipeline.jar --mode pcaoptimum --in genotypedir --out outdir --inexp traitfile --inexpplatform platformidentifier --inexpannot annotationfile --gte genotypephenotypecoupling --cissnps cissnpfile --transsnps transsnpfile
```

If you want to run this analysis specifically for cis-QTLs, you can omit the `--transsnps transsnpfile` part of the command. Conversely, if you only want to run the trans-QTL analysis, you can omit the `--cissnps cissnpfile` part of the command.

To run the same analysis, taking the genetic association of PCs with SNPs into account (and to create the phenotype files that have been corrected with this approach), append the command above with the command line switch `--pcqtl`. 

If you are running the software in a cluster environment, you can specificy the number of threads to use (`nrthreads`) by appending the command above with the following command line switch `--threads nrthreads` (nrthreads should be an integer).

By default, the software uses 10 permutations to determine the False Discovery Rate (FDR) p-value threshold during the *cis*-eQTL mapping step. If you want to change the number of permutations (`nrperm`), you can append the command above with the following command line switch `--perm nrperm` (nrperm should be an integer).


###Check your data
After running the pcaoptimum command (both variants), the `outdir` will contain a number of directories from the performed QTL analyse. These folders will be named **Cis-nPCAsRemoved-GeneticVectorsNotRemoved** and **Trans-nPCAsRemoved-GeneticVectorsNotRemoved** (one folder per iteration of n PCs removed). The contents of these directories are detailed here: [QTL mapping output](link). If you have run the `--pcqtl` variant of this method, an additional folder will be created in the `outdir`, containing the QTL mapping on the PC eigenvectors. 

Additionally, at the end of the analyses, the program produces a table with the number of significant *cis*- and *trans*-QTLs for each increment of PC removal, the number of shared QTLs, and the number of QTLs with a different allelic direction compared to the number of detected QTLs using your input `traitfile` (this output will be printed on your screen). 

Furthermore, the program will produce a number of scatterplots in the `outdir` (x-axis: z-score for eQTLs after number of PCs removed, y-axis: z-score of eQTLs when 0 PCs removed) and some summary files. In excel, you can easily plot those PCA optimum numbers for both analyzes. 

Apart from QTL mapping results, the `--pcqtl` procedure will also produce new phenotype files in the same directory as your initial `traitfile`. These files have the suffix **-GeneticVectorsNotRemoved.txt.gz**. The format of these files is identical to the format of the input `traitfile` and can subsequently be used for QTL mapping.

##Step 6 - Perform the final QTL analysis
In a single command, the final QTL mapping can be performed. Standard settings for both *cis*- and *trans*-QTL mapping are: HWEP > 0.0001, MAF > 0.05, and Call Rate > 0.95. For *cis*-QTL mapping, the maximum distance between the SNP and the middle of the probe is 250.000bp. For *trans*-QTL mapping, the minimum distance between the SNP and the middle of the probe equals 5.000.000bp. To control for multiple testing, we perform 10 permutations, thereby shuffling sample labels, to calculate the false discovery rate (FDR) at 0.05 for both the *cis*- and *trans*-analysis.

###Preparations
1. Note down the full path to your TriTyper genotype data directory. We will refer to this directory `genotypedir`. 
2. Determine the full path of the trait data you want to use, and make sure this data is normalized (e.g. use Quantile Normalized, Log<sub>2</sub> transformed Illumina gene expression data, or any of the files produced in [Step 5](link)). We will refer to this path as `traitfile`. 
3. Locate your phenotype annotation file: `annotationfile`. Also note down the platform identifier `platformidentifier`.
4. Locate your `genotypephenotypecoupling` if you have such a file. You can also use this file to test specific combinations of genotype and phenotype individuals. 
5. Find a location on your hard drive to store the output. We will refer to this directory as `outputdir`.
6. (Optional) You can confine your analysis in several ways. For example, you can create a text-file with a list of snps (referred to as `snplist`; one SNP per line, single column), or a text file containing combination between SNPs and probes/traits/genes (referred to as `snpprobelist`; tab-separated, SNP identifier on first column, probe/gene/trait on second column).
7. (Optional) We have found that removing *cis* effects greatly enhances the power to detect *trans* effects. Consequently, our software provides a way to correct your `traitfile` data for *cis* effects. Note down the full path of the file containing the *cis*-QTL effects to be removed. We will refer to this file as `qtlfile`. The format of this file is identical to the [eQTLs.txt.gz file](link). 
    - **Note:** Removing QTL effects only properly works on imputed genotype data. The program will not remove QTL effects when imputation dosages are not available.

###Commands to be issued

The default command for a *cis*-QTL analysis:

```
java –jar eQTLMappingPipeline.jar --mode metaqtl --in genotypedir --out outdir --inexp traitfile --inexpplatform platformidentifier --inexpannot annotationfile --cis
```

By replacing the `--cis` command line switch with `--trans` in the above command, a genome-wide *trans* analysis is performed. 

If you want to run both *cis* and *trans* analyses at once, you can supply the software with both command line switches:

```
java –jar eQTLMappingPipeline.jar --mode metaqtl --in genotypedir --out outdir --inexp traitfile --inexpplatform platformidentifier --inexpannot annotationfile --cis --trans
```

If you want to confine your analysis to a certain set of SNPs, you can append the command with `--snps snplist`. Alternatively, you can supply a combination of SNPs and traits/probes/genes by appending `--snpprobe snpprobelist`.

If you want to remove QTL effects before your analysis, append the command with `--reqressouteqtls qtlfile`.

If you are running the software in a cluster environment, you can specificy the number of threads to use (`nrthreads`) by appending the command above with the following command line switch `--threads nrthreads` (nrthreads should be an integer).

By default, the software uses 10 permutations to determine the False Discovery Rate (FDR) p-value threshold during the *cis*-eQTL mapping step. If you want to change the number of permutations (`nrperm`), you can append the command above with the following command line switch `--perm nrperm` (nrperm should be an integer).

You can set the maximum number of results returned (`nrresults`) in the eQTLs.txt.gz file by appending the command above with `--maxresults nrresults`. The default is 500,000. Please note that increasing this number also increases the memory usage of the program and the file size of the output files.

By default, the program outputs results in a text-based format. However, for meta-analysis purposes, a binary format is also provided. You can switch to the binary format by appending `--binary` to the above commands. Tools to meta-analyze these binary files will be released at a later stage. You can also produce output in both text-based and binary formats, by appending `--text --binary` to the commands.

### QTL Mapping output - Text mode
This section describes the default output of a QTL mapping analysis. Depending on settings, one or more of the files below may not be present in your `outdir`.

|File|Description|
|----------------|
|DotPlot-FDR0.05.png|The x-axis depicts SNP chromosome position, the y-axis depicts the probe chromosome position. Each dot represents the position of a significant association.|
|eQTLProbesFDR0.05.txt|Contains eQTL results after correction for multiple testing. This file lists the strongest effect for each gene expression probe. FDR threshold value = 0.05.|
|eQTLSNPsFDR0.05.txt|Contains eQTL results after correction for multiple testing. This file lists the strongest effect for each SNP. FDR threshold value = 0.05.|
|eQTLs.txt|This file contains raw eQTL results on real data, not corrected for multiple testing: only the top 150,000 cis-eQTL results are shown, sorted by p-value.| 
|eQTLsFDR0.05.txt|This file contains all significant cis-eQTL effects. FDR threshold value = 0.05.|
|PermutedEQTLsPermutationRoundn.txt|This file contains eQTLs as determined based upon n permutation rounds. These files are used to calculate the FDR distribution.|
|eQTLsFDR0.05QQPlot.pdf|QQ plot of test statistics.|
|excludedSNPsBySNPFilter.txt|These SNPs are excluded by the user (by command)|
|excludedSNPsBySNPProbeCombinationFilter.txt|There are no gene expression probes located within 250.000bp of this SNP position.|

### QTL Mapping output - Binary mode
This section describes the binary output files generated using the `--binary` command line switch.

|File|Description|
|-----|
|Dataset.ProbeSummary.dat|This file contains the basic description of the tested probes: chromosome, chromosome position, probename and gene symbol. This file defines the order of probes in the ZscoreMatrix described below.|
|Dataset.SNPSummary.dat|This file contains basic information on the tested SNPs: chromosome, chromosome position, minor allele frequency, HWE, call-rate, minor allele, the assessed allele, the number of samples for which there were genotypes available, and the position in the ZscoreMatrix|
|Dataset.ZscoreMatrix.dat|This file contains a matrix the size of the number of probes in Dataset.ProbeSummary x Dataset.SnpSummary. The binary contents (1 floating point value per probe) are gzipped per SNP, which makes fast random access possible.|
|Dataset-PermutationRound-n.SNPSummary.dat|SNP Summary for the permuted data analysis|
|Dataset-PermutationRound-n.ZScoreMatrix.dat|Zscore matrix for the permuted data analysis|
|excludedSNPsBySNPFilter.txt.gz|These SNPs are excluded by the user.|
|excludedSNPsBySNPProbeCombinationFilter.txt.gz|There are no gene expression probes located within 250.000bp of this SNP position.|
|ProbeQCLog.txt.gz|This file describes which probes have been excluded from analysis.|
|SNPQCLog.txt.gz|This file contains QC information for all tested SNPs, including MAF, HWE-pvalue and call-rate (which is always 1.0 for imputed SNPs).|


#Additional analyses and advanced settings



##Multiple linear regression using interaction model

##Conditional analysis

##Meta-analysis and settings file
The command line interface of this software allows for basic QTL analyses. However, our software has many more capabilities that are not accesible via the command line. In these cases, an XML file is required that describes the different settings (full path referred to as `settingsfile`). [An example `settingsfile` is provided in the repository](link). Using a settings file allows you to quickly rerun certain analyses and to perform on-the-fly meta-analyses. A copy of the `settingsfile` will always be copied to your `outdir`. 

Currently, `settingsfile` can only be used in the `--mode metaqtl` mode. You should note that a `settingsfile` overrides all command line switches. The `settingsfile` can be used as follows:

```
java –jar eQTLMappingPipeline.jar --mode metaqtl --settings settingsfile
```

###Available options in settingsfile
XML is, like HTML, a hierarchical markup language, which works with so-called markup tags. An example of such tags (describing two QC settings) is below:

```
<defaults>
    <qc>
        <maf>0.05</maf>
        <hwep>0.001</hwep>    
    </qc>
</defaults>
```

Please note that if you open a tag "`<maf>`" you also need to close it: "`</maf>`". Also note that the `<maf>` tag is part of both `<qc>` and `<defaults>`, as an example of the hierarchy of XML files. If you have issues with the settings file, check whether all tags are opened and closed properly, and whether the hierarchy is correct. Also note that these tags are case-sensitive. 

Because of the readability of the table below, we reduce the above hierarchy to `defaults.qc.maf` and `defaults.qc.hwep`, respectively. 

|Setting|Value|Description|
|-------------------|
|sett.ing|double|description|

#File formats
This section lists the different file formats used by the software package.
##Probe annotation file
A probe annotation file is required when running a *cis*-eQTL analysis. This file describes where each probe/trait/gene is located on the genome. This file is a simple tab-separated text file, with a line for each probe, and a header.

**Please note:** for reliable meta-analysis, probe annotations should be identical between datasets. Therefore, we recommend use of the array address for your platform as a probe identifier if you are using Illumina array based expression data.

###File example
<pre>
Platform    HT12v4-ArrayAddress Symbol	Chr		ChrStart ChrEnd Probe     Seq
HT12v4      00001               GeneX	1       1504        1554        0       CGCTCCCCTTATAACTT-etc.
HT12v4      00002               GeneY	11      19900       19950       1       GGATCCCAGATTCCCT-etc.
HT12v4      00003               GeneZ	23      101         151         2       TTCTCCAGAGTCGAGC-etc.
</pre>


##Phenotype file, covariate file
Phenotype and covariate files have the same basic format. We use a tab separated text-based table, with individuals on columns and probes or covariates on the rows. 

###File example
<pre>
ProbeArrayAddress    Sample1	Sample2		Sample3
00001		    	1.640		1.553		1.441
00002		    	1.671		0.201		5.321
00003		    	1.126		1.710		2.569
00004		    	1.129		1.002		1.313
</pre>

**Please note:** if you are exporting your data from R, note that the first column often does not have a column identifier. Our software expects n+1 columns in both the header as well as the data (where n equals the number of individuals in your dataset).

##Genotype - phenotype coupling
Sometimes, the sample identifiers used in your phenotype data may not be identical to the identifiers used in your genotype data. Our software allows linking such identifiers through an external file we call a genotype-phenotype coupling file. This file also allows you to test specific combinations of genotype and phenotype individuals or can be used to select certain individuals to test by inclusion/exclusion.

The format of this file is a simple tab-delimited text file, one sample pair per row. The file should not contain duplicate entries. This file has no header. Please make sure the individual identifiers in the coupling file are identical to those in your genotype data and phenotype data.

###File example
<pre>
genotypesample1     phenotypesample1
genotypesample2     phenotypesample2
genotypesample3     phenotypesample3
</pre>

##TriTyper genotype data
The TriTyper format consists of several files, each describing an aspect of the genotype data. 

|File name|Required|Description|
|----------|-------|---------|
|**GenotypeMatrix.dat**|YES|Binary file containing genotype data. GenotypeMatrix.dat has the following file size: (number of SNPs x 2) x number of individuals. The ImputedDosageMatrix.dat should be half this size.|
|**ImputedDosageMatrix.dat**|NO|Binary file containing imputed genotype dosage values. The ImputedDosageMatrix.dat should be half this size in bytes compared to the GenotypeMatrix.dat.|
|**SNPs.txt**|YES|The list of SNPs that are encoded within the GenotypeMatrix.dat file. One line per SNP.|
|**SNPMappings.txt**|YES|The list of SNPs that are encoded within the GenotypeMatrix.dat file. One line per SNP, tab-separated: first column contains the chromosome number, second column contains the SNP position, and third column contains the SNPID (rs ID).|
|**Individuals.txt**|YES|The list of individuals that are encoded within the GenotypeMatrix.dat file. One line per individual. **Do not change the order of the individuals in this file, or the number of individuals in this file. You can change the individual identifiers, although duplicates are not allowed.**|
|**PhenotypeInformation.txt**|YES|This file describes the phenotypes of the individuals. One line per individual, 4 columns per individual: individual ID, case/control status, include/exclude a certain individual, gender (female/male). **This file does not have to contain all individuals contained in Individuals.txt and can be used to exclude certain individuals from the analysis**|

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
10    62765	rs11511647
10    84172	rs12218882
10    84426	rs10904045	
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
Sample1    control    include    female
Sample2    control    include    male
Sample3    control    include    female
</pre>



