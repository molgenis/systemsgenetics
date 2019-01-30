#QTL mapping pipeline information
This software allows for QTL mapping using linear models and direct meta-analysis of such data through weighted Z-score analysis.

##Questions or suggestions?
You can contact the authors of this software at westra.harmjan @ gmail.com, or lude @ ludesign.nl for questions or suggestions regarding the software or the manual.

##Manual contents

1. [Downloading the software](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#downloading-the-software)
2. [Before you start](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#before-you-start)
3. [Step by step QTL analysis](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-by-step-eqtl-analysis)
    * [Step 1 - Preparation phenotype data](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-1---preparation-phenotype-data)
    * [Step 2 - Preparation genotype data](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-2---preparation-of-genotype-data)
    * [File Checklist](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#file-checklist)
    * [Step 3 - Phenotype data normalization](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-3---phenotype-data-normalization)
    * [Step 4 - MixupMapper](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-4---mixupmapper)
    * [Step 5 - Determining the optimum number of PCs to remove](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-5---the-optimum-number-of-pcs-to-remove)
    * [Step 6 - Perform the final QTL analysis](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-6---perform-the-final-qtl-analysis)
4. [Additional analyses and advanced settings](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#additional-analyses-and-advanced-settings)
    * [Cell type specificity analysis](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#cell-type-specificity-analysis)
    * [Meta-analysis and settings file](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#meta-analysis-and-settings-file)
5. [File formats](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#file-formats)
    * [Probe annotation file](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#probe-annotation-file)
    * [Phenotype file, covariate file](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file)
    * [Genotype - phenotype coupling](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#genotype---phenotype-coupling)
    * [TriTyper genotype data](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#trityper-genotype-data)
6. [General software information](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#general-software-information) 
    * [QTL Mapping](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#qtl-mapping)
    * [Meta-analysis](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#meta-analysis)
    * [Multiple testing correction](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#multiple-testing-correction)
    * [MixupMapper](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#mixupmapper)
7. [Frequently asked questions]()

##Downloading the software
You can download the latest version of the software here: [Latest version](http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$eqtl-mapping-pipeline/lastBuild/).

Make sure to download `eqtl-mapping-pipeline-*.*.*-SNAPSHOT-dist.zip` or `eqtl-mapping-pipeline-*.*.*-SNAPSHOT-dist.tar.gz`

Please note that the manual refers to eqtl-mapping-pipeline.jar, while the name of the package described above may be different (because of different version numbers etc).

##Before you start

###Path definitions and commands
Please note that our software expects full paths, although shorter relative paths wil also work in most cases. So, if you are on Windows, a full path to a genotype directory would be similar to `c:\path\to\genotype\dir\` and a full path to a file would be `c:\path\to\genotype\directory\file.txt`. Linux and Mac OS use different path separators. On these systems, these paths would be similar to the following `/path/to/genotype/dir/` and `/path/to/genotype/dir/file.txt`. Our main point here is that when pointing to a directory, use a 'trailing slash'

This manual will combine references to paths with commands that need to be issued for a certain task. For example, at some point in this manual we refer to your phenotype data as `traitfile`, which will be printed in a grey box. Commands will also be in grey boxes, and can make references to paths defined earlier (to keep the manual readable), as follows:

````
java -jar eqtl-mapping-pipeline.jar --mode metaqtl --inexp traitfile
````

To run the command in the example above, you have to replace the string `traitfile` with the full path of your `traitfile` after the command line switch `--inexp`. So if the full path to your `traitfile` would be `/path/to/traitfile.txt`, the final command would be:

````
java -jar eqtl-mapping-pipeline.jar --mode metaqtl --inexp /path/to/traitfile.txt
````

###Java Virtual Machine
Our software is written in Java, which makes the software both fast and portable across multiple operating systems. Executing a java program is similar to executing a normal program or app, although there are some considerations:

* QTL mapping heavily relies on available memory. Make sure your machine is 64-bit and has lots of memory installed (at least 4Gb). 

* Please make sure your version of java is up-to-date. Use at least the 64-bit version of Java version 6 (also called version 1.6). The software should also work on version 7 and higher. If you are running on Windows, you can download the so-called Java Runtime Environment (JRE) from: http://www.java.com. If you are running Linux, the virtual machine may be present in the proprietary section of your package manager, or may be available via http://www.java.com.

* Java executables are called jar-files. The eQTL mapping pipeline is such a jar-file. You can execute it by using the following command from a terminal / console:
```    
    java –jar eqtl-mapping-pipeline.jar 
```

* You need to specify the maximal amount of memory available to the program using the command-line switch `–Xmx` (case-sensitive!). It is also wise to set an initial amount of memory available to the program which can be specified with the `-Xms` option (case-sensitive!). The amount of memory can be specified in megabytes (using an m suffix) or in gigabytes (using a g suffix). To be sure your computer is running java at 64-bit, please add the switch `–d64`. These java VM switches (`–Xmx`, `-Xms`, `–d64` and others) should be called prior to the `–jar` switch. To be sure you have enough space to put SNP and probe information in you should also set `-XX:StringTableSize=`, chosing a prime number which is slighly higher than the amount of SNPs and probes combined will yield optimal performance. An example command where the program is started with an allocation of 2gb of initialy memmory, a maximum of 4gb of memory is allowed and roughly 1000000 SNPs and probes are tested will look like this:

```    
    java –d64 -XX:StringTableSize=1000003 -Xms2g –Xmx4g –jar eqtl-mapping-pipeline.jar 
```

* Try to increase the `–Xmx` amount when you get Out-Of-Memory-Errors or errors involving ‘heap space’
* If the program return an " java.lang.OutOfMemoryError: PermGen space" error try adding `-XX:MaxPermSize=512m`, before the `-jar`. In the later releases of the software we use another representation of SNPs which needs another memory setting to be specified.

**IMPORTANT NOTE:** In this manual, we assume you understand the principle that you need to allocate sufficient amounts of RAM and therefore we excluded the VM switches from the example commands. Please be aware that you should use it, as some of the commands may require a substantial amount of memory (depending on your dataset and settings)!

###General information about software
* The eQTL mapping pipeline is a command line program, which makes the user interface not very intuitive. In order to help you a bit, an overview of available switch options is displayed when a command is incomplete or incorrect. Furthermore, each mode of the program also has its own overview of available switches with a small description of its functionality. For example: ```java –jar eqtl-mapping-pipeline.jar``` produces a list of available modes, while ```java –jar eqtl-mapping-pipeline.jar  --mode metaqtl``` produces a list of all available options for metaqtl.
* The software is able to process GZipped text files for most of the input files (not files in .tar archives however), which allows you to save some space on your hard drive.


#Step by step eQTL analysis

**This is not the anayslis plan for the eQTL meta analysis. Those instructions can be found here: https://github.com/molgenis/systemsgenetics/wiki/eQTL-mapping-analysis-cookbook**

This is a step by step guide which will guide you through the QTL mapping process using our software. Please note that this step by step guide illustrates only a part of the capabilities of our software. However, this guide does explain the different command line switches. [Additional analyses and advanced settings](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#additional-analyses-and-advanced-settings). Our general method consists of six steps described below:

* [Step 1 - Preparation phenotype data](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-1---preparation-phenotype-data)
* [Step 2 - Preparation genotype data](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-2---preparation-of-genotype-data)
* [File Checklist](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#file-checklist)
* [Step 3 - Phenotype data normalization](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-3---phenotype-data-normalization)
* [Step 4 - MixupMapper](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-4---mixupmapper)
* [Step 5 - Determining the optimum number of PCs to remove](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-5---the-optimum-number-of-pcs-to-remove)
* [Step 6 - Perform the final QTL analysis](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-6---perform-the-final-qtl-analysis)

##Definitions
Througout the manual, references to different full paths will be made. Here is an overview of these paths:

* The phenotype file will be referred to as `traitfile`
* The probe/trait/gene annotation file will be referred to as `annotationfile`
* The full path of your genotype data will be referred to as `genotypedir`
* The file linking phenotype individuals to genotype individuals will be referred to as `genotypephenotypecoupling`
* The file containing covariates will be defined as `covariatefile`

Descriptions of each of these files and their usage is detailed below, and their formats are described in the data formats section of this manual.

##Step 1 - Preparation phenotype data
Because our software uses a nonparametric test by default, you can use virtually any continuous data as trait values to map a variety of QTL effects. However, currently the normalization tools provided with this package are focused on array based methylation data, array based (Illumina) expression data, and preprocessed RNA-seq data (e.g. transcript level quantified data) or (GC)-RMA processed Affymetrix data. 

Please format your phenotype data in a simple tab separated text file. The format of this file is described here: [Phenotype file, covariate file](https://github.com/harmjanwestra/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file).

Some analyses will require an annotation of the probes/traits/genes in your `traitfile`. We store this annotation in a separate file. The format of this file is described here: [Probe annotation File](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#probe-annotation-file). 

###Illumina array based expression data

Because Illumina has several ways to annotate and preprocess their gene expression arrays, we provide a detailed instruction below, that will yield raw, untouched gene expression levels.

Use the GenomeStudio files of your expression arrays. **Note:** Please do not use any normalization, transformation, imputation, or background correction method that is offered by Illumina’s GenomeStudio! Just load your raw data, without doing any normalization, transformation, or background correction.

* Export a so-called ‘Final Report’ using GenomeStudio, including array_address_id (as the probe identifier) and the average expression signals per sample for all probes (click 'Analysis' in the top menu bar, then 'Reports').
* If you create more than one Final Report, merge the Final Reports, so that all data will be combined in one Expression Matrix File
* Rows should contain the different probes, and columns should contain the different sample IDs
* Remove any header information that Genome studio might produce: the header of the matrix is expected at the first row
* Create a new header (expected at the first row and column) for both rows and columns, describing the sample names (columns) and the probe IDs (rows). Note that you use the probe array address for your platform as the probe name.
* The final file should look like this format: [Phenotype file, covariate file](https://github.com/harmjanwestra/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file).


##Step 2 - Preparation of genotype data
Our software is able to use both unimputed called genotypes, as well as imputed genotypes and their dosage values. The software, however, requires these files to be in [TriTyper](http://genenetwork.nl/wordpress/trityper/) format. We provide [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) to harmonize and convert your genotype files. However, we first need to calculate MDS components to correct the expression data for population stratification.

###Correction for population stratification.

You will need to correct the expression data for 4 MDS components during the gene expression normalization step. To obtain these MDS components, you can use a simple command using the PLINK tool to create a multidimensional scaling matrix. The command line will look like this:

`plink --file mydata --cluster --mds-plot 4`

This command takes a single parameter, the number of dimensions to be extracted, 4 in this case.

This command creates the file plink.mds which contains one row per individual with the following fields:

FID		Family ID;
IID		Individual ID;
SOL		Assigned solution code;
C1		Position on first dimension;
C2 		Position on second dimension;
C3		Position on third dimension;
C4		Position on fourth dimension;

You can find more information on [PLINK website](http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml)

###Convert data to the TriTyper fileformat.
For the QTL mapping we need to data to  be in TriTyper format. Using [GenotypeHarmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype%20Harmonizer%20Download) one can harmonize and convert genotype formats. In this step we will directly harmonize, filter and convert the genotype data. The data is harmonized to all be matched to the GIANT release of 1000G. You can either run this step per CHR, or for all chromosomes at once.

`java -jar ./GenotypeHarmonizer.jar -i {locationOfInputData} -I {InputType} -o{Outputlocation} -r {locationOf100G-GaintVcfs} --refType VCF_FOLDER --update-id --keep -cf 0.95 -hf 0.0001 -mf 0.05`

Details on the flags: `-i` location of the input data, `-I` input type, `-o` output location, `-r` location of the reference data, `--refType` is the reference data, `-cf` is used for the callrate filter, `-hf` is used for the hardy weinberg equilibrium threshold,`-mf` is used for minor allele frequency threshold. See the Genotype Harmonization manual for further [details](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer#arguments-overview) on the input flags. Please make sure you use these exact parameters for the eQTL meta-analysis (-cf 0.95 -hf 0.0001 -mf 0.05).

**Check and send the log of Genotype Harmonizer!**

If you ran the GH per chromosome you now need to merge all TriTyper folders per CHR to one TriTyper folder containig data from all chromosomes. By using the following comman you merge the individual TriTyper folders.

`java -jar ./eqtl-mapping-pipeline.jar --imputationtool --mode concat --in {folder1;folder2;folder3;ETC} --out {OutputFolder} `
 
`--in` is a tab-separated list of input TriTyper Folders with information per chromosome, `--out` is the output location of the merged TriTyper data.

You now have your genotype data ready to go!

###Check your data
After converting your genotype data to TriTyper format, a number of files should have been created in your output directory. For a description of the different files, please refer to [TriTyper genotype data](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#trityper-genotype-data)

#File Checklist
Before you continue in this manual, this is a good time check whether your files are in the correct format and whether you have all the required files ready.

* Check whether all required files are in your `genotypedir`.
* Check whether the `traitfile` is properly normalized.
* Check whether you have an `annotationfile`. The format is described here [Probe annotation File](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#probe-annotation-file). 
    - Some QTL data may not have annotation (for example if you want to do a GWAS on lipid levels). In such cases an `annotationfile` is not required for QTL mapping. Please see [Step 6 - Perform the final QTL analysis](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-6---perform-the-final-qtl-analysis) for more details. 
    - For the other steps, the annotation file is required, although you can fool the software by setting your genotype genomic locations to 1 (both chromosome and chromosome position) and doing the same for your phenotype data using the `annotationfile`.
* If the sample identifiers differ between genotype and phenotype data, you have to create a file that links these identifiers together. This format is described here: [Genotype - phenotype coupling](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#genotype---phenotype-coupling). This file is optional, although we will mention this file in this manual as `genotypephenotypecoupling`

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
Note down the full path to your `traitfile`. The output of the normalization tool will be written to the same directory by default. Optionally, note down the full path to the file containing covariates. We will refer to the covariate file as `covariatefile`. The format of this file is identical to the `traitfile`: [Phenotype file, covariate file](https://github.com/harmjanwestra/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file). 

###Commands to be issued
To run the general normalization strategy described above, you can run the following command:

```
java –jar eqtl-mapping-pipeline.jar --mode normalize --in traitfile
```

You can specify an output directory with the following command (specifying an `outdir`):

```
java –jar eqtl-mapping-pipeline.jar --mode normalize --in traitfile --out outdir
```

To run the general normalization strategy described above, and correct for covariates:

```
java –jar eqtl-mapping-pipeline.jar --mode normalize --in traitfile --adjustcovariates --cov covariatefile
```

Individual elements of the normalization strategy can also be separately executed. For example to only run Quantile normalization, and Log<sub>2</sub> transformation run the following command:

```
java –jar eqtl-mapping-pipeline.jar --mode normalize --in traitfile --qqnorm --logtransform
```

**Note:**
Several other parameters can be set to customize your normalization strategy (e.g. which procedures to run, number of PCs to remove, step size for PC removal, handling of missing values, etc). However, the order of procedures is fixed (Quantile Normalize > Log<sub>2</sub> transform > covariate adjustment > centering and scaling > PCA adjustment), irregardless of the order of each command line switch. To review the available options for normalization, issue the following command:

```
java –jar eqtl-mapping-pipeline.jar --mode normalize
```

###Check your data
Running the general normalization procedure yields a number of files in the directory of your `traitfile`, or in the `outdir` if you specified one. The default procedure will generate files suffixes listed below. Suffixes will be appended in the default order as described above. Selecting multiple normalization methods will add multiple suffixes. 
 

|Suffix|Description|
|----|-----------|
|**QuantileNormalized**|Quantile Normalized trait data|
|**Log2Transformed**|Log<sub>2</sub> Transformed trait data.|
|**ProbesCentered**|Probes were centered.|
|**SamplesZTransformed**|Samples were Z-transformed.|
|**CovariatesRemoved**|Gene expression was adjusted for covariates.|
|**PCAOverSamplesEigenvalues**|Eigenvalues created during eigenvalue decomposition of the gene expression sample correlation matrix (created from the Quantile Normalized, Log2 Transformed, Z-transformed data)|
|**PCAOverSamplesEigenvectors**|Eigenvectors created during eigenvalue decomposition of the gene expression sample correlation matrix (created from the Quantile Normalized, Log2 Transformed, Z-transformed data)|
|**PCAOverSamplesEigenvectorsTransposed**|Eigenvectors transposed|
|**PCAOverSamplesPrincipalComponents**|Principal Components describing the sample correlation matrix (created from the Quantile Normalized, Log2 Transformed, Z-transformed data)|
|**nPCAsOverSamplesRemoved**|Expression data, Quantile Normalized, Log2 Transformed, Z-transformed, with n Principal Components regressed out.|


##Step 4 - MixupMapper
We have shown in a paper published in Bioinformatics ([Westra et al.: *MixupMapper: correcting sample mix-ups in genome-wide datasets increases power to detect small genetic effects*](http://bioinformatics.oxfordjournals.org/content/27/15/2104)), that sample mix-ups often occur in genetical genomics datasets (i.e. datasets with both genotype and gene expression data). Therefore, we developed a method called *MixupMapper*, which is implemented in the eQTL Mapping Pipeline. This program performs the following steps:
1. At first a *cis*-eQTL analysis is conducted on the dataset:
    -	using a 250 kb window between the SNP and the mid-probe position
    -	performing 10 permutations to control the false discovery rate (FDR) at 0.05
2. Calculate how well a gene expression array matches a genotype array. 
For details how this exactly works, please have a look at the paper or read the attachment. 

###Preparations/variables
1. Note down the full path to your TriTyper genotype data directory. We will refer to this directory `genotypedir`. 
2. Determine the full path of the trait data you want to use, and make sure this data is normalized (e.g. use Quantile Normalized, Log<sub>2</sub> transformed Illumina gene expression data, that has been corrected for eventual covariates.). We will refer to this path as `traitfile`. 
3. Locate your phenotype annotation file: `annotationfile`. Also note down the platform identifier `platformidentifier`.
4. (Optional) Locate your `genotypephenotypecoupling` if you have such a file. You can also use this file to test specific combinations of genotype and phenotype individuals. 
5. Find a location on your hard drive to store the output. We will refer to this directory as `outputdir`.

###Commands to be issued
The *MixupMapper* analysis can be run using the following command:

```
java –jar eqtl-mapping-pipeline.jar --mode mixupmapper --in genotypedir --out outdir --inexp traitfile --inexpplatform platformidentifier --inexpannot annotationfile --gte genotypephenotypecoupling
```

By default, the software tests all SNPs in your genotype data, having a minor allele frequency of > 0.05, a Hardy-Weinberg P-value > 0.001 and a call-rate > 0.95. If you want to test a subset of SNPs, create a text file (`snpfile`) with one column, one SNP identifier per row, and append the command above with the following command line switch `--snps snpfile` (remember to use the full path).

By default, the software uses 10 permutations to determine the False Discovery Rate (FDR) p-value threshold during the *cis*-eQTL mapping step. If you want to change the number of permutations (`nrperm`), you can append the command above with the following command line switch `--perm nrperm` (nrperm should be an integer).

If you are running the software in a cluster environment, you can specificy the number of threads to use (`nrthreads`) by appending the command above with the following command line switch `--threads nrthreads` (nrthreads should be an integer).

If you want to test all possible combinations in your dataset, you can append the command above using the following command line switch `--testall`.

Note that the `--gte` switch is optional, and only applies if you are using a `genotypephenotypecoupling` file.

If you want to use a set of QTLs that you have previously calculated (`eqtlfile`, for example from another dataset, which we do not recommend because of technical and biological differences between datasets), you can append the command above using the following command line switch `--eqtls eqtlfile` (remember to use the full path). The format of such file is described here: [QTL file](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#qtl-file).

###Check your data
*MixupMapper* is a two stage approach. As such, the default procedure creates two directories in the `outdir` you specified: *cis*-eQTLs and MixupMapping. Both folders contain a different set of output files, described below.

####*cis*-eQTLs directory
This directory contains output from a default *cis*-eQTL mapping approach. The contents of this directory are detailed here: [QTL Mapping output - Text mode](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#qtl-mapping-output---text-mode).

####MixupMapping directory


|File|Description|
|----|-----------|
|Heatmap.pdf|Visualization of overall Z-scores per assessed pair of samples. The genotyped samples are plotted on the X-axis, and the gene expression samples are plotted on the Y-axis. The brightness of each box corresponds to the height of the overall Z-score, with lower values having brighter colours. Samples are sorted alphabetically on both axes.|
|BestMatchPerGenotype.txt|This file shows the best matching trait samples per genotype: the result matrix (MixupMapperScores.txt) is not symmetrical. As such, scanning for the best sample per genotype may yield other results than scanning for the best sample per trait.|
|BestMatchPerTrait.txt|This file shows the best matching genotype sample per trait sample.|
|MixupMapperScores.txt|A matrix showing the scores per pair of samples (combinations of traits and genotypes). |

In the BestMatchPerGenotype.txt and BestMatchPerTrait.txt files, you can find the best matching trait sample for each genotyped sample and vice versa:

*   1st column = genotyped sample ID, or trait sample ID dependent on file chosen (see above)
*	2nd column = trait sample originally linked to genotype sample ID in column 1, or genotype sample originally linked with trait sample in column 1
*	3rd column = the MixupMapper Z-score for the link between the samples in column 1 and 2
*	4rd column = best matching trait (for BestMatchPerGenotype.txt) or best matching genotype (for BestMatchPerTrait.txt)
*	5th column = the MixupMapper Z-score for the link between the samples in column 1 and 4
*	6th column = this column determines whether the best matching trait or genotype is identical to the sample found in column 2

**Example of BestMatchPerGenotype.txt**
<pre>
Trait   OriginalLinkedGenotype  OriginalLinkedGenotypeScore     BestMatchingGenotype    BestMatchingGenotypeScore       Mixup
Sample1		Sample1		-11.357			Sample1		-11.357	false
Sample2		Sample2		-15.232			Sample2		-15.232	false
Sample3		Sample3		-3.774			Sample5		-6.341	true
Sample4		Sample4		-3.892			Sample6		-12.263	true
</pre>

###Resolving Sample mix-ups
As described above, the BestMatchPerGenotype.txt describes possible sample mix-ups in your data. Resolving sample mix-ups is however a bit of a puzzle: many things could have happened during hybridization to the genotype and gene expression chips. For example, samples could have been duplicated (hybridized to the array twice), contaminated, swapped or not hybridized at all. The MixupMapper tries to resolve these issues automatically, although you should always check in your logbook whether the proposed sample mix-ups actually make sense (for example, are the mixed-up samples located on the same row or column of the chip, was a complete row inverted, or was the RNA quality poor). For each match, a z-score is presented which describes how much the samples are alike. You should interpret the z-score as a distance measure, so the lower the z-score, the better the match.

You should read the output of the program as follows: for every genotype sample, a single gene expression sample is matched (it is the gene expression sample with the lowest z-score for the genotype sample). The assessed genotype sample ID is located in column 1, the matching gene expression sample is in column 4. The second column describes what assignment you gave to the gene expression sample (eg: to which genotype sample you think the gene expression sample belongs). Now, if the original assigned genotype sample in column 4 is identical to the genotype sample in column 1, there is no indication of a sample mix-up, and column 6 will be FALSE. However, if the sample names in column 1 and 4 are different, column 6 is TRUE and something might be going on with the samples.

Consider the following example:
<pre>
GT-1    Ex-1	-3.4	Ex-2	-10.6	TRUE
GT-2	Ex-2	-2.60	Ex-1	-9.5	TRUE
</pre>	

The example described here is to be considered a classical sample swap: you observe that Ex-2 matches GT-1 best, and that gene Ex-1 matches GT-2 best. In this case, you can see that column 6 is also TRUE for both samples: for GT-1, this means that not only Ex-2 is the best match for GT-1, but also GT-1 is the best match for Ex1 (eg: the relationship is bidirectional). If we now observe that for GT-2 this relationship is also bidirectional, and the z-scores in column 5 are very low (eg: below -4, although this depends on the dataset), we get a strong indication that these samples are swapped.  

Now consider the following example:
<pre>
GT-1    Ex-1	-3.4	Ex-2	-10.6	TRUE
GT-2	Ex-2	-9.5	Ex-2	-9.5	TRUE
</pre>

In this case, Ex-2 is matched to two genotype samples. This means that either GT-1 and GT-2 are identical, or Ex-2 is contaminated (eg: a mix of RNA of both GT-1 and GT-2). The problem here is to decide which sample to include and which sample to exclude. In such a case, the best choice will often be to stick with the original assignment provided by you, and exclude GT-1 (see column 7) even though the z-score is lower for the GT-1-Ex-2 match. The choice to exclude a sample should however not be made by a program: like described above, you should check whether things makes sense from the lab. The z-score can however give you an indication of how  much a gene expression sample resembles the genotype.

Finally, consider the following example:
<pre>
GT-1	Ex-1	-4.3	Ex-6	-10.4	TRUE
GT-2	Ex-2	-4.2	Ex-5	-9.60	TRUE
GT-3	Ex-3	-5.3	Ex-4	-10.4	TRUE
GT-4	Ex-4	-4.3	Ex-3	-9.60	TRUE
GT-5	Ex-5	-7.5	Ex-2	-10.4	TRUE
GT-6	Ex-6	-2.3	Ex-1	-9.60	TRUE
</pre>

The above example shows you an example of what a row inversion on a chip would look like: GT-1 matches Ex-6, GT-2 matches Ex-5, etcetera.

After checking the BestMatchPerGenotype.txt or BestMatchPerTrait.txt files, some samples can be identified as sample mix-ups. You can easily replace the mixed IDs in the `genotypephenotypecoupling` file, without the need to change the ExpressionData files (adjusted or not) themselves. Please check in your plate- and array layouts whether it was possible to make these sample mix-ups. If your layouts don’t give you any information on why a sample could have been mixed-up our best practice is to remove the sample pair from further analysis.

If you want to remove a sample after for example the mix-up step, remove the sample by either deleting it from your `genotypephenotypecoupling`, by setting the sample to ‘exclude’ in the ‘PhenotypeInformation.txt’, or by removing the sample from the gene expression data.
**Make sure that you never remove lines from the Individuals.txt in your TriTyper folder, as this will result in erroneous genotypes for the remainder of the samples.**

Always rerun the Sample Mix-up Mapper after removing or changing sample IDs and check whether the results become better. 

##Step 5 - The optimum number of PCs to remove

**This is not the anayslis plan for the eQTL meta analysis. Those instructions can be found here: https://github.com/molgenis/systemsgenetics/wiki/eQTL-mapping-analysis-cookbook**

Prior to eQTL mapping, we would like to determine whether removing PCs increases power to detect *cis*- and *trans*-QTLs. For each PC removing step (during normalization), this method runs both *cis*- and *trans*-eQTLs mapping(to refresh your memory: [Step 3 - Phenotype data normalization](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-3---phenotype-data-normalization)). 

1. To be able to determine the optimum number of PCs to remove, and to reduce calculation time, we ordinarily run the *cis*-analysis on a selection of about 300.000 SNPs (Illumina HumanHap 300K content) and the *trans*-analysis on a selection of about 5.000 SNPs (content of the GWAS database: Genome.gov/GWAS/). However, other selections of SNPs can be used as well as detailed below. At the end of the analyses, the program will produce a list with the number of significant *cis*- and *trans*-QTLs, for each increment of PC removal. Based on this list, a table will be created showing the optimum number of PCs to remove. **This analysis is optional**.
2. Because we have seen that some PCs explain genetic variation, we perform an additional analysis to only adjust for PCs that have no genetic association. To identify the PCs that have no genetic association, we perform a QTL mapping on the principal component eigenvector matrix (PCAOverSamplesEigenvectorsTransposed.txt.gz). We call PCs genetically associated when they have a FDR of 0 (thus selecting truly significantly affected components only). Subsequently, we repeat the *cis*- and *trans*-QTL analyses, although this time we do not remove PCs that are genetically associated. Additionally, this second step also writes new PC Corrected phenotype files, although this time, excluding those compontents with a genetic association.

###Preparations
1. Note down the full path to your TriTyper genotype data directory. We will refer to this directory `genotypedir`. 
2. Determine the full path of the trait data you want to use, and make sure this data is normalized (e.g. use Quantile Normalized, Log<sub>2</sub> transformed Illumina gene expression data). We will refer to this path as `traitfile`. 
3. Locate your phenotype annotation file: `annotationfile`. Also note down the platform identifier `platformidentifier`.
4. (Optional) Locate your `genotypephenotypecoupling` if you have such a file. You can also use this file to test specific combinations of genotype and phenotype individuals. 
5. Find a location on your hard drive to store the output. We will refer to this directory as `outputdir`.
6. Create (or download) a list of SNPs to use for the *cis*- and *trans*-QTL analyses, and save this in (a) text-file(s). Use a single column, and one line per SNP identifier. We will refer to these files as `cissnpfile` and `transsnpfile`. As with *MixupMapper*, SNPs will only be tested with a minor allele frequency of > 0.05, a Hardy-Weinberg P-value > 0.001 and a call-rate > 0.95. Examples of these files can be found here: [resources](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline/src/main/resources)

###Commands to be issued
To run the analysis, not taking into account the genetic association of PCs with SNPs, use the following command:

```
java –jar eqtl-mapping-pipeline.jar --mode pcaoptimum --in genotypedir --out outdir --inexp traitfile --inexpplatform platformidentifier --inexpannot annotationfile --gte genotypephenotypecoupling --cissnps cissnpfile --transsnps transsnpfile
```

If you want to run this analysis specifically for *cis*-QTLs, you can omit the `--transsnps transsnpfile` part of the command. Conversely, if you only want to run the *trans*-QTL analysis, you can omit the `--cissnps cissnpfile` part of the command.

To run the same analysis, taking the genetic association of PCs with SNPs into account (and to create the phenotype files that have been corrected with this approach), append the command above with the command line switch `--pcqtl`. 

If you are running the software in a cluster environment, you can specificy the number of threads to use (`nrthreads`) by appending the command above with the following command line switch `--threads nrthreads` (nrthreads should be an integer).

By default, the software uses 10 permutations to determine the False Discovery Rate (FDR) p-value threshold during the *cis*-eQTL mapping step. If you want to change the number of permutations (`nrperm`), you can append the command above with the following command line switch `--perm nrperm` (nrperm should be an integer).

Note that the `--gte` switch only applies if you are using a `genotypephenotype` coupling file.

###Check your data
After running the pcaoptimum command (both variants), the `outdir` will contain a number of directories from the performed QTL analyse. These folders will be named **Cis-nPCAsRemoved-GeneticVectorsNotRemoved** and **Trans-nPCAsRemoved-GeneticVectorsNotRemoved** (one folder per iteration of n PCs removed). The contents of these directories are detailed here: [QTL Mapping output - Text mode](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#qtl-mapping-output---text-mode). If you have run the `--pcqtl` variant of this method, an additional folder will be created in the `outdir`, containing the QTL mapping on the PC eigenvectors. 

Additionally, at the end of the analyses, the program produces a table with the number of significant *cis*- and *trans*-QTLs for each increment of PC removal, the number of shared QTLs, and the number of QTLs with a different allelic direction compared to the number of detected QTLs using your input `traitfile` (this output will be printed on your screen). 

Furthermore, the program will produce a number of scatterplots in the `outdir` (x-axis: z-score for eQTLs after number of PCs removed, y-axis: z-score of eQTLs when 0 PCs removed) and some summary files. In excel, you can easily plot those PCA optimum numbers for both analyzes. 

Apart from QTL mapping results, the `--pcqtl` procedure will also produce new phenotype files in the same directory as your initial `traitfile`. These files have the suffix **-GeneticVectorsNotRemoved.txt.gz**. The format of these files is identical to the format of the input `traitfile` and can subsequently be used for QTL mapping.

##Step 6 - Perform the final QTL analysis

**This is not the anayslis plan for the eQTL meta analysis. Those instructions can be found here: https://github.com/molgenis/systemsgenetics/wiki/eQTL-mapping-analysis-cookbook**

In a single command, the final QTL mapping can be performed. Standard settings for both *cis*- and *trans*-QTL mapping are: HWEP > 0.0001, MAF > 0.05, and Call Rate > 0.95. For *cis*-QTL mapping, the maximum distance between the SNP and the middle of the probe is 250.000bp. For *trans*-QTL mapping, the minimum distance between the SNP and the middle of the probe equals 5.000.000bp. To control for multiple testing, we perform 10 permutations, thereby shuffling sample labels, to calculate the false discovery rate (FDR) at 0.05 for both the *cis*- and *trans*-analysis.

###Preparations
1. Note down the full path to your TriTyper genotype data directory. We will refer to this directory `genotypedir`. 
2. Determine the full path of the trait data you want to use, and make sure this data is normalized (e.g. use Quantile Normalized, Log<sub>2</sub> transformed Illumina gene expression data, or any of the files produced in [Step 5 - Determining the optimum number of PCs to remove](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-5---the-optimum-number-of-pcs-to-remove)). We will refer to this path as `traitfile`. 
3. Locate your phenotype annotation file: `annotationfile`. Also note down the platform identifier `platformidentifier`.
4. (Optional) Locate your `genotypephenotypecoupling` if you have such a file. You can also use this file to test specific combinations of genotype and phenotype individuals. 
5. Find a location on your hard drive to store the output. We will refer to this directory as `outputdir`.
6. (Optional) You can confine your analysis in several ways. For example, you can create a text-file with a list of snps (referred to as `snplist`; one SNP per line, single column), or a text file containing combination between SNPs and probes/traits/genes (referred to as `snpprobelist`; tab-separated, SNP identifier on first column, probe/gene/trait on second column).
7. (Optional) We have found that removing *cis* effects greatly enhances the power to detect *trans* effects. Consequently, our software provides a way to correct your `traitfile` data for *cis* effects. Note down the full path of the file containing the *cis*-QTL effects to be removed. We will refer to this file as `qtlfile`. The format of this file is identical to the [QTL file](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#qtl-file). 
    - **Note:** Removing QTL effects only properly works on imputed genotype data. The program will not remove QTL effects when imputation dosages are not available.

###Commands to be issued

The default command for a *cis*-QTL analysis:

```
java –jar eqtl-mapping-pipeline.jar --mode metaqtl --in genotypedir --out outdir --inexp traitfile --inexpplatform platformidentifier --inexpannot annotationfile --cis
```

By replacing the `--cis` command line switch with `--trans` in the above command, a genome-wide *trans* analysis is performed. 

If you want to run both *cis* and *trans* analyses at once, you can supply the software with both command line switches:

```
java –jar eqtl-mapping-pipeline.jar --mode metaqtl --in genotypedir --out outdir --inexp traitfile --inexpplatform platformidentifier --inexpannot annotationfile --cis --trans
```

If you want to confine your analysis to a certain set of SNPs, you can append the command with `--snps snplist`. Alternatively, you can supply a combination of SNPs and traits/probes/genes by appending `--snpprobe snpprobelist`.

If you want to remove QTL effects before your analysis, append the command with `--reqressouteqtls qtlfile`.

If you are running the software in a cluster environment, you can specificy the number of threads to use (`nrthreads`) by appending the command above with the following command line switch `--threads nrthreads` (nrthreads should be an integer).

By default, the software uses 10 permutations to determine the False Discovery Rate (FDR) p-value threshold during the *cis*-eQTL mapping step. If you want to change the number of permutations (`nrperm`), you can append the command above with the following command line switch `--perm nrperm` (nrperm should be an integer).

You can set the maximum number of results returned (`nrresults`) in the eQTLs.txt.gz file by appending the command above with `--maxresults nrresults`. The default is 500,000. Please note that increasing this number also increases the memory usage of the program and the file size of the output files.

Note that you should apply the `--gte` switch if you are using a `genotypephenotype` coupling file.

By default, the program outputs results in a text-based format. However, for meta-analysis purposes, a binary format is also provided. You can switch to the binary format by appending `--binary` to the above commands. Tools to meta-analyze these binary files will be released at a later stage. You can also produce output in both text-based and binary formats, by appending `--text --binary` to the commands.

### QTL Mapping output - Text mode
This section describes the default output of a QTL mapping analysis. Depending on settings, one or more of the files below may not be present in your `outdir`.

|File|Description|
|----|-----------|
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
|-----|----------|
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



##Cell type specificity analysis
###Background
Many eQTL datasets are created by measuring gene expression tissues that consist of many different cell types. As about 40% of the trait-associated SNPs show a cis-eQTL effect in these compound tissues, it is often unclear what the causal cell-type is for the disease. The method described here is able to use an a priori defined list of cell-type specific genes to determine whether a cis-eQTL is specific for a cell type.

This method is part of a manuscript titled:
**Cell-type specific eQTL analysis without the need to sort cells**

###Method
Provided a list of genes that show a high correlation with the cell-count for the cell-type in question, we can determine a proxy phenotype in another gene expression dataset that does not have cell counts. Our method first creates a correlation matrix for the list of cell-type specific genes using the raw gene expression data (quantile normalized, log2 transformed, corrected for MDS components). Using principal component analysis on the correlation matrix, we obtain principal components (PC) that describe the variation among these genes. The first PC (PC1) attributes the largest amount of variation, and can thus be seen as a proxy-phenotype for the cell-type. PC1 can then be used as an independent variable in the linear model for the cis-eQTL. This means that apart from the genotype effect, we can now also determine the effect of PC1 on gene expression. Apart from these two effects, we also determine the interaction term between genotype and PC1. This interaction term describes the dependence of the cis-eQTL on PC1, and thus effectively the interaction between the cis-eQTL and the cell type (see the figure on the next page). We choose to only use PC1 here, because PC1 describes the majority of the variation amongst the cell-type specific genes and also to limit the number of independent variables in our model.
A typical cis-eQTL linear model is the following:

`y ~ β*g + e`

where y is the gene expression of the gene, and β is the slope of the linear model, g is the genotype and e is the intersect with the y-axis. Including PC1 as an independent variable, the model is as follows:

`y ~ β0*g + β1*p + β2*p*g + e`

where p is the PC1 covariate and p*g is the interaction term between PC1 and the genotype. Note that this model fits three linear models at once (three different slopes). Apart from the cell count proxy, we also determine the interaction term for all genes that show a cis-eQTL effect.


###Method overview
The method is a two-step process:
One step performs initial normalization of the gene expression data and adds an extra quality control step (by correlating the first PC over the sample correlation matrix, samples with poor RNA quality can be determined. We use a correlation threshold of 0.9 with PC1 to remove poor quality samples). Then, the program calculates the first PC using the cell-type specific probe correlation matrix to use as a proxy phenotype in the second step of the program. The second step performs the actual eQTL mapping using the Ordinary Least Squares model with two independent variables and an interaction term.

####Normalization - Step 1 - Prepare your data
Locate and/or download the following files (to avoid confusion, **use full paths** when supplying these files to the software):

-	ExpressionData.txt.gz: the raw gene expression data (not normalized, not corrected for MDS components, not corrected for any *cis*-eQTL effects). We will call this file `traitfile`.
-	Find a location on your hard drive to store the output. We will refer to this directory as `outputdir`.
-	CellTypeSpecificProbeList.txt: the list of probes that correlate with the cell type of interest. One probe per line. Make sure the probe identifiers (or a subset thereof) are consistent with those in your `traitfile` We will call this file `probelist`.
-	MDSComponents.txt: this (optional) file contains the four MDS components that were previously used for trait data normalization (which were calculated using pruned genotypes). We will call this file `mdscomponents`. This format of this file is described here [Phenotype file, covariate file](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file).
-	GenotypeToExpressionCoupling.txt: the file linking genotypes to gene expression sample identifiers. This file is also optional (for example if your sample identifiers in your gene expression data correspond to the genotype sample identifiers). We will call this file `genotypetotraitcoupling`. The format is described here: [Genotype - phenotype coupling](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#genotype---phenotype-coupling).
-	CellCounts.txt: in order to determine whether PC1 actually reflects cell-type differences, we may ask you to supply this optional file. The program will determine and output the correlation between the proxy phenotype and the actual phenotype if this file is supplied. The format is identical to the `mdscomponents` file. We will call this file `cellcounts`

####Normalization - Step 2 - Run the normalization and QC
``java –jar eqtl-mapping-pipeline.jar --mode celltypespecific --step normalize --inexpraw traitfile --out outputdir --celltypespecificprobes probelist --mdscomponents mdscomponents --gte genotypetotraitcoupling --cellcounts cellcounts``

**Please note that the --gte, --mdscomponents and --cellcount switches are optional** 

####Normalization - Step 3 - Check the output
After completion, the normalization step will have generated a number of files in your `outputdir`:

|File|Description|
|----|-----------|
|CellTypeProxyFile.txt|This file contains the first principal component, calculated over the cell type specific probe correlation matrix. This file is required for step 2.|
|CellTypeSpecificProbeCorrelationMatrix.txt.gz|This file contains the correlation matrix for the cell type specific probes|
|CellTypeSpecificProbePCA.PCAoverSamplesEigenvalues.txt.gz|This file contains the eigenvalues for the PCA analysis over the cell type specific probes|
|CellTypeSpecificProbePCA.PCAoverSamplesEigenvectors.txt.gz|This file contains the eigenvectors for the PCA analysis over the cell type specific probes|
|CellTypeSpecificProbePCA.PCAoverSamplesEigenvectorsTransposed.txt.gz|This is the transposed version of the eigenvector file above.|
|CellTypeSpecificProbePCA.PCAoverSamplesPrincipalComponents.txt.gz|This file contains the actual principal components for the cell type specific probe PCA|
|SampleCorrelationMatrix.txt|The correlation matrix over the samples.|
|SamplePC1Correlations.txt|A PCA analysis has been performed using the sample correlation matrix. This file contains the correlations of all samples with the first principal component, as a quality control measure. Samples with a correlation < 0.9 with this first PC have been excluded to calculate the PCA that was used for CellTypeProxyFile.txt|
|ComparisonToCellCount.txt|If you have used the --cellcounts option, this file contains the correlation of PC1 with your actual cell count.|
|plot.pdf|If you have used the --cellcounts option, this file plots the cell-counts on the y-axis and the PC1 values on the x-axis.|

Additionally, the normalization step will have created a subdirectory in your `outputdir` containing the following files:

|File|Description|
|----|-----------|
|CellTypeSpecificProbeExpression.txt.gz|This file contains the (Quantile normalized, log2 transformed) gene expression data for the cell type specific probes.|
|ExpressionData-QNormLog2Transformed.CovariatesRemoved.txt.gz OR
ExpressionData-QNormLog2Transformed.txt.gz|This file contains the gene expression data you used as inexp, although now it is quantile normalized and log2 transformed. Please note that if you did not correct for MDS components, this file will not end with ‘CovariatesRemoved’.| 
|ExpressionDataPCQC-QNormLog2Transformed.CovariatesRemoved.txt.gz OR
ExpressionDataPCQC-QNormLog2Transformed.txt.gz|As described before, PC1 of the sample correlation matrix is used as a quality control measure. Samples with a correlation of < 0.9 with PC1 will be excluded from further analysis. This file contains only those samples that pass this threshold. This file is required for step 2.|
|PCAResults.PCAoverSamplesEigenvalues.txt.gz|PCA eigenvalues for the quality control PCA over the sample correlation matrix|
|PCAResults.PCAoverSamplesEigenvectors.txt.gz|PCA eigenvectors for the quality control PCA over the sample correlation matrix|
|PCAResults.PCAoverSamplesEigenvectorsTransposed.txt.gz|Transposed PCA eigenvectors for the quality control PCA over the sample correlation matrix|
|PCAResults.PCAoverSamplesPrincipalComponents.txt.gz|Principal components for the quality control PCA over the sample correlation matrix|

####Association analysis - Step 1 - Prepare your data
Locate and/or download the following files (to avoid confusion, **use full paths** when supplying these files to the software): 

-	The directory containing your TriTyper genotype data. We will call this directory `genotypedir`.
-	ExpressionDataPCQC-QNormLog2Transformed.CovariatesRemoved.txt.gz OR ExpressionDataPCQC-QNormLog2Transformed.txt.gz: the raw gene expression data as generated in the normalization step (described above). We will call this file `covariates`
-	ExpressionData.txt.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz: trait data generated by [Step 3 - Phenotype data normalization](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#step-3---phenotype-data-normalization) of the eQTL mapping pipeline, corrected for 40PCs (please make sure this file was not corrected for cis-eQTL effects but was corrected for MDS components). We will call this file `traitfile`
-	Find a location on your hard drive to store the output. We will refer to this directory as `outputdir`.
-	SNPProbe.txt: the file containing the *cis*-eQTL SNP probe combinations to test. One combination per line, with the first column showing the genotype variant identifier. Please use a file that is suitable for your platform. We will call this file `snpprobefile`.
-	GenotypeToExpressionCoupling.txt: the file linking genotypes to gene expression sample identifiers. This file is also optional (for example if your sample identifiers in your gene expression data correspond to the genotype sample identifiers). We will call this file `genotypetotraitcoupling`. The format is described here: [Genotype - phenotype coupling](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#genotype---phenotype-coupling).
-	CellCounts.txt: the file containing the proxy-phenotype generated by the normalization step described above. The format is that of a [Phenotype file, covariate file](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file). We will call this file `cellcountproxy`
-	Determine the number of available cores/processors in your machine (optional). We will call this `nrThreads`

####Association analysis - Step 2 - Run the association analysis

``java –jar eqtl-mapping-pipeline.jar --mode celltypespecific --step mapeqtls --inexp traitfile --covariates covariates --cellcounts cellcountproxy --in genotypedir --out outputdir --snpprobe snpprobefile --threads nrThreads``

**Please note that the --gte and --threads switches are optional**

####Association analysis - Step 3 - Check the output
After finalizing, the association analysis will have generated a couple of files in your `outputdir`:

|File|Description|
|----|-----------|
|CellTypeSpecificEQTLEffects.txt|This file contains the Z-scores for both the cis-eQTL as well as the interaction term with the (proxy) cell count| 
|CellTypeSpecificityMatrix.binary.columns.ser|This binary file contains the column description for the CellTypeSpecificityMatrix|
|CellTypeSpecificityMatrix.binary.dat|This binary matrix contains Z-scores: each column is a cis-eQTL, each row is a probe interaction Z-score, for each probe that was used as a covariate in the model. Also, this matrix contains two summary statistics (beta + Z-score) for the (proxy) cell count. This file can be used for eventual meta-analysis of cell-type specific effects. You can convert this file to a text-based matrix using the utility described in Step 3.|
|CellTypeSpecificityMatrix.binary.rows.ser|This binary file contains the row description of the CellTypeSpecificityMatrix|
|eQTLsNotPassingQC.txt|The set of eQTLs that did not pass QC|
|SNPSummaryStatistics.txt|This file contains summary statistics on the SNPs: MAF, HWE, Call-rate, Alleleles, minor allele, etc.|

####Association analysis - Step 4 - Convert the binary output table to text (optional)
If you want to investigate the data stored within the CellTypeSpecificityMatrix.binary.dat file, you can use the following command to convert the binary file to a text-based, tab-separated matrix:

``java –jar eqtl-mapping-pipeline.jar --mode util --convertbinarymatrix --in /path/to/ CellTypeSpecificityMatrix.binary.dat --out /path/to/textoutput.txt``

##Meta-analysis and settings file
The command line interface of this software allows for basic QTL analyses. However, our software has many more capabilities that are not accesible via the command line. In these cases, an XML file is required that describes the different settings (full path referred to as `settingsfile`). [An example `settingsfile` is provided in the repository](https://github.com/molgenis/systemsgenetics/blob/master/eqtl-mapping-pipeline/src/main/resources/settings.xml). Using a settings file allows you to quickly rerun certain analyses and to perform on-the-fly meta-analyses. A copy of the `settingsfile` will always be copied to your `outdir`. 

Currently, `settingsfile` can only be used in the `--mode metaqtl` mode. You should note that a `settingsfile` overrides all command line switches. The `settingsfile` can be used as follows:

```
java –jar eqtl-mapping-pipeline.jar --mode metaqtl --settings settingsfile
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
|-------|-----|-----------|
|sett.ing|double|description|
|sett.ing2|double|description2|

#File formats
This section lists the different file formats used by the software package.

##Trait file
This file is required by most parts of the program, since it contains the actual quantitative trait measurements. The format of this file is a simple text-based matrix, which may be gzipped. If you are exporting your data from R, make sure that the first line starts with an empty tab (otherwise you are bound to get column-shift problems).

###File example
<pre>
Probe   Sample1 Sample2	Sample3
0001    0.2     0.5     0.6
0002    0.8     0.6     0.6
0003    0.9    -1.5     7.9
</pre>


##Probe annotation file
A probe annotation file is required when running a *cis*-eQTL analysis. This file describes where each probe/trait/gene is located on the genome. This file is a simple tab-separated text file, with a line for each probe, and a header.

**Please note:** for reliable meta-analysis, probe annotations should be identical between datasets. Therefore, we recommend use of the array address for your platform as a probe identifier if you are using Illumina array based expression data.

###File example
<pre>
Platform    HT12v4-ArrayAddress Symbol	Chr		ChrStart    ChrEnd     Probe     Seq
HT12v4      00001               GeneX	1       1504        1554       0         CGCTCCCCTTATAACTT-etc.
HT12v4      00002               GeneY	11      19900       19950      1         GGATCCCAGATTCCCT-etc.
HT12v4      00003               GeneZ	23      101         151        2         TTCTCCAGAGTCGAGC-etc.
</pre>


##Phenotype file, covariate file
Phenotype and covariate files have the same basic format. We use a tab separated text-based table, with individuals on columns and probes/traits/genes or covariates on the rows. 

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

#General software information

##QTL File
QTL file description

##QTL mapping
MetaQTL is the part of the eQTLMappingPipeline that actually performs the QTL mapping. This part is actually used by multiple parts of the program, such as during the pcaoptimum method and the MixupMapper. The program is able to run genome-wide analyses, but is also able to run mapping on a selected number of SNPs, probes or a combination thereof. For a given SNP, the program first determines which probes it should be tested against. Generally for *cis*-eQTLs, this implies gene expression probes that are on the same chromosome and for which the genomic position of the middle of the probe is within 250kb of the SNP position. For *trans*-eQTLs, however, the minimum distance between the SNP and the gene expression probe should be at least 5Mb, or the gene expression probe should be located on a different chromosome. During eQTL mapping, we use a Hardy-Weinberg p-value threshold of 0.0001, a minor allele frequency threshold of 0.05 and a call-rate threshold of 95% as a quality control on the SNPs we test.

Once a list of probes has been assigned to a SNP that passes quality control, we calculate a number of statistics. For each combination of SNP and gene expression probe, we calculate the correlation with the (imputed) genotype with the ranked (normalized) gene expression data for the gene expression probe using Spearman’s correlation. From this calculation, we can derive a t-statistic, using the following formula, where r is the correlation and n is the number of samples for which there is genotype data:


t = r / sqrt(1-r<sup>2</sup>) / n - 2

From the t-statistic, we can calculate the z-score, from the cumulative t-distribution. Finally, using an inverse normal distribution, we can determine the p-value describing the significance of the association. However, for the meta-analysis, we don’t use the p-values from the individual cohorts, but instead, we store the z-scores.

##Meta-analysis
Our software can perform meta-analysis using a weighted Z-score method, described by Whitlock et al: for a given SNP j and gene expression probe k from dataset i and the number of individuals n for which genotype data was available, the QTL Z-score is calculated as follows:

Z<sub>weighted<sub>Dataset<sub>i</sub>, SNP<sub>j</sub>, Probe<sub>k</sub></sub> 
= sqrt( n <sub>dataset<sub>i</sub>, SNP<sub>j</sub></sub>) x Z<sub>Dataset<sub>i</sub>, SNP<sub>j</sub>, Probe<sub>k</sub></sub>

The meta-analysis Z-score over all datasets is calculated as follows:

Z<sub>sum<sub>SNP<sub>j</sub>, Probe<sub>k</sub></sub> = sum( Z<sub>weighted<sub>Dataset<sub>i</sub>, SNP<sub>j</sub>, Probe<sub>k</sub></sub></sub> )

Finally, the summed Z-score is weighted for the total sample size ( N ):

Z<sub>meta<sub>SNP<sub>j</sub>, Probe<sub>k</sub></sub> = Z<sub>sum<sub>SNP<sub>j</sub>, Probe<sub>k</sub></sub></sub> / sqrt( N )

Since Zmeta follows a normal distribution, we can calculate the eventual meta p-value from the normal cumulative distribution.

##Multiple testing correction
When performing a statistical test, often a threshold is set at a p-value of 0.05 to declare significance of a test. During QTL mapping, sometimes millions of tests can be performed. However, depending on the data, when randomizing the data you will observe that many p-values will actually be below this threshold. As a consequence, often the threshold is divided by the number of performed tests (Bonferroni correction), which yields a new threshold for significance. However, because many of the genotypes and probes/traits/genes may be correlated, applying the Bonferroni correction may be too stringent (this multiple testing burden is especialy present when performing *trans*-analyses).

Because of the correlation structure in both the genotype and phenotype data, our software applies a permutation strategy: the software runs the analysis on the real data, and afterwards repeats the same analysis on the data, but this time shuffles the links between genotype and phenotype. This means we do not touch the genotype or phenotype data itself, which keeps the correlations between SNPs and between probes intact. From both the permuted data as well as the real data, we can create two distributions, which can then be used to determine the  False Discovery Rate (FDR). 

In other words: by randomly shuffling the sample labels, we can compute the sampling distribution for any test statistic, under the strong null hypothesis that a set of genetic variants has absolutely no effect on the outcome. FDR controls the expected proportion of false positives among the results based on the real data analysis.  

##MixupMapper
To identify sample mix-ups, MixupMapper uses each significantly detected *cis*-eQTL in the dataset. For each of these *cis*-eQTLs, the mean (μAA, μAB and μBB) and standard deviation (σAA, σAB and σBB) of the gene expression values were determined for each of the three genotypes (AA, AB and BB). For each pair of genotype and gene expression array we determined the SNP genotype (g), and calculated the number of standard deviations that the gene expression value (e) differed from the expected value associated with the SNP genotype using an absolute Z-score. For each sample pair the absolute Z-scores of all significant cis-eQTLs were summed and the average Z-score for each sample pair was determined to account for differences in the number of tested eQTLs per sample pair due to missing SNP genotypes.

MixupMapper normalizes the Z-scores by subtracting the average of the overall Z-scores for the expression sample and divides it by the standard deviation of the overall Z-scores for this expression sample. Similarly, the Z-scores were normalized by subtracting the average of the overall Z-scores for the genotype sample and were divided by the standard deviation of the overall Z-scores for the genotype sample. After these normalizations MixupMapper determines what the expression array was with the lowest overall normalized Z-score for each genotyped sample. This expression sample was considered to reflect the particular genotyped sample. Once the best matching expression sample had been identified for each genotyped sample, it will be compared to what had been initially defined. 

