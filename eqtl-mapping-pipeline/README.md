EQTL mapping pipeline manual
============================

This software allows for (e)QTL mapping using linear models and direct meta-analysis of such data.

Downloading the software
------------------------
You can download the latest version of the software here: [Latest version](http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$eqtl-mapping-pipeline/lastBuild/).
Make sure to download the stand alone jar: `eqtl-mapping-pipeline-******-jar-with-dependencies.jar`

Please note that the manual refers to eQTLMappingPipeline.jar, while the name of the package described above may be different (because of different version numbers etc)

Documentation Change Log
----------------------------------------------
Changes in the documentation can be found here.

Before you start
================

Java Virtual Machine
--------------------
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

Step by step eQTL analysis
==========================
This is a step by step guide which will guide you through the QTL mapping process using our software. Please note that this step by step guide illustrates a  Our method consists of six steps described below:

<ol>
<li>Step 1 - Preparation phenotype data</li>
<li>Step 2 - Preparation genotype data</li>
<li>Step 3 - Phenotype data normalization</li>
<li>Step 4 - Mapping mix-ups</li>
<li>Step 5 - Determining the optimum number of PCs to remove</li>
<li>Step 6 - Perform the final QTL analysis</li>
</ol>

Step 1 - Preparation phenotype data
==========================
Because our software uses a nonparametric test by default, you can use virtually any continuous data as trait values to map a variety of QTL effects. However, currently the normalization tools provided with this package are focused on array based methylation data, array based (Illumina) expression data, and preprocessed RNA-seq data (e.g. transcript level quantified data) or (GC)-RMA processed Affymetrix data. 
Please format your phenotype data in a simple tab separated text file. The format of this file is described here: [Data Formats - Phenotype data](https://github.com/harmjanwestra/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file).

**Throughout the manual, we will refer to your phenotype file as `traitfile`**

Illumina array based expression data
------------------------------------
Because Illumina has several ways to annotate and preprocess their gene expression arrays, we provide a detailed instruction below, that will yield raw, untouched gene expression levels.

1.	Use the GenomeStudio files of your expression arrays. **Please do not use any normalization, transformation, imputation, or background correction method that is offered by Illumina’s GenomeStudio! Just load your raw data, without doing any normalization, transformation, or background correction.** 
	-	Export a so-called ‘Final Report’ using GenomeStudio, including array_address_id (as the probe identifier) and the average expression signals per sample for all probes (click 'Analysis' in the top menu bar, then 'Reports').
	-	If you create more than one Final Report, merge the Final Reports, so that all data will be combined in one Expression Matrix File
	-	Rows should contain the different probes, and columns should contain the different sample IDs
	-	Remove any header information that Genome studio might produce: the header of the matrix is expected at the first row
	-	Create a new header (expected at the first row and column) for both rows and columns, describing the sample names (columns) and the probe IDs (rows). Note that you use the probe array address for your platform as the probe name.
    -   The final file should look like this format: [Data Formats - Phenotype data](https://github.com/harmjanwestra/systemsgenetics/tree/master/eqtl-mapping-pipeline#phenotype-file-covariate-file).


Step 2 - Preparation of genotype data
==========================
Our software is able to use both unimputed called genotypes, as well as imputed genotypes and their dosage values. However, currently our software can only interpret files that are in the [TriTyper](http://genenetwork.nl/wordpress/trityper/) format. We are currently working on a generic method of genotype input through [Genotype IO](https://github.com/harmjanwestra/systemsgenetics/tree/master/Genotype-IO). In the mean time, users of the eQTL mapping pipeline should convert their data to TriTyper using ImputationTool. This tool is also integrated in the eQTL mapping pipeline, and can be called using the following command:

````
java –d64 –Xmx4g –jar eQTLMappingPipeline.jar 
````


Step 3 - Phenotype data normalization
==========================
Step 4 - Mapping mix-ups
==========================
Step 5 - Determining the optimum number of PCs to remove
==========================
Step 6 - Perform the final QTL analysis
==========================


File formats
============

Probe annotation file
----------------------

Phenotype file, covariate file
------------------------------