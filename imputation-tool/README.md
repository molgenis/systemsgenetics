#Imputation-tool manual
This tool provides several ways of processing (imputed) TriTyper genotype data. Additionally this program provides converter programs to convert your datasets into [TriTyper](link) format. If you have any questions or remarks regarding this software or the manual, please contact me at westra.harmjan@outlook.com

##Downloading the software
You can download the latest version of the software here: [Latest version](http://www.molgenis.org/jenkins/job/systemsgenetics/lastStableBuild/nl.systemsgenetics$imputation-tool/).
Make sure to download the stand alone jar: `imputation-tool-******-jar-with-dependencies.jar`

The genotype harmonization capabilities of this software have been replaced by the [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/tree/master/Genotype-Harmonizer), and will therefore not be documented in this manual (the focus will be on converting your genotype data to TriTyper format). Genotype Harmonizer is capable of processing most genotype formats natively, through the implementation of the [Genotype IO](https://github.com/molgenis/systemsgenetics/tree/master/Genotype-IO) libraries.

Please note that the manual refers to ImputationTool.jar, while the name of the package described above may be different (because of different version numbers etc).

##Path definitions and commands
Please note that our software expects full paths, although shorter relative paths wil also work in most cases. So, if you are on Windows, a full path to a genotype directory would be similar to `c:\path\to\genotype\dir\` and a full path to a file would be `c:\path\to\genotype\directory\file.txt`. Linux and Mac OS use different path separators. On these systems, these paths would be similar to the following `/path/to/genotype/dir/` and `/path/to/genotype/dir/file.txt`. Our main point here is that when pointing to a directory, use a 'trailing slash'

This manual will combine references to paths with commands that need to be issued for a certain task. For example, at some point in this manual we refer to your genotype data as `genotypedir`, which will be printed in a grey box. Commands will also be in grey boxes, and can make references to paths defined earlier (to keep the manual readable), as follows:

````
java -jar ImputationTool.jar --mode convert --in genotypedir
````

##Java Virtual Machine
Our software is written in Java, which makes the software both fast and portable across multiple operating systems. Executing a java program is similar to executing a normal program or app, although there are some considerations:

* Please make sure your version of java is up-to-date. Use at least the 64-bit version of Java version 6 (also called version 1.6). The software should also work on version 7 and higher. If you are running on Windows, you can download the so-called Java Runtime Environment (JRE) from: http://www.java.com. If you are running Linux, the virtual machine may be present in the proprietary section of your package manager, or may be available via http://www.java.com.

* Java executables are called jar-files. The eQTL mapping pipeline is such a jar-file. You can execute it by using the following command from a terminal / console:
```    
    java –jar ImputationTool.jar 
```

* You need to specify the amount of memory available to the program using the command-line switch –Xmx (case-sensitive!). The amount of memory can be specified in megabytes (using an m suffix) or in gigabytes (using a g suffix). To be sure your computer is running java at 64-bit, please add the switch –d64. Both switches (–Xmx and –d64) should be called prior to the –jar switch. An example where the program is allowed to use 4gb of memory:
```    
    java –d64 –Xmx4g –jar ImputationTool.jar 
```

* Try to increase the –Xmx amount when you get Out-Of-Memory-Errors or errors involving ‘heap space’

**IMPORTANT NOTE:** In this manual, we assume you understand the principle that you need to allocate sufficient amounts of RAM and therefore we excluded the –Xmx switch from the example commands. Please be aware that you should use it, as some of the commands may require a substantial amount of memory (depending on your dataset and settings)!


##General information about software
* The ImputationTool is a command line program, which makes the user interface not very intuitive. In order to help you a bit, an overview of available switch options is displayed when a command is incomplete or incorrect. Furthermore, each mode of the program also has its own overview of available switches with a small description of its functionality. For example: ```java –jar ImputationTool.jar``` produces a list of available modes.
* This software is also part of our [QTL Mapping pipeline](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline). You can call the imputationtool from within the QTL Mapping pipeline as follows:

```
java –jar eQTLMappingPipeline.jar --imputationtool
```

#Converting Illumina FinalReport files to TriTyper
After genotype calling by Illumina's Genome Studio, so-called FinalReportFiles can be generated. Genotype calls in these FinalReportFiles can be exported using several formats. We have found that the TOP format produces the most consistent results across GenomeStudio releases. Note down the full path to your FinalReportFile `finalreportfile`, and a directory to store the converted data `trityperoutputdir`. 

Then, issue the following command:
```
java –jar ImputationTool.jar --mode ftt --in finalreportfile --out trityperoutputdir
```



#Converting PLiNK PED and MAP files to TriTyper
ImputationTool accepts the standard PLINK text-based [PED and MAP files](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml) as input for conversion to TriTyper. We currently do not have a method of converting Binary PED files (which have the .bed extension). For information on interconverting PED, MAP, and BED files, please refer to the PLINK website. For some projects, you may have your genotype data split over multiple files. If this is the case for your dataset, ImputationTool can convert these files all at once. *Do make sure however, that the map files have exactly the same file name as the ped files (i.e. having only a different extension). Also make sure you have a MAP file for all PED files. Finally, make sure that if your PED files are GZipped, your MAP files are as well.*

Place your PED and MAP file(s) in a folder, the full path of which we will refer to as `pedmapdir`. Then, find a spot on your harddrive to store the converted TriTyper output, the full path of which we will refer to as `trityperoutputdir`. Then, issue the following command:

```
java –jar ImputationTool.jar --mode pmtt --in pedmapdir --out trityperoutputdir
```

#Converting MACH imputed files to TriTyper
ImputationTool can also convert data imputed with [MACH](http://www.sph.umich.edu/csg/abecasis/MACH/tour/imputation.html). The program will look for files with the following extensions: .mlgeno, .mlinfo, .mlgeno and .mldose. Currently, this tool does not support GZipped files. As with the PED and MAP data, please make sure that if you have multiple files, you use the same file name for each matching set of files (i.e. only change the extension). 

Place your MACH output file(s) in a folder, the full path of which we will refer to as `machdir`. Then, find a spot on your harddrive to store the converted TriTyper output, the full path of which we will refer to as `trityperoutputdir`. Then, issue the following command:

```
java –jar ImputationTool.jar --mode mtt --in machdir --out trityperoutputdir
```

#Converting Minimac imputed files to TriTyper
ImputationTool can also convert data imputed with [Minimac](http://genome.sph.umich.edu/wiki/Minimac). The program will look for files with the following extensions: .info and .dose. Currently, this tool does not support GZipped files. As with the PED and MAP data, please make sure that if you have multiple files, you use the same file name for each matching set of files (i.e. only change the extension). 

Place your Minimac output file(s) in a folder, the full path of which we will refer to as `minimacdir`. Then, find a spot on your harddrive to store the converted TriTyper output, the full path of which we will refer to as `trityperoutputdir`. Then, issue the following command:

```
java –jar ImputationTool.jar --mode mmtt --in minimacdir --out trityperoutputdir
```

#Converting BEAGLE imputed files to TriTyper
[Beagle](http://faculty.washington.edu/browning/beagle/beagle.html) is also a frequently used imputation program. ImputationTool can also convert files created by this program, which includes files ending with .gprobs and .r2. We have not tested this tool with the neweste Beagle version (Beagle v4). ImputationTool expects your .gprobs files to be GZipped, and your r2 files to be unzipped (which is the default Beagle output format). The description below assumes procedures as implemented in our imputation pipeline. If you want to use this tool, you probably have to convert your filenames to a matching scheme.

Start by placing all your Beagle output files in a single directory. We'll refer to this directory as `beagledir`. Also, find a spot on your harddrive to store the converted TriTyper dataset `trityperoutputdir`. 

Beagle imputed data is often divided in several batches of samples. Count the number of batches in your dataset `nrbatches` (should be an integer). ImputationTool can convert all these batches at once, although it does require you to define your file naming scheme. In our imputation pipeline, we used batches of 300 samples. Batches were named using the following scheme: aa, ab, ac etc. Our imputation pipeline then ran imputation per batch of samples, per chromosome. As a consequence, the Beagle output files would have a filename similar to this:

```
Chr-1.aa-Chr-1.aa.gprobs.gz
Chr-1.ab-Chr-1.ab.gprobs.gz
Chr-1.ac-Chr-1.ac.gprobs.gz
Chr-2.aa-Chr-2.aa.gprobs.gz
Chr-2.ab-Chr-2.ab.gprobs.gz
Chr-2.ac-Chr-2.ac.gprobs.gz
```

From this we then created a `filenametemplate`, by replacing the chromosome number with the text CHROMOSOME and by replacing the batch identifier with BATCH. In our example above, the `filenametemplate` would be:

`Chr-CHROMOSOME.BATCH-Chr-CHROMOSOME.BATCH`

The above imputation output can then be converted using the following command:

```
java –jar ImputationTool.jar --mode bttb --in beagledir --tpl filenametemplate --out trityperoutputdir --size nrbatches
``` 

If you want to convert a single chromosome (`chromosomenr`, should be integer), append the above command with the following command line switch

```
--chr chromosomenr
```

If you want to convert a subset of chromosomes, append the above command using the following command line switch

```
--chrstart startchromosomenr --chrend endchromosomenr
```


#Converting IMPUTE imputed files to TriTyper
Files generated after imputation using the [Impute program](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html) can be converted to TriTyper using the ImputationTool. Imputation using Impute is often performed in batches of a couple megabases. Our software allows for the conversion of these batches all at once. Data may be GZipped. Just put all your imputed files in a folder `imputedir`. However, due to the format of Impute output, you will also have to supply the number of samples `nrsamples` (should be an integer) that was imputed using Impute. Finally, locate a spot on your harddrive to store the converted TriTyper dataset `trityperoutputdir`. *Make sure your `imputedir` only contains Impute output files.*

Then, issue the following command:

```
java –jar ImputationTool.jar --mode itt --in imputedir --out trityperoutputdir --nrSamples nrsamples
```

The Impute format does not store sample identifiers. If you want to use the original sample identifiers in your TriTyper dataset, create a file containing one sample per line, and note down the full path to this file. We will call the full path to this file `samplefile`. Then append the command above using the following command line switch:


```
--samples samplefile
```

Impute output files can be huge. If you want to convert only a selection of samples, create a text file with one sample identifier per line, and note down the full path to this file (`samplestoincludefile`). Then append the command above using the following command line switch:

```
--samplestoinclude samplestoincludefile
```

Additionally, if you want to convert a subset of SNPs, create a text file with one SNP identifier per line, and note down the full path to this file (`snpstoincludefile`). Then append the command above using the following command line switch:

```
--snps snpstoincludefile
```

By default, the program looks for files containg either chr#, chr_# or chr-# (where # is the chromosome number). If your dataset follows a different numbering scheme, you can use a regular expression to specify another patterh `regexpattern` by appending the following command line switch: 

```
--fileMatchRegex regexpattern
```



#Converting VCF files to TriTyper
1000 genomes data is often stored in the [VCF format](http://www.1000genomes.org/node/101). These files can also be converted to Trityper, and may be GZipped. Simply put your VCF files in a folder `vcfdir`, then, find a spot on your harddrive to store the converted TriTyper output, the full path of which we will refer to as `trityperoutputdir`. Then, issue the following command:

```
java –jar ImputationTool.jar --mode mmtt --in vcfdir --out trityperoutputdir
```

By default, the program will consider all lines in the VCF files as being biallelic genotypes. However, sometimes a VCF file may also contain other genotypic variants. In such case, the type of variant is declared in the INFO field. In such a case, you can specify how your file defines SNPs using a regular expression pattern. For example, if your SNPs are defined using the following string: `BiAllelicSNP`, you can tell the program to look specifically for these variants using the following command:

```
java –jar ImputationTool.jar --mode mmtt --in vcfdir --out trityperoutputdir --pattern BiAllelicSNP
```

