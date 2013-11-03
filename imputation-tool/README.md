#Imputation-tool manual
This tool provides several ways of processing (imputed) TriTyper genotype data. Additionally this program provides converter programs to convert your datasets into [TriTyper](link) format.

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
* The software is able to process GZipped text files for most of the input files (not files in .tar archives however), which allows you to save some space on your hard drive.

#Converting PLiNK PED and MAP files to TriTyper

#Converting MACH imputed files to TriTyper

#Converting BEAGLE imputed files to TriTyper

#Converting IMPUTE imputed files to TriTyper

