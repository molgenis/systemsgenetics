EQTL mapping pipeline
========================

This software allows for (e)QTL mapping using linear models and direct meta-analysis of such data.

Downloading the software
------------------------
You can download the latest version of the software here:
http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$eqtl-mapping-pipeline/lastBuild/
Make sure to download the stand alone jar: `eqtl-mapping-pipeline-******-jar-with-dependencies.jar`


Before you start
================

Java Virtual Machine
--------------------
Our software is written in Java, which makes the software both fast and portable across multiple operating systems. Executing a java program is similar to executing a normal program or app, although there are some considerations:
*	eQTL mapping heavily relies on available memory. Make sure your machine is 64-bit and has lots of memory installed (at least 4Gb). 
*	Please make sure your version of java is up-to-date. Use at least the 64-bit version of Java version 6 (also called version 1.6). You can download the so-called Java Runtime Environment (JRE) from: http://www.java.com. If you are running Linux, the virtual machine may be present in the proprietary section of your package manager, or may be available via http://www.java.com
*	Java executables are called jar-files. The eQTL mapping pipeline is such a jar-file. You can execute it by using the following command from a terminal / console:
```	java –jar eQTLMappingPipeline.jar ```
*	You need to specify the amount of memory available to the program using the command-line switch –Xmx. The amount of memory can be specified in megabytes (using an m suffix) or in gigabytes (using a g suffix). To be sure your computer is running java at 64-bit, please add the switch –d64. Both switches (–Xmx and –d64) should be called prior to the –jar switch. An example where the program is allowed to use 4gb of memory:
```	java –d64 –Xmx4g –jar eQTLMappingPipeline.jar  ```
*	Try to increase the –Xmx amount when you get Out-Of-Memory-Errors or errors involving ‘heap space’

**IMPORTANT NOTE: In this manual, we assume you understand the principle that you need to allocate sufficient amounts of RAM and therefore we excluded the –Xmx switch from the example commands. Please be aware that you should use it, as most of the commands require a substantial amount of memory!**
**The eQTL mapping pipeline is a command line program, which makes the user interface not very intuitive. In order to help you a bit, an overview of available switch options is displayed when a command is incomplete or incorrect. Furthermore, each mode of the program also has its own overview of available switches with a small description of its functionality. For example: “java –jar eQTLMappingPipeline.jar” produces a list of available modes, while “java –jar eQTLMappingPipeline.jar  --mode metaqtl” produces a list of all available options for metaqtl. 
You can also find all possible switches in the manual section “command line options”.**


