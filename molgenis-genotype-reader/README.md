Molgenis Genotype Reader
========================

The molgenis genotype reader is a java library that allows accessing genotype data form difference sources in a uniform and fast way.

Downloading latested jar
------------------------

Most user can simply download the latested jar and use it in there programs.

http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$molgenis-genotype-reader/lastBuild/

Make sure to download the stand alone jar: `molgenis-genotype-reader-******-stand-alone.jar`


Features
------------------------
### Readers

* Plink PED/MAP
* Pink binary format
* VCF
* TriTyper

### Writers

* Plink PED/MAP
* IMPUTE2 phased haplotypes

### Memory efficient

With the execption of Plink PED/MAP files only indices are stored in memory. Sample genotypes are only loaded when accessed. The loading of genotypes is cached preventing unnecessary disk IO.

### Sample and Variant filters

Easy methods to select a subset of samples or variants based on identifiers or features.

### Variants

* SNP's and indel's
* Basic statistics; MAF, HWE, callrate
* LD calculation between variants
* Both dosage and called genotype support

### Editing data

In memory support to modify data. Writers can be used to save modifications.

### High unit testing coverage

Usage notes
------------------------

### VCF

Only VCF files that are compressed using bgzip and indexed with a tabix are supported. This prevents having to read all data into memory.

#### Installing tabix

```bash
wget http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
tar -jxvf tabix-0.2.6.tar.bz2
cd tabix-0.2.6/
make
```

#### Preparing a VCF file

```bash
bgzip example.vcf > example.vcf.gz
tabix -p vcf example.vcf.gz
```

Examples
------------------------

The `org.molgenis.genotype.examples` package contains some basic examples. This should get you started


Working with the code
------------------------
* Install Maven (http://maven.apache.org/) if you haven't done so already
* Follow instructions in the global systemsgenetics readme.md on how to instal local jars
* Build molgenis-genotype-reader (`mvn clean install`)

### How to open in Eclipse

* Install the m2e Maven plugin; This is standard installed in the more recent versions of Eclipse, if you got an older version you probably need to install it.

* Import the project in Eclipse with `File/Import -> Existing Maven Projects` and select the molgenis-genotype-reader root folder



