Molgenis Genotype Reader
========================

The molgenis genotype reader is a java library that allows accessing genotype data form difference sources in a uniform and fast manner.

Downloading latested jar
------------------------

Most user can simply download the latested jar and use it in there programs.

http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$molgenis-genotype-reader/lastBuild/

Make sure to download the stand alone jar: `molgenis-genotype-reader-******-jar-with-dependencies`


Features
------------------------
### Readers

* Plink PED/MAP
* Plink binary format BED/BIM/FAM
* VCF
* TriTyper

### Writers

* Plink PED/MAP
* Plink binary format BED/BIM/FAM
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
bgzip -c example.vcf > example.vcf.gz
tabix -p vcf example.vcf.gz
```

Examples
------------------------

### Basic genotype data

The normal genotype data interface allows the iteration over variants

````java
String datasetPath = "/some/path/file"; //Note omition of extentions

//Create a binary plink genotype data object
GenotypeData genotypeData = null;
try {
	genotypeData = new BedBimFamGenotypeData(datasetPath, 1000);
} catch (IOException ex) {
	LOGGER.fatal("IO error: " + ex.getMessage());
	System.exit(1);
}

for(GeneticVariant variant : genotypeData){
	//Iterate over all variants
	System.out.println(variant.getPrimaryVariantId());
}

for(Sample sample : genotypeData.getSamples()){
	//Note sex is an enum, to string outputs a human readable form of the gender
	System.out.println("Sample ID: " + sample.getId() + " sex: " + sample.getSex());
}
````

###Random access genotype data
In most cases the random access genotype reader is more usefull. 3 examples are shown on how they can be loaded.

````java
RandomAccessGenotypeData randomAccessGenotypeData = null;

try {

	//Here we read the dataset of the specified type using the auto loader
	randomAccessGenotypeData = new BedBimFamGenotypeData(datasetPath);

	//Or
	randomAccessGenotypeData = RandomAccessGenotypeDataReaderFormats.PLINK_BED.createGenotypeData(datasetPath, 1000);

	//Or
	randomAccessGenotypeData = RandomAccessGenotypeDataReaderFormats.valueOf("PLINK_BED").createGenotypeData(datasetPath, 1000);

} catch (IOException ex) {
	LOGGER.fatal("IO error: " + ex.getMessage());
	System.exit(1);
} catch (GenotypeDataException ex) {
	LOGGER.fatal("Genotype data error: " + ex.getMessage());
	System.exit(1);
}

for(String sequenceNames : randomAccessGenotypeData.getSeqNames()){
	//No a sequence can be chr or contig
	System.out.println("Seq/Chr: " + sequenceNames);
}

for(GeneticVariant variant : randomAccessGenotypeData.getVariantsByRange("1", 1, 10)){
	//get variants from chr 1 pos 1 to 10 (exclusive)
	System.out.println("Variant ID: " + variant.getPrimaryVariantId());
}

//Create hashmap based on primary ID. 
//Note this takes lots of memory. 
//Variants without an ID will always be ignored.
//If multiple variants have the same ID an arbritary variant is selected.
HashMap<String, GeneticVariant> variantHashMap = randomAccessGenotypeData.getVariantIdMap();
````

###Genetic variant

Some examples on what can be done with genetic variants. 

````java
GeneticVariant snp = randomAccessGenotypeData.getSnpVariantByPos("1", 1);
//Get snp variant at this position. null if not present

snp.isSnp(); //must be true since we dit getSnpVariantByPos

for(Alleles sampleAlleles : snp.getSampleVariants()){
	//Iterate over the alleles form all samples.
	System.out.println(sampleAlleles.getAllelesAsString());
}

for(byte dosage : snp.getSampleCalledDosages()){
	//Iterate over the dosage values of snps. (0,1,2)
	System.out.println(dosage);
}

for(float dosage : snp.getSampleDosages()){
	//Iterate over the dosage values of snps. range 0 - 2
	System.out.println(dosage);
}

snp.getMinorAllele();
snp.getMinorAlleleFrequency();

snp.getCallRate();

snp.getHwePvalue();

snp.isBiallelic();

//chr != 0 && pos != 0
snp.isMapped();

snp.isAtOrGcSnp();
````

###Alleles and Allele

The Alleles object stores a collection of Allele objects

````java
Alleles snpAlleles =  snp.getVariantAlleles();
		
for(Allele a : snpAlleles){
	//Print the alles found for this variant
	System.out.println("Allele: " + a);
}

Alleles alleles1 = Alleles.createAlleles(Allele.A, Allele.C);
Alleles alleles2 = Alleles.createAlleles(Allele.C, Allele.A);
Alleles alleles3 = Alleles.createAlleles(Allele.A, Allele.C);

if(alleles1 == alleles2){}; // false because order is different
if(alleles1 == alleles3){}; // true because it is guaranteed that if alleles and order is endentical it is the same object
if(alleles1.sameAlleles(alleles2)){}; // true
````

###Linkage Disequilibrium

We allow easy LD calculation

````java
Ld ld = null;
try {
	ld = snp.calculateLd(snp2);
} catch (LdCalculatorException ex) {
	LOGGER.fatal("Error in LD calculation: " + ex.getMessage(), ex);
	System.exit(1);
}
ld.getR2();
ld.getDPrime();
````

###Filtering random access genotype data

Filtering a dataset on variants or samples will create a new RandomAccessGenotypeData that when accessed will only reviel the included variants and samples as if the others do not even exist. Data is not copied, instead the filtered datasets are backed by the orignal

####Filtering variants

````java
HashSet<String> includedSnps = new HashSet<String>();
includedSnps.add("rs1");
includedSnps.add("rs2");

//Create filter
VariantFilter variantFilter = new VariantIdIncludeFilter(includedSnps);

//Filter to only include selected variants
RandomAccessGenotypeData genotypeData1 = new VariantFilterableGenotypeDataDecorator(randomAccessGenotypeData, variantFilter);
````

####Filtering samples

````java
HashSet<String> includedSsamples = new HashSet<String>();
includedSsamples.add("pietje");
includedSsamples.add("jansje");

SampleFilter sampleFilter = new SampleIdIncludeFilter(includedSsamples);

//Filter the data with the selected variants to only include the selected samples
RandomAccessGenotypeData randomAccessGenotypeDataSampleFilter = new SampleFilterableGenotypeDataDecorator(randomAccessGenotypeDataVariantFilter1, sampleFilter);

//Here we define a combined variant filter consting of 2 filters. This can contain as many filters as requered
VariantFilter combinedFilter = new VariantCombinedFilter(new VariantQcChecker(0.05f, 0.95f, 0.001d), new VariantFilterBiAllelic());

//Now we also do some QC of the variants after the samples are filtered. maf 0.05, call rate 0.95 and hwe-p 0.001
RandomAccessGenotypeData randomAccessGenotypeDataVariantFilter2 = new VariantFilterableGenotypeDataDecorator(randomAccessGenotypeDataSampleFilter, combinedFilter);
````

####Autoloader and filters
It is also possible to use the auto loader from the basicUsage examples in combination with filters. This is faster and more memory efficient for some filetypes, like TriTyper and in the futere (or now if this is not updated) also for ped_map and binary plink. I is posible to use this method also on the other files types altough it does not give increased performance. It is possible to set either the variant filter or the sample filter to `null` include all variants or samples.

Note that samples are filtered first and then the variant filter is applied, so a when filtering on the minor allele frequency then it will be caculated only using the included samples.

```java
RandomAccessGenotypeDataReaderFormats.VCF.createFilteredGenotypeData(datasetPath, 1000, combinedFilter, sampleFilter);
```

###Using Genotype-IO in rJava
Reading in genotype information in R can be a big problem due to file size and fileformat support. Using rJava and the Genotype-IO API this is no problem anymore. Using the example method below you can easaly read in any supported fileformat. Just put in the basepath of the input and fileformat. If necessary one can change the cachesize and include filters for samples and variants.


```S
loadGenotypeData <- function( basePath, dataType, cacheSize=1000, variantFilter = .jnull(class = "org/molgenis/genotype/variantFilter/VariantFilter"), sampleFilter = .jnull("org/molgenis/genotype/sampleFilter/SampleFilter")){
  dataType <- toupper(dataType)
  genotypeDataFormat <- .jcall("org/molgenis/genotype/RandomAccessGenotypeDataReaderFormats", "Lorg/molgenis/genotype/RandomAccessGenotypeDataReaderFormats;","valueOf", dataType)
  return(.jcall(genotypeDataFormat, "Lorg/molgenis/genotype/RandomAccessGenotypeData;", "createFilteredGenotypeData", basePath, as.integer(cacheSize), variantFilter, sampleFilter))
}
```

In the code below examples for the variant filter and sample filter are given. "includedSnps" and "includedSamples" are R character arrays.

```S
variantFilter <-  .jcast(.jnew("org/molgenis/genotype/variantFilter/VariantIdIncludeFilter",includedSnps),"org/molgenis/genotype/variantFilter/VariantFilter")
sampleFilter <- .jcast(.jnew("org/molgenis/genotype/sampleFilter/SampleIdIncludeFilter",includedSamples), "org/molgenis/genotype/sampleFilter/SampleFilter")
```

A short example of a use case of this code is for instance when trying to read in genotype data, select a specific variant and printing the histogram of dosages.

```S
genotypeData <- loadGenotypeData("PathToFiles", "Plink_BED")
snp <- .jcall(genotypeData, "Lorg/molgenis/genotype/variant/GeneticVariant;", "getSnpVariantByPos", "8", as.integer(18257854))
hist(as.numeric(.jcall(snp, "[B", "getSampleCalledDosages")))
```

Use the rJava vignet for more information on rJava.


###More examples
The `org.molgenis.genotype.examples` package contains these and other basic examples.


Working with the code
------------------------
* Install Maven (http://maven.apache.org/) if you haven't done so already
* Build molgenis-genotype-reader (`mvn clean install`)

### How to open in Eclipse

* Install the m2e Maven plugin; This is standard installed in the more recent versions of Eclipse, if you got an older version you probably need to install it.

* Import the project in Eclipse with `File/Import -> Existing Maven Projects` and select the molgenis-genotype-reader root folder



