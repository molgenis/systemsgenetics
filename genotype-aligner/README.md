Genotype Aligner
================

The Genotype Aligner is an easy to use commandline tool that allows harmonization of genotype data 
stored using different fileformats with different and potentially unknown strands. 

LD patterns are used to determine the correct strand GC and AT SNPs and by using 
the Molgenis Genotype Reader we can import and export differnt file format.

Getting started
----------------

### Downloading the Genotype Aligner
The last build from the genotype aligner can be downloaded here:
http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$genotype-aligner/lastBuild/

Click on `genotype-aligner-***-jar-with-dependencies.jar` to download

In case of succesfull build there will a blue circel before the `Build #`. 
It is possible that you visit the website when a new build is in progress, please try again in a few minutes or feel free to download an older build http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$genotype-aligner/changes

### Running the Genotype Aligner
type `java -jar genotype-aligner-***-jar-with-dependencies.jar` to run. You will now get an overview of the different commandline options

### Basic usage
In the most basic usage scenario you need to define:

* A dataset that you want to align and the type of this dataset
* A dataset that you want to use as refernce and the type of this dataset
* The output path and type where you want to write the aliged data to

You command will look like this:
```
Java -jar genotype-aligner-***-jar-with-dependencies.jar \
	--input /data/demoInputData \
	--inputType PED_MAP \
	--ref /data/demoRefData \
	--refType VCF \
	--output /data/demoOuput \
	--outputType PED_MAP \
```

Note: this is a single commandline command. The `\` is only for readabily.

In case of this example the programm exprects that `/data/demoInputData.map` and `/data/demoInputData.ped` exist and that ` /data/demoRefData.vcf.gz` and `/data/demoRefData.vcf.gz.tbi` exist.

`/data/demoOuput.map`, `/data/demoOuput.ped` and `/data/demoOuput.log` will be created.

### Advanced options
TODO

### Using VCF files
Before VCF files can be used they need to be compressed using bgzip and indexed with a tabix. This prevents having to read all data into memory yet still allows quick access.

#### Installing tabix and bgzip

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

### Shapeit2 output
TODO

Arguments overview
----------------

TODO

Test data
----------------
This chapter is not relevant for the usage of the program but allows reproducibility of the test data

The genotype aligener contains test data. For the genotype data to align we use HapMap3 data and as a reference we use 1000G data. 

This dataset is always tested when building the project and by our Jenkins server (http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$genotype-aligner/)

### HapMap3 data

The following tools are needed for this script:

* plink
* ucsc liftover + chain hg18ToHg19

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/latest_phaseIII_ncbi_b36/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/latest_phaseIII_ncbi_b36/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/latest_phaseIII_ncbi_b36/plink_format/relationships_w_pops_121708.txt

tar -zxvf hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
tar -zxvf hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2

#Create list of CEU sampels to extract
awk '$7 == "CEU" {print $1,$2}' relationships_w_pops_041510.txt > ceuSamples.txt

#Extract first 6Mb of chr20 for CEU samples
plink --noweb --chr 20 --file ../hapmap3_r3_b36_fwd.consensus.qc.poly --out hapmap3CeuChr20B36Mb6 --from-mb 0 --to-mb 6  --recode --keep ceuSamples.txt

# - liftover to b37 - 

#Create bed 
awk '{$5=$2;$2=$4;$3=$4+1;$1="chr"$1;print $1,$2,$3,$5}' OFS="\t" hapmap3CeuChr20B36Mb6.map > hapmap3CeuChr20B36Mb6b36.bed

#Update mapping
liftOver -bedPlus=4 hapmap3CeuChr20B36Mb6b36.bed hg18ToHg19.over.chain hapmap3CeuChr20B36Mb6b37.bed hapmap3CeuChr20B36Mb6unmapped.txt

#All snps are mapped. Normally we would have to account for this

#Create mapping update list used by plink
awk '{print $4, $2}' OFS="\t" hapmap3CeuChr20B36Mb6b37.bed > hapmap3CeuChr20B36Mb6b37.txt

#Update plink mappings
plink --noweb --file hapmap3CeuChr20B36Mb6 --recode --out hapmap3CeuChr20B37Mb6 --update-map hapmap3CeuChr20B36Mb6b37.txt

#No we have to again create a plink file to make sure the implied order is correct after liftover.
plink --noweb --file hapmap3CeuChr20B37Mb6 --out hapmap3CeuChr20B37Mb6 --make-bed
```

We have now created a subset of the hapmap3 which is all in forward strand. We are no going to swap a large number variants. The aliger can identify these swapped variants and flip them back to forward strand using the 1000G data.

```Bash
#Create swap list 50% of SNPs
awk '
  BEGIN { srand(1)} 
  { if (rand() <= .5) print $2}
' < hapmap3CeuChr20B37Mb6.bim > flipList.txt

plink --noweb --bfile hapmap3CeuChr20B37Mb6 --make-bed --flip flipList.txt --out hapmap3CeuChr20B37Mb6RandomStrand

```

### 1000G

The following tools are needed for this script:

* vcftools
* tabix

```bash
wget ftp://share.sph.umich.edu/1000genomes/fullProject/2012.03.14/phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz.tgz
tar xvzf phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz.tgz

#Create subset of data
vcftools --gzvcf chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz --out 1000gCeuChr20Mb6 --chr 20 --from-bp 0 --to-bp 6000000 --recode --remove-indels --remove-filtered-all
#Comprese subset using bgzip (part of tabix package)
gzip 1000gCeuChr20Mb6.recode.vcf > 1000gCeuChr20Mb6.vcf.gz
#Create index using tabix
tabix -p vcf 1000gCeuChr20Mb6.vcf.gz
```
