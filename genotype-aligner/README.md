Genotype Aligner
================

The Genotype Aligner is an easy to use commandline tool that allows harmonization of genotype data 
stored using different fileformats with different and potentially unknown strands. 

LD patterns are used to determine the correct strand GC and AT SNPs and by using 
the [Molgenis Genotype Reader](https://github.com/PatrickDeelen/systemsgenetics/tree/master/molgenis-genotype-reader) we can import and export different file format.

Getting started
----------------

### Downloading the Genotype Aligner
The last build from the genotype aligner can be downloaded here:
http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$genotype-aligner/lastBuild/

Click on `genotype-aligner-*.*.*-dist.zip` or `genotype-aligner-*.*.*-dist.tar.gz` to download

In case of succesfull build there will a green circel before the `Build #`. 
It is possible that you visit the website when a new build is in progress, please try again in a few minutes.

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
	--inputType PLINK_BED \
	--ref /data/demoRefData \
	--refType VCF \
	--output /data/demoOuput \
	--outputType PLINK_BED \
```

Note: this is a single commandline command. The `\` is only for readabily.

In case of this example the programm exprects that `/data/demoInputData.bed`, `/data/demoInputData.bim`  and `/data/demoInputData.fam` exist and that ` /data/demoRefData.vcf.gz` and `/data/demoRefData.vcf.gz.tbi` exist.

`/data/demoOuput.bed`, `/data/demoOuput.bim`, `/data/demoOuput.fam` and `/data/demoOuput.log` will be created.

**Note: it is important to make sure that both study and refernece are using the same genome build**

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

### Using SHAPEIT2 output

The output format of SHAPEIT2 is documented on their website: http://www.shapeit.fr/pages/m02_formats/hapssample.html. However the actual output does not contain the chromosome on this column. `--forceChr` option forces the input data chromosome to the specified value. Note this is only valid if all variants are present on this chromosome. This feature is currently only implemented for the input data and not the reference data. Feel free to raise a [new issue](https://github.com/molgenis/systemsgenetics/issues/new) on our github project to request this.

Typical usage scenarios
----------------

### Preparing data for genotype imputation

When imputing genotype data the strand of both the study data to impute and the reference data used for imputation need to be identical. Some imputation tools can swap the strand of non-ambigous SNPs but this is not possible for AT and GC SNPs. AT and GC can be swapped using minor allele frequency but this is not reliable, especially for variants with a high minor allele frequency. The Genotype Aligner solves these problems by using LD structure of nearby variants. 

In combination with the `--update-id` option the Genotype Aligner is a convineant preperation of genotype data before imputation.

Aligning prephased data, generated using tools like SHAPEIT2, is particulary usefull. A dataset only needs to be prephased once and can then be aligned and imputed using different reference set or different versions of reference sets.

### Merging data from different genotyping platforms

When combining different genotype datasets, either samples ran on multiple genotyping chips or different batches of samples, it is important to have identical strands. Merge data in plink will give a warning when it detects strand issues in non AT or non GC SNPs but can not automaticly correct this. AT and GC SNP swaps are not automaticly detected. 

The `--keep` option is particulary usefull here to keep the SNPs not shared by both datasets. The `--update-id` will also make merging using plink or other tools more easy.

Arguments overview
----------------

| Short | Long           | Description                                                                                                                          |
|-------|----------------|--------------------------------------------------------------------------------------------------------------------------------------|
| -i    | --input        | The base path of the data to align. The extensions are determined based on the input data type.|
| -I    | --inputType    | The input data type. (see inputeType options) |
| -r    | --ref          | The base path of the reference data used for alignment. The extensions are determined based on the input data type.|
| -R    | --refType      | The input data type. (see refType options) |
| -o    | --output       | The base path of the output data. |
| -O    | --outputType   | The output data type. (--outputType options) |
| -id   | --update-id    | Update the variant identifiers using the reference data. The identifiers of the output data will be the same as the reference data |
| -l    | --min-ld       | The minimum LD (r2) between the variant to align and potential supporting variants |
| -m    | --min-variants | The minimum number of supporting variant before before we can do an alignment |
| -s    | --variants     | Number of flanking variants to consider |
| -f    | --forceChr     | SHAPEIT2 does not output the sequence name in the first column of the haplotype file. Use this option to force the chromosome for all variants. This option is only valid incombination with `--inputType SHAPEIT2`
| -c    | --check-ld     | Also check the LD structure of non AT and non GC variants. Variants that do not pass the check are excluded. |
| -d    | --debug        | Activate debug mode. This will result in a more verbose log file |

####--inputeType /--refType options

Base path refers to either --input or --ref 

* PED_MAP
 * Expects plink PED file at: `${base path}.ped`
 * Expects plink MAP file at: `${base path}.map`
 * Note: it is recommened to use PLINK_BED due to the large memory usage of the current implementation
* VCF
 * Expects VCF file at: `${base path}.vcf.gz`
 * Must be compresed usign bgzip. (see chapter: Preparing a VCF file)
 * Expects tabix file at: `${base path}.vcf.gz.tbi` (see chapter: Preparing a VCF file)
* PLINK_BED
 * Expects plink BED file at: `${base path}.bed`
 * Expects plink BIM file at: `${base path}.bim`
 * Expects plink FAM file at: `${base path}.fam`
 * Must be in SNP major mode. This is the default of  plink 1.07.
* VCFFOLDER
 * Matches all vcf.gz files in the folder specified with the bash path.
* SHAPEIT2
 * Expects haps file at: `${base path}.haps`
 * Expects sample file at: `${base path}.sample`
 * See also chapter on SHAPEIT2 output.

#####--outputType options

Base path refers to --output

Regardless of the output type a log file will always be created at: `${base path}.log`

* PED_MAP
 * Writes plink PED file to: `${base path}.ped`
 * Writes plink MAP file to: `${base path}.map`
 * Note: it is recommened to use PLINK_BED since the writing is much faster
* PLINK_BED
 * Writes plink BED file to: `${base path}.bed`
 * Writes plink BIM file to: `${base path}.bim`
 * Writes plink FAM file to: `${base path}.fam`
 * Data is written in SNP major mode
* SHAPEIT2
 * Writers haps file at: `${base path}.haps`
 * Writers sample file at: `${base path}.sample`
 
Bug reports / feature requests
----------------

Please use the GitHub issue tracker to post feature request or to report bugs. https://github.com/molgenis/systemsgenetics/issues/new.

Test data
----------------
This chapter is not relevant for the usage of the program but allows reproducibility of the test data

The genotype aligener contains test data. For the genotype data to align we use HapMap3 data and as a reference we use 1000G data. 

This dataset is always tested when building the project and by our Jenkins server (http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$genotype-aligner/)

### HapMap3 data

The following tools are needed for this script:

* plink
* ucsc liftover + chain hg18ToHg19
* SHAPEIT2
* Genetic Map (b37)

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/latest_phaseIII_ncbi_b36/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/latest_phaseIII_ncbi_b36/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/latest_phaseIII_ncbi_b36/plink_format/relationships_w_pops_121708.txt

tar -zxvf hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
tar -zxvf hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2

#Create list of CEU sampels to extract
awk '$7 == "CEU" {print $1,$2}' relationships_w_pops_041510.txt > ceuSamples.txt

#Extract first 6Mb of chr20 for CEU samples
plink --noweb --chr 20 --file ../hapmap3_r3_b36_fwd.consensus.qc.poly --out hapmap3CeuChr20B36Mb6 --from-mb 0 --to-mb 6 --recode --keep ceuSamples.txt

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

We also want this data phased using SHAPEIT2, for aditional testing.

```Bash
shapeit.v2.r644.linux.x86_64 --input-bed ./hapmap3CeuChr20B37Mb6RandomStrand -M genetic_map_chr20_combined_b37.txt --output-max ./hapmap3CeuChr20B37Mb6RandomStrand --noped --thread 4
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
bgzip 1000gCeuChr20Mb6.recode.vcf > 1000gCeuChr20Mb6.vcf.gz
#Create index using tabix
tabix -p vcf 1000gCeuChr20Mb6.vcf.gz
```
