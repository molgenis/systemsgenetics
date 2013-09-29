<!--

  It is recommended to view this file using a markdown viewer
  or view this readme online: https://github.com/PatrickDeelen/systemsgenetics/edit/master/genotype-harmonizer/README.md

-->
Genotype Harmonizer
================

The Genotype Harmonizer is an easy to use command-line tool that allows harmonization of genotype data 
stored using different file formats with different and potentially unknown strands. 

Linkage disequilibrium (LD) patterns are used to determine the correct strand GC and AT SNPs and by using 
the [Genotype IO](https://github.com/PatrickDeelen/systemsgenetics/tree/master/Genotype-IO) package we can import and export different file format.

Getting started
----------------

### Downloading the Genotype Harmonizer
The last build from the genotype harmonizer can be downloaded here:
http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$Genotype-Harmonizer/lastBuild/

Click on `genotype-harmonizer-*.*.*-dist.zip` or `genotype-harmonizer-*.*.*-dist.tar.gz` to download the Genotype Harmonizer, test data and 2 example scripts.

In case of successfully build there will a green circle before the `Build #`. 
It is possible that you visit the website when a new build is in progress (the circle will be blinking), please try again in a few minutes.

### Running the Genotype Harmonizer
type `GenotypeHarmonizer.sh`, `GenotypeHarmonizer.bat` or `java -jar GenotypeHarmonizer.jar` to run. You will now get an overview of the different command-line options.

In the case of an heapspace or out of memory error you need allocate more memory to run the Genotype Harmonizer. If this should happen use this command to run: `Java -jar -Xms##g -Xmx##g -jar GenotypeHarmonizer.jar`. Replace ## with the number of gigabytes of memory you want to allocate.

### Basic usage
In the most basic usage scenario you need to define:

* A dataset that you want to align and the type of this dataset
* A dataset that you want to use as reference and the type of this dataset
* The output path and type where you want to write the aligned data to

Your command will look like this:
```
GenotypeHarmonizer.sh \
	--input path_to_study_data \
	--output path_of_output \
	--ref path_to_reference
```

Note: this is a single command line command. The `\` is only for readability.

You can find more examples as a script in the root of the distribution.

**Note: it is important to make sure that both study and reference are using the same genome build**

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
bgzip -c example.vcf > example.vcf.gz
tabix -p vcf example.vcf.gz
```

### Using SHAPEIT2 output

The output format of SHAPEIT2 is documented on their website: http://www.shapeit.fr/pages/m02_formats/hapssample.html. However, the actual output does not contain the chromosome on the first column. The `--forceChr` option forces the input data to the chromosome which is specified. Note this is only valid if all variants are indeed on this chromosome. This feature is currently only implemented for the input data and not for the reference data. Feel free to raise a [new issue](https://github.com/molgenis/systemsgenetics/issues/new) on our github project to request this.

Typical usage scenarios
----------------

### Preparing data for genotype imputation

When imputing genotype data the strand of both the study data to impute and the reference data used for imputation need to be identical. Some imputation tools can swap the strand of non-ambiguous SNPs but this is not possible for AT and GC SNPs. AT and GC can be swapped using minor allele frequency but this is not reliable, especially for variants with a high minor allele frequency. The Genotype Harmonizer solves these problems by using LD structure of nearby variants. 

In combination with the `--update-id` option the Genotype Harmonizer is a convenient tool for preparation of genotype data before imputation.

Aligning pre-phased data, generated using tools like SHAPEIT2, is particularly useful. A dataset only needs to be pre-phased once and can then be aligned and imputed using different reference set or different versions of reference sets.

### Merging data from different genotyping platforms

When combining different genotype datasets, either samples ran on multiple genotyping chips or different batches of samples, it is important to have identical strands. Merge data in Plink will give a warning when it detects strand issues in non AT or non GC SNPs but can not automatically correct this. AT and GC SNP swaps are not automatically detected. 

The `--keep` option is particularly useful here to keep the SNPs not shared by both datasets. The `--update-id` will also make merging using Plink or other tools more easy.

### Conversion

Although there are many tools that can convert from one file type to an other it can also easily be done using our Genotype Harmonizer by simply not specifying a reference. In such a case the input data will simply be converted to the specified output type.

Arguments overview
----------------

| Short | Long           | Description                                 |
|-------|----------------|---------------------------------------------|
| -i    | --input        | \* The base path of the data to align. The extensions are determined based on the input data type.|
| -I    | --inputType    | The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path (see inputeType options) |
| -r    | --ref          | The base path of the reference data used for alignment. The extensions are determined based on the input data type. If not specified the input data is simply converted to the specified output type.|
| -R    | --refType      | The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path (see refType options) |
| -o    | --output       | \* The base path of the output data. |
| -O    | --outputType   | The output data type. Defaults to --inputType or to PLINK_BED if there is no writer for the impute type. (--outputType options) |
| -id   | --update-id    | Update the variant identifiers using the reference data. The identifiers of the output data will be the same as the reference data |
| -l    | --min-ld       | The minimum LD (r2) between the variant to align and potential supporting variants. Defaults to 0.3 |
| -m    | --min-variants | The minimum number of supporting variant before before we can do an alignment. Defaults to 3 |
| -v    | --variants     | Number of flanking variants to consider. Defaults to 100 |
| -f    | --forceChr     | SHAPEIT2 does not output the sequence name in the first column of the haplotype file. Use this option to force the chromosome for all variants. This option is only valid in combination with `--inputType SHAPEIT2`
| -c    | --check-ld     | Also check the LD structure of non AT and non GC variants. Variants that do not pass the check are excluded. |
| -d    | --debug        | Activate debug mode. This will result in a more verbose log file |
| -ma   | --mafAlign     | If there are not enough variants in LD and the minor allele frequency (MAF) of a variant <= the specified value in both study as in reference then the minor allele can be used as a backup for alignment. Defaults to 0 |

\* = required

####--inputeType /--refType options

Base path refers to either --input or --ref 

* PED_MAP
 * Expects Plink PED file at: `${base path}.ped`
 * Expects Plink MAP file at: `${base path}.map`
 * Note: it is recommend to use PLINK_BED due to the large memory usage of the current implementation
* VCF
 * Expects VCF file at: `${base path}.vcf.gz`
 * Must be compressed using bgzip. (see chapter: Preparing a VCF file)
 * Expects tabix file at: `${base path}.vcf.gz.tbi` (see chapter: Preparing a VCF file)
* PLINK_BED
 * Expects Plink BED file at: `${base path}.bed`
 * Expects Plink BIM file at: `${base path}.bim`
 * Expects Plink FAM file at: `${base path}.fam`
 * Must be in SNP major mode. This is the default of Plink 1.07.
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
 * Writes Plink PED file to: `${base path}.ped`
 * Writes Plink MAP file to: `${base path}.map`
 * Note: it is recommend to use PLINK_BED since the writing is much faster
* PLINK_BED
 * Writes Plink BED file to: `${base path}.bed`
 * Writes Plink BIM file to: `${base path}.bim`
 * Writes Plink FAM file to: `${base path}.fam`
 * Data is written in SNP major mode
* SHAPEIT2
 * Writers haps file at: `${base path}.haps`
 * Writers sample file at: `${base path}.sample`
 
####Tweaking the alignment using the advanced options

It can be worthwhile to tweak the alignment algorithm using the advanced options (`--min-ld`, `--min-variants`, `--variants` & `--mafAlign`) to improve the reliability of the alignment and to reduce the number of excluded variant that could not be aligned. The default values are quite conservative and sooner exclude variant than falsely swap the strand.

We have found that for datasets with a small number of samples it is wiser to have a high value for `--min-ld` (as is the default). This resulted in a reliable aliment of the small example datasets. However for larger datasets a lower value can be used to retain more variants, we have good experiences with a cut-off 0.1 for datasets with more than 1000 samples. Increasing `--min-variants` safely allows lower values of `--min-ld` this can improve the number of aligned variants.

When imputing the HLA region in densely genotyped dataset containing a large number of rare variants it helped to increase the `--variants` although this did effect the speed performance it did allow us to align more variants. This resulted in better imputation when using [snp2hla](https://www.broadinstitute.org/mpg/snp2hla/) from the Broad Institute.

Finally, it is also possible to align using the minor allele when the LD alignment fails using the `--mafAlign` option. By default this option is off (the maximum maf is set to 0) so it is never done. However in a lot of cases it is save to use the minor when there are not enough in LD with a variant. A save value would be 0.1, in that case only variant with a maf <= to 0.1 will be aligned by matching the minor allele.
 
Bug reports / feature requests
----------------

Please use the GitHub issue tracker to post feature request or to report bugs. https://github.com/molgenis/systemsgenetics/issues/new.

Test data
----------------
This chapter is not relevant for the usage of the program but allows reproducibility of the test data.

The genotype harmonizer contains test data. For the genotype data to align we use HapMap3 data and as a reference we use 1000G data. 

This dataset is always tested when building the project and by our Jenkins server (http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$genotype-harmonizer/). It is also supplied in the Genotype Harmonizer package to get you started.

### HapMap3 data

The following tools are needed for this script:

* Plink
* UCSC liftover + chain hg18ToHg19
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

#All SNPs are mapped. Normally we would have to account for this

#Create mapping update list used by Plink
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

We also want this data phased using SHAPEIT2, for additional testing.

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
