# MbQTL - QTL analysis with beta-distributed null
This repository contains a version of MetaQTL (https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline)
that uses a different multiple testing strategy than the original MetaQTL that was used previously used in the eQTLgen and sc-eQTLgen studies.
The main differences between MetaQTL and this version are that this version no longer makes use of the TriTyper genotype format, but rather uses
tabix indexed VCF files as input, and that this version implements a different multiple testing correction strategy. 
While the previous version of MetaQTL determined a null distribution over the top associated variants of all genes, the method
implemented here is more akin to FastQTL (http://fastqtl.sourceforge.net/) or QTLTools (https://qtltools.github.io/qtltools/).
The main difference between FastQTL and this QTL software, is that while FastQTL assesses associations over all included samples,
this software allows on the fly meta-analysis of different datasets. Additionally, while FastQTL generally uses Pearson correlations to
determine relationships between phenotype and genotype, this software uses Spearman correlations by default.

# Compilation
This software depends on the 'genetica-libaries' repository (https://github.com/molgenis/systemsgenetics/tree/master/genetica-libraries) for various functions.
Download the full https://github.com/molgenis/systemsgenetics repository, and import as a maven project in e.g. IntelliJ Idea. Then, build the genetica-libraries module first,
and then this module.

# Download
A binary version of this software is available here: </br>
https://jenkins.harmjanwestra.nl/job/systemsgenetics_hjw/lastBuild/nl.systemsgenetics$MbQTL/

# Usage

The full list options can be viewed by executing the application without any parameters:

```console
[cluster-user@cluster software]$ java -jar MbQTL-1.4.2-SNAPSHOT-jar-with-dependencies.jar 
Command line parse exception: Missing required option: m
usage: mbqtl.jar
 -a,--annotation <PATH>         Gene annotation file
    --chr <INT>                 Chromosome number
    --ciswindow <INT>           Cis window size [default: 1mb]
    --cr <FLOAT>                Call-rate threshold [default: 0.95]
 -e,--exp <PATH>                Expression matrix (can be gzipped)
 -eg,--expgroups <PATH>         File defining groups of phenotypes
    --eqtlfile <FILE>           [determineldgwas] - eQTL file - txt.gz
                                file
    --eqtlset <FILE>            [determineldgwas] - List of eQTL snps to
                                test - txt.gz file
 -g,--gte <PATH>                Genotype to expression to dataset linkfile
                                (tab separated)
 -gl,--genelimit <PATH>         Gene limit file (one line per gene ID)
    --gwasassocfile <FILE>      [determineldgwas] - GWAS association file
                                - txt.gz file
    --gwaslistfile <FILE>       [determineldgwas] - GWAS study file -
                                txt.gz file
    --gwasset <FILE>            [determineldgwas] - List of GWAS snps to
                                test - txt.gz file
    --hwep <FLOAT>              Hardy-Weinberg p-value threshold [default:
                                0.0001]
    --input <arg>               Input file
    --input2 <arg>              Input file 2
    --ldthreshold <FLOAT>       [determineldgwas] - LD threshold [default:
                                0.8]
 -m,--mode <STRING>             Mode:
                                [metaqtl|mbqtl|mbqtlsingleds|mbqtlplot|reg
                                ressqtl|sortfile|determineld|determineldgw
                                as|concatconditional]
    --maf <FLOAT>               Minor allele frequency threshold [default:
                                0.01]
    --matchbyrsid               [determineldgwas] - Match by RsId in stead
                                of full variant id
    --minobservations <FLOAT>   Require at least this many observations
                                per dataset (i.e. non-NaN
                                genotypes/phenotypes) [default: 10]
    --norank                    Do not rank expression data
    --nrdatasets <INT>          Minimum number of datasets required in
                                meta-analysis [default: 2]
    --nriters <arg>             Nr iters to concatenate using
                                concatconditional
 -o,--out <PATH>                Output prefix
    --outputall                 Output all associations, not just top
                                association per gene
    --outputallpermutations     Output all permuted associations, not just
                                top association per gene
    --perm <INT>                Number of permutations [default: 1000]
    --prunedistance <INT>       [determineldgwas] - Pruning distance
                                [default: 1mb]
    --prunethreshold <FLOAT>    [determineldgwas] - Pruning LD threshold
                                [default: 0.2]
    --replacemissinggenotypes   Replace missing genotypes with average
                                genotype: use this when both genotypes and
                                expression data have missing values and
                                perm > 0
    --seed <LONG>               Random seed [default: 123456789]
 -sgl,--snpgenelimit <PATH>     SNP-gene limit file (one line per snp-gene
                                combination, tab separated)
    --skipchr6                  [determineldgwas] - Skip chr6 when
                                calculating LD
 -sl,--snplimit <PATH>          SNP limit file (one line per gene ID)
    --snpannotation <FILE>      SNP annotation file, tab separated: SNPID
                                chr pos
    --snplog                    Output SNP summary stats per snp/gene
                                pair.
    --sortbyz                   Sort by Z-score
    --testnonparseablechr       [mbqtl] - Test variants and genes that map
                                to non-oarseable chromosomes (e.g. patch
                                chromosomes)
 -v,--vcf <PATH>                Tabix indexed VCF
```

A number of parameters will be explained below

### -a,--annotation

This is a tab-separated file containing the annotations for the features included in the expression file

| Platform | ArrayAddress | Symbol | Chr | ChrStart | ChrEnd | Probe | Strand
| :------- | :------------------------------ | :------------------------------ | - | :------ | :------ | :------------------------------ | -
| 10x_v3.1 | ENSG00000078369_ENSG00000198804 | ENSG00000078369_ENSG00000198804 | 1 | 1785285 | 1891117 | ENSG00000078369_ENSG00000198804 | -
| 10x_v3.1 | ENSG00000078369_ENSG00000198938 | ENSG00000078369_ENSG00000198938 | 1 | 1785285 | 1891117 | ENSG00000078369_ENSG00000198938 | -
| 10x_v3.1 | ENSG00000078369_ENSG00000198712 | ENSG00000078369_ENSG00000198712 | 1 | 1785285 | 1891117 | ENSG00000078369_ENSG00000198712 | -
| 10x_v3.1 | ENSG00000078369_ENSG00000198899 | ENSG00000078369_ENSG00000198899 | 1 | 1785285 | 1891117 | ENSG00000078369_ENSG00000198899 | -
| 10x_v3.1 | ENSG00000078369_ENSG00000251562 | ENSG00000078369_ENSG00000251562 | 1 | 1785285 | 1891117 | ENSG00000078369_ENSG00000251562 | -
| 10x_v3.1 | ENSG00000078369_ENSG00000198727 | ENSG00000078369_ENSG00000198727 | 1 | 1785285 | 1891117 | ENSG00000078369_ENSG00000198727 | -
| 10x_v3.1 | ENSG00000078369_ENSG00000198840 | ENSG00000078369_ENSG00000198840 | 1 | 1785285 | 1891117 | ENSG00000078369_ENSG00000198840 | -
| 10x_v3.1 | ENSG00000078369_ENSG00000198886 | ENSG00000078369_ENSG00000198886 | 1 | 1785285 | 1891117 | ENSG00000078369_ENSG00000198886 | -
| 10x_v3.1 | ENSG00000078369_ENSG00000211592 | ENSG00000078369_ENSG00000211592 | 1 | 1785285 | 1891117 | ENSG00000078369_ENSG00000211592 | -


### -e,--exp

This is a tab-separated file containing the expression of each feature per sample. The first column contains the features names. The rest of the columns are the sample names

| phenotype_id | SAMPLE0101 | SAMPLE0102 | SAMPLE0103 | SAMPLE0104
| :-------------------------------| :------------------ | :------------------- | :----------------- | :------------------
| ENSG00000078369_ENSG00000198804 | -0.0522097510600244 | -0.0295754428440798  | 0.202974140896468  | -0.0276880518949276
| ENSG00000078369_ENSG00000198938 | 0.0450745639059172  | -0.0672694684512958  | 0.0184521946269516 | 0.022512806829165
| ENSG00000078369_ENSG00000198712 | -0.0141565967202825 | -0.00490039746665645 | 0                  | -0.0406239420986158
| ENSG00000078369_ENSG00000198899 | 0.165412275032376   | 0.0580198320270886   | 0.276782919404274  | 0.0528713582666007
| ENSG00000078369_ENSG00000251562 | 0.344754320610721   | 0.32330837105614     | 0.538005119742025  | 0.0261570703361926
| ENSG00000078369_ENSG00000198727 | 0.0386521803204386  | 0.0996988465401328   | 0.138989056460888  | 0.0080459526571041
| ENSG00000078369_ENSG00000198840 | -0.0160774862150101 | 0.0617262806145658   | 0.171037928759085  | -0.0923119406889147
| ENSG00000078369_ENSG00000198886 | 0.173666906863794   | 0.139841698948845    | 0.0144995229332674 | -0.0299424428882636
| ENSG00000078369_ENSG00000211592 | -0.441021383760673  | -0.124575974108114   | -0.203135231509838 | 0.0179638736664903


### -eg,--expgroups

This is a tab-separated file containing the expression groups. During the permutation step, where the minimum permuted P-value is calculated, the minimum is taken per group. This parameter is optional, and if not supplied, the groups are each feature. The first column is the feature, and the second column is the feature group. There is no header.

ENSG00000078369_ENSG00000198804 | ENSG00000078369
| :-----------------------------| :--------------
ENSG00000078369_ENSG00000198938 | ENSG00000078369
ENSG00000078369_ENSG00000198712 | ENSG00000078369
ENSG00000078369_ENSG00000198899 | ENSG00000078369
ENSG00000078369_ENSG00000251562 | ENSG00000078369
ENSG00000078369_ENSG00000198727 | ENSG00000078369
ENSG00000078369_ENSG00000198840 | ENSG00000078369
ENSG00000078369_ENSG00000198886 | ENSG00000078369
ENSG00000078369_ENSG00000211592 | ENSG00000078369
ENSG00000078369_ENSG00000198888 | ENSG00000078369


### -g,--gte

This is a tab-separated file mapping the genotype identifiers to expression identifiers, and the dataset the samples are from. The first column is the genotype ID, the second the expression ID, and the third column the dataset (name). There is no header

SAMPLE0101 | SAMPLE0101 | Dataset1
| :--------| :----------| :-------
SAMPLE0102 | SAMPLE0102 | Dataset1
SAMPLE0103 | SAMPLE0103 | Dataset1
SAMPLE0104 | SAMPLE0104 | Dataset1
SAMPLE0105 | SAMPLE0105 | Dataset1


### -v,--vcf

This is the genotype file to be used. Format should be (a g-zipped) VCF file. This file needs to be tabix-indexed.


### --perm

This is the number of permutations that is performed. The default is set to 1000, though in most cases a number of ten would already suffice.


### -m,--mode

Main mode to run the application in. In most cases this would be QTL-mapping, which is the mbQTL mode.
