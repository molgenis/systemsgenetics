<!--

  It is recommended to view this file using a markdown viewer
  or view this readme online: github.com/Molgenis/systemsgenetics/blob/master/GeneticRiskScoreCalculator/README.md

-->
Genetic Risk Score Calculator
================

Readme using GitHub markdown: 

* https://help.github.com/articles/markdown-basics
* https://help.github.com/articles/github-flavored-markdown

Editor: http://jbt.github.io/markdown-editor/

Getting started
----------------

Arguments overview
----------------


| Short | Long           | Description                                 |
|-------|----------------|---------------------------------------------|
| -i    | --input        | \* The base path of the data to align. The extensions are determined based on the input data type.|
| -I    | --inputType    | The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path (see inputeType options) |
| -o    | --output       | \* The base path of the output data. |

\* = required

####--inputeType options

Base path refers to --input value

* PED_MAP
 * Expects Plink PED file at: `${base path}.ped`
 * Expects Plink MAP file at: `${base path}.map`
 * Note: it is strongly recommend to use PLINK_BED due to the large memory usage of the current implementation. When dealing with large datasets it is recommend to first use plink to convert the data binary plink using this command `plink --file pedmapData --make-bed --out binaryData`
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
 * See also chapter on 'Using SHAPEIT2 Output and Oxford Gen format'
* GEN
 * Expects gen file at: `${base path}.gen` or without the extentions `${base path}`
 * Expects sample file at: `${base path}.sample`
 * See also chapter on 'Using SHAPEIT2 Output and Oxford Gen format'
* TRITYPER
  * Expects folder with TriTyper data. Can handle dosage matrix

 
Bug reports / feature requests
----------------

Please use the GitHub issue tracker to post feature request or to report bugs. https://github.com/molgenis/systemsgenetics/issues/new.

