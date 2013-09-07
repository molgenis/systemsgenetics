#!/bin/bash

#Script to run the genotype aligner


genotypeAlignerDir=`dirname $0`

java -Xmx1g -jar ${genotypeAlignerDir}/GenotypeAligner.jar $*