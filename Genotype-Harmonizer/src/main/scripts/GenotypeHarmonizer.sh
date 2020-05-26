#!/bin/bash

#Script to run the genotype aligner

genotypeHarmonizerDir=`dirname $0`

java -Xmx1g -jar ${genotypeHarmonizerDir}/GenotypeHarmonizer.jar $*