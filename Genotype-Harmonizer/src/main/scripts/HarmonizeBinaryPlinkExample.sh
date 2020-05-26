#!/bin/bash

GenotypeHarmonizerDir=`dirname $0`

cd ${GenotypeHarmonizerDir}

sh GenotypeHarmonizer.sh \
	--inputType PLINK_BED \
	--input ./exampleData/hapmap3CeuChr20B37Mb6RandomStrand \
	--update-id \
	--outputType PLINK_BED \
	--output ./exampleOutput/binaryPlinkExampleOut \
	--refType VCF \
	--ref ./exampleData/1000gCeuChr20Mb6


echo ""
echo "Demo output saved in: ${GenotypeHarmonizerDir}/exampleOutput/"
