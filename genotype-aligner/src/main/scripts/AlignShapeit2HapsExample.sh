#!/bin/bash

GenotypeAlignDir=`dirname $0`

cd ${GenotypeAlignDir}

sh GenotypeAligner.sh \
	--inputType SHAPEIT2 \
	--input ./exampleData/hapmap3CeuChr20B37Mb6RandomStrand \
	--update-id \
	--outputType SHAPEIT2 \
	--output ./exampleOutput./shapeit2ExampleOut \
	--refType VCF \
	--ref ./exampleData/1000gCeuChr20Mb6 \
	--forceChr 20


echo ""
echo "Demo output saved in: ${GenotypeAlignDir}/exampleOutput"
