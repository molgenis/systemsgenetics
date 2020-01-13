#!/bin/bash

# Usage
# ./Depict2Run.sh <traitfile> <output folder> <parameterfile>
# Path to the depict2 compatible zscore matrix
# Parameter file is a bashscript containing the paramters voor een run

# If any exit code in this script is not 0, quit
set -e
# Load java
ml Java

START_TIME=$(date +%s)

# User defined settings for this run
INPUT="$1"
OUTPUT="$2"
PARAMS="$3"
PHENOTYPE="$(basename $INPUT)"

# Load the parameter file
source $PARAMS

# DEPICT2 settings
TYPE="RUN"
DEPICT2="${BUNDLE_DIR}/depict2/Depict2-2.0.${VERSION}-SNAPSHOT/Depict2.jar"

if [ -f "${DEPICT2}" ]; then
	echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] Using version: ${VERSION}"
else 
	echo "[ERROR - $(date '+%Y-%m-%d %H:%M:%S')] Depict2 jar not found at: ${DEPICT2}"
	exit 1;
fi

if [ "${USE_LOCAL_TMP}" == true ]; then
	# TMPDIR on local storage for speedy IO
	# Copy input zscores to tmpdir
	mkdir -p ${TMPDIR}/${PHENOTYPE}_${VERSION}
	cp ${INPUT}* ${TMPDIR}/${PHENOTYPE}_${VERSION}
	# Write output to local storage
	TMP_IN=${TMPDIR}/${PHENOTYPE}_${VERSION}/${PHENOTYPE}
	TMP_OUT="${TMPDIR}/${PHENOTYPE}_${VERSION}"
	echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] copied input to tmpdir"
else 
	TMP_IN=${INPUT}
	TMP_OUT=${OUTPUT}/${PHENOTYPE}_${VERSION}
fi

mkdir -p ${TMP_OUT}

# Construct depict command
CMD="java -Xmx${MEMORY} -XX:ParallelGCThreads=${CORES} -jar ${DEPICT2} \
-g ${TMP_IN} \
-m ${TYPE} \
-o ${TMP_OUT}/${PHENOTYPE} \
-ge ${REFERENCE_ENSEMBL} \
-r ${REFERENCE_GENOTYPES} \
-rs ${REFERENCE_SAMPLES} \
-R ${REFERENCE_TYPE} \
-t ${CORES} \
-v ${VARIANT_PRUNING} \
-w ${GENE_WINDOW} \
-gpr ${GENE_PRUNING} \
-p ${PERMUTATIONS} \
-pr ${PERMUTATION_RESCUE} \
-pgc ${PERMUTATION_GENECOR} \
-ppe ${PERMUTATION_PATHWAY} \
-gcw ${GENE_COR_WINDOW} \
${ENRICHMENTS}"

# Boolean flags
if [ "$CORRECT_LAMBDA" == true ]; then
	CMD="$CMD --correctLambda"
fi

if [ "$EXLUDE_HLA" == true ]; then
	CMD="$CMD --excludeHla"
fi

if [ "$GENE_FORCE_NORMAL" == true ]; then
	CMD="$CMD --forceNormalGenePvalues"
fi

if [ "$PATH_FORCE_NORMAL" == true ]; then
	CMD="$CMD --forceNormalPathwayPvalues"
fi

echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] Starting Depict2 run"
echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] command = ${CMD}"
echo "=========================================================================="

eval $CMD

# Copy output to network
if [ "${USE_LOCAL_TMP}" == true ]; then
	# Create output dir if not exists
	mkdir -p ${OUTPUT}/${PHENOTYPE}_${VERSION}
	# Copying output to network drive
	cp -r ${TMP_OUT}/* ${OUTPUT}/${PHENOTYPE}_${VERSION}/
fi

echo "=========================================================================="
echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')]ran for $(( ($(date +%s) - ${START_TIME}) / 60)) minutes"

