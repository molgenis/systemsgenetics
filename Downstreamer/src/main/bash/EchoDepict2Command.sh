#!/bin/bash

# Usage
# ./Depict2Run.sh <traitfile> <output folder> <parameterfile>
# Path to the depict2 compatible zscore matrix or gene pvalues when using RUN2
# Parameter file is a bashscript containing the paramters voor een run

# If any exit code in this script is not 0, quit
set -e

# Load java
ml Java
echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] Using java version:	$(which java)"

START_TIME=$(date +%s)

# User defined settings for this run
INPUT="$1"
OUTPUT="$2"
PARAMS="$3"
PHENOTYPE="$(basename $INPUT)"

# Load the parameter file
source $PARAMS

# Utility needed for making the symlinks credit to:
# https://stackoverflow.com/questions/3915040/bash-fish-command-to-print-absolute-path-to-a-file
get_abs_filename() {
  # $1 : relative filename
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

if [ -f "${DEPICT2}" ]; then
	echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] Using version: ${VERSION}"
else 
	echo "[ERROR - $(date '+%Y-%m-%d %H:%M:%S')] Depict2 jar not found at: ${DEPICT2}"
	exit 1;
fi

TMP_IN=${INPUT}
TMP_OUT=${OUTPUT}/${PHENOTYPE}_${VERSION}${OUTDIR_SUFFIX}

# Construct depict command
# Add 4G overhead on the memory to avoid out of memory errors during spikes
# Reserve 2 cores for GC
CMD="java -Xms$((MEMORY - 4))G -Xmx$((MEMORY - 4))G -XX:ParallelGCThreads=2 -jar ${DEPICT2} \
-g ${TMP_IN} \
-m ${TYPE} \
-o ${TMP_OUT}/${PHENOTYPE} \
-ge ${REFERENCE_ENSEMBL} \
-r ${REFERENCE_GENOTYPES} \
-rs ${REFERENCE_SAMPLES} \
-R ${REFERENCE_TYPE} \
-t $((CORES - 2)) \
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
if [ "$FILTER_MAF" == true ]; then
        CMD="$CMD --maf ${MAF}"
fi

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

if [ "$SAVE_EXCEL" == true ]; then
        CMD="$CMD --saveExcel"
fi


echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] Starting Depict2 run"
echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] command = ${CMD}"
echo "=========================================================================="

echo $CMD


