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

# If the outputdir already exists, make a new one with the date
#if [ -d "${OUTPUT}/${PHENOTYPE}_${VERSION}" ]; then
#	OUTDIR_SUFFIX="_$(date +%m-%d-%Y-%H-%M)"
#else
#	OUTDIR_SUFFIX=""
#fi

# If using a local tmp dir on the node, create the folder structure and copy files
if [ "${USE_LOCAL_TMP}" == true ]; then
	# TMPDIR on local storage for speedy IO
	# Copy input zscores to tmpdir
	mkdir -p ${TMPDIR}/${PHENOTYPE}_${VERSION}${OUTDIR_SUFFIX}
	cp ${INPUT}* ${TMPDIR}/${PHENOTYPE}_${VERSION}${OUTDIR_SUFFIX}

	# Write output to local storage
	TMP_IN=${TMPDIR}/${PHENOTYPE}_${VERSION}${OUTDIR_SUFFIX}/${PHENOTYPE}
	TMP_OUT="${TMPDIR}/${PHENOTYPE}_${VERSION}${OUTDIR_SUFFIX}"
	echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] copied input to tmpdir"
else 
	TMP_IN=${INPUT}
	TMP_OUT=${OUTPUT}/${PHENOTYPE}_${VERSION}${OUTDIR_SUFFIX}
fi

mkdir -p ${TMP_OUT}

# When running RUN2 make symlynks to the new output folder if needed
if [ ! "${OUTDIR_SUFFIX}" == "" ]; then
	if [ "${TYPE}" == "RUN2" ]; then
		echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] Creating symlinks for gene pvalues"
		#echo "[WARN - $(date '+%Y-%m-%d %H:%M:%S')] Only works with relative paths ATM"
		ABS_INPUT_PATH=$(get_abs_filename "${INPUT}")
		ln -s ${ABS_INPUT_PATH}_gene* ${TMP_OUT}/	
	fi
fi

if [ "${TYPE}" == "RUN2" ]; then
	TMP_IN="${TMP_IN}_genePvalues"
fi

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

if [ "$DEBUG" == true ]; then
        CMD="$CMD -d"
fi


echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] Starting Depict2 run"
echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] command = ${CMD}"
echo "=========================================================================="
eval $CMD

# Copy output to network drive
if [ "${USE_LOCAL_TMP}" == true ]; then
	# Create output dir if not exists
	mkdir -p ${OUTPUT}/${PHENOTYPE}_${VERSION}${OUTDIR_SUFFIX}
	# Copying output to network drive
	cp -r ${TMP_OUT}/* ${OUTPUT}/${PHENOTYPE}_${VERSION}${OUTDIR_SUFFIX}/
fi

echo "=========================================================================="
echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')]ran for $(( ($(date +%s) - ${START_TIME}) / 60)) minutes"


