#!/bin/bash

# Usage
# ./SubmitDepict2Run.sh <traitfile> <output dir> <parameterfile>
# Traitfile contains per row a path to the depict2 compatible zscore matrix
# Parameter file is a bashscript containing the paramters voor een run

# TODO: Make these named
TRAITS="$1"
OUTPUT_DIR="$2"
PARAMS="$3"

source ${PARAMS}

mkdir -p ${OUTPUT_DIR}/logs

# If the outputdir already exists, add the date as a suffix
#if [ -d "${OUTPUT}/${PHENOTYPE}_${VERSION}" ]; then
#	echo "[INFO - $(date '+%Y-%m-%d %H:%M:%S')] Output folder already exists, appending date to the end"
#	OUTDIR_SUFFIX="_$(date +%m-%d-%Y-%H-%M)"
#else
#	OUTDIR_SUFFIX=""
#fi

# Construct sbatch command
CMD="sbatch --time=${TIME} --mem=${MEMORY}G --ntasks=${CORES} --nodes=1"

if [ "${USE_LOCAL_TMP}" == true ]; then
	CMD="${CMD} --tmp=${TMPDIR_SIZE} "
fi

# Submit all the runs in traitfile
while read traitfile; do

NAME="$(basename $traitfile)"
CMD_CUR="${CMD} \
--job-name=Depict2_v${VERSION}_${NAME}${OUTDIR_SUFFIX} \
--output=${OUTPUT_DIR}/logs/Depict2_v${VERSION}_${NAME}${OUTDIR_SUFFIX}.out \
--err=${OUTPUT_DIR}/logs/Depict2_v${VERSION}_${NAME}${OUTDIR_SUFFIX}.err \
--exclude=gs-vcompute03 \
Depict2Run.sh ${traitfile} ${OUTPUT_DIR} ${PARAMS}" 

echo $CMD_CUR
eval $CMD_CUR

done < $TRAITS; 

