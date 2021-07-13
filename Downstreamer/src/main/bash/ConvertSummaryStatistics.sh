#!/bin/bash
#SBATCH --time=05:59:00
#SBATCH --mem=100G
# Converts SNP <tab> P formatted summary statistics to a binary zscore matrix for use with Depict 2

set -e
ml Java

# Cluster settings
TIME="23:59:00"
MEMORY="100G"
CORES=1

# DEPICT2 global settings
BUNDLE_DIR="/groups/umcg-wijmenga/tmp04/projects/depict2/depict2_bundle"
VERSION="51"
TYPE="RUN"
DEPICT2="${BUNDLE_DIR}/depict2/Depict2-2.0.${VERSION}-SNAPSHOT/Depict2.jar"

while read line;
do

INPUT_FILE="$(echo $line | awk '{print $1}')"
OUTPUT_FILE="$(echo $line | awk '{print $2}')"

echo "[INFO] Proccesing ${INPUT_FILE}"

CMD="java -Xmx${MEMORY} -XX:ParallelGCThreads=${CORES} -jar ${DEPICT2} \
-m CONVERT_TXT \
-p2z \
-g ${BUNDLE_DIR}/summary_statistics/plaintext/SNP_P_FORMATTED/${INPUT_FILE} \
-o ${BUNDLE_DIR}/summary_statistics/binary_matrix/${OUTPUT_FILE}"

eval $CMD

done < traitfiles/files_to_convert.tsv
