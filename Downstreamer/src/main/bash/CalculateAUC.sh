#!/bin/bash
#SBATCH --time="05:59:00"
#SBATCH --memory="32G"

# java -jar ../../../depict2/Depict2-2.0.51-SNAPSHOT/Depict2.jar \
#--mode CORE_GENE_AUC \
#-g educational_attainment_2018_30038396_hg19_genePvalues \
#-o auc.txt \
#-pd hpo=../../../reference_datasets/human_b37/hpo/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix


set -e
ml Java


MEMORY="30G"
CORES=1

# DEPICT2 global settings
BUNDLE_DIR="/groups/umcg-wijmenga/tmp04/projects/depict2/depict2_bundle"
VERSION="51"
TYPE="CORE_GENE_AUC"
DEPICT2="${BUNDLE_DIR}/depict2/Depict2-2.0.${VERSION}-SNAPSHOT/Depict2.jar"

INPUT_FILE=$1

echo "[INFO] Proccesing ${INPUT_FILE}"

CMD="java -Xmx${MEMORY} -XX:ParallelGCThreads=${CORES} -jar ${DEPICT2} \
-m ${TYPE} \
-g ${BUNDLE_DIR}/output/height_paper/${INPUT_FILE}_${VERSION}/${INPUT_FILE}_intermediates/Coregulation_Enrichment_zscoreExHla \
-o ${BUNDLE_DIR}/output/height_paper/${INPUT_FILE}_${VERSION}/${INPUT_FILE}_intermediates/${INPUT_FILE}_coreGene_hpoAUC \
-pd hpo=${BUNDLE_DIR}/reference_datasets/human_b37/hpo/PhenotypeToGenes.txt_matrix.txt"

#-g ${BUNDLE_DIR}/output/maf_filtered/${INPUT_FILE}_${VERSION}/${INPUT_FILE}_intermediates/Coregulation_Enrichment_zscoreExHla \
#-pd hpo=${BUNDLE_DIR}/reference_datasets/human_b37/hpo/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix"
#-pd hpo=${BUNDLE_DIR}/reference_datasets/human_b37/hpo/PhenotypeToGenes.txt_matrix.txt"

#echo $CMD
eval $CMD


CMD="java -Xmx${MEMORY} -XX:ParallelGCThreads=${CORES} -jar ${DEPICT2} \
-m ${TYPE} \
-g ${BUNDLE_DIR}/output/height_paper/${INPUT_FILE}_${VERSION}/${INPUT_FILE}_genePvalues \
-o ${BUNDLE_DIR}/output/height_paper/${INPUT_FILE}_${VERSION}/${INPUT_FILE}_intermediates/${INPUT_FILE}_genePvalue_hpoAUC \
-pd hpo=${BUNDLE_DIR}/reference_datasets/human_b37/hpo/PhenotypeToGenes.txt_matrix.txt"

eval $CMD


# Make gene pvalues to txt file
CMD="java -Xmx${MEMORY} -XX:ParallelGCThreads=${CORES} -jar ${DEPICT2} \
-m CONVERT_BIN \
-g ${BUNDLE_DIR}/output/height_paper/${INPUT_FILE}_${VERSION}/${INPUT_FILE}_genePvalues \
-o ${BUNDLE_DIR}/output/height_paper/${INPUT_FILE}_${VERSION}/${INPUT_FILE}_genePvalues"

eval $CMD


