# Cluster settings
TIME="23:59:00"
MEMORY="100G"
CORES=8
TMPDIR_SIZE=16G
USE_LOCAL_TMP=true

# DEPICT2 global settings
BUNDLE_DIR="/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle"
VERSION="37"
TYPE="RUN2"
DEPICT2="${BUNDLE_DIR}/depict2/Depict2-2.0.${VERSION}-SNAPSHOT/Depict2.jar"

# When using default bundle these can be left default
REFERENCE_SET="${BUNDLE_DIR}/reference_datasets/human_b37/"
REFERENCE_GENOTYPES="${REFERENCE_SET}/1000G_EUR_noFIN/HarmoniserAllChrs"
REFERENCE_SAMPLES="${REFERENCE_SET}/1000G_EUR_noFIN/phase3_Europeans.txt"
REFERENCE_TYPE="TRITYPER"
REFERENCE_ENSEMBL="${REFERENCE_SET}/ensgR75_skewness025.txt"
PATHWAY_DATABASES="${REFERENCE_SET}/pathway_databases/"

# Technical parameters
VARIANT_PRUNING=0.95
GENE_PRUNING=0.9
GENE_WINDOW=50000
GENE_FORCE_NORMAL=false
PATH_FORCE_NORMAL=false
EXLUDE_HLA=true
CORRECT_LAMBDA=false
GENE_COR_WINDOW=2000000

# Permutation settings
PERMUTATIONS=100000
PERMUTATION_RESCUE=1000000000
PERMUTATION_PATHWAY=10000
PERMUTATION_GENECOR=10000

# Enrichments
ENRICHMENTS="\
-pd Coregulation=${PATHWAY_DATABASES}/geneCoregulationZscore \
-pd Reactome=${PATHWAY_DATABASES}/reactome_prepared_predictions \
-pd GO_P=${PATHWAY_DATABASES}/go_P_prepared_predictions \
-pd GO_C=${PATHWAY_DATABASES}/go_C_prepared_predictions \
-pd GO_F=${PATHWAY_DATABASES}/go_F_prepared_predictions \
-pd KEGG=${PATHWAY_DATABASES}/kegg_prepared_predictions \
-pd HPO=${PATHWAY_DATABASES}/hpo_prepared_predictions \
-pd gtex=${PATHWAY_DATABASES}/gtexExpression \
-pd eigen=${PATHWAY_DATABASES}/eigenvectors_1588 \
-pd eigen_eQTLGen=${PATHWAY_DATABASES}/eigenvectors_eQTLGen"
