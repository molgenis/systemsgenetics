VERSION="37"
BUNDLE_DIR="/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_depict2/depict2_bundle"
MEMORY="100G"
CORES=1

TYPE="CONVERT_TXT"

CMD="java -Xmx${MEMORY} -XX:ParallelGCThreads=${CORES} -jar ${DEPICT2} \
-g ${TMP_IN} \
-m ${TYPE} \
-o ${TMP_OUT}/${PHENOTYPE} \