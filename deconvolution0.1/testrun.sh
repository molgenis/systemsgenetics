java -Xms8000M -Xmx10000M -jar deconvolution.jar \
    -o /Users/NPK/UMCG/tmp/decon/ \
    -c /Users/NPK/UMCG/tmp/decon/simulatedCellcount.txt \
    -e /Users/NPK/UMCG/tmp/decon/simulatedExpressionData.txt \
    -sn /Users/NPK/UMCG/tmp/decon/simulatedGeneSNPpair.txt \
    -g /Users/NPK/UMCG/tmp/decon/simulatedGenotypeData.txt \
    --output_coefficients \
    --whole_blood_qtl 
