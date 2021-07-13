java -jar Depict2.jar \
	-g demoData/demoGwas \
	-ge demoData/demoGenes.txt \
	-m RUN \
	-o demoResult/demo \
	-p 1000 \
	-r demoData/demoGenotypes.vcf.gz \
	-rs demoData/demoSamples.txt \
	-R VCF \
	-t 1 \
	-v 0.95 \
	-w 1000 \
	--debug

