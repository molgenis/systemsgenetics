awk -v OFS="\t" '$1=$1' $1 > $1.tab
bgzip $1.tab
tabix -b 3 -e 3 $1.tab.gz