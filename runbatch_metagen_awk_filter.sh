#!/bin/bash
# Colin Davenport, August 2018

# Keep all lines Where col3 is greater than X (currently 200). Good for idxstats files eg metagenomics

# sorted descending in col3 for within experiment clarity
for i in `find . -name "*.bam.txt"`
        do
	awk -F "\t" '$3>=20' $i | sort -t$'\t' -k3 -nr > $i.filt.sort.csv
done

# unsorted to maintain row comparability over experiments
#for i in `find . -name "*.bam.txt"`
#        do
#        awk -F "\t" '$3>=20' $i > $i.filt.unsort.csv
#done


