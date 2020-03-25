#!/bin/bash
# Colin Davenport, August 2018

# Keep all lines Where col3 is greater than X (here probably 20, or 200). Good for idxstats files eg metagenomics

# Sorted descending in column 3 for within experiment clarity
for i in `find . -name "*.bam.txt"`
        do
	awk -F "\t" '$3>=20' $i | sort -t$'\t' -k3 -nr > $i.filt.sort.csv
done

# Unsorted to maintain row comparability over experiments.
# Paste output together with paste -d"" unsort1.csv unsort2.csv > summary.csv
#for i in `find . -name "*.bam.txt"`
#        do
#        awk -F "\t" '$3>=20' $i > $i.filt.unsort.csv
#done
