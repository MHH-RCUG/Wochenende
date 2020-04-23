#!/bin/bash
# Colin Davenport, Jan 2020

# Operates on window output files from runbatch_sambamba_depth.sh
# Keep all lines without Where col3 is greater than X (currently 200). Good for idxstats files eg metagenomics

# sorted descending in col5 for within experiment clarity and to find regions with excessive counts 
# eg false positives, illumina adapters in assemblies etc

for i in `ls *window.txt`
        do
	# filter out human reads
	grep -v "1_1_1"  "$i" > "$i".filt.csv

	# filter out human reads and sort by total reads per window
	grep -v "1_1_1"  "$i" | sort -t$'\t' -k5 -nr > "$i".filt.sort.csv
done

