#!/bin/bash
# Colin Davenport, April 2020 - Jan 2021
# Run multiqc report
# Collect mapping stats from flagstat
# Run filter: Keep all lines in the bam.txt where column 3 (reads aligned)
# is greater than X (here probably 20). Good for idxstats files i.e. bam.txt files from Wochenende



# Run multiqc
echo "INFO:  Running multiqc"
#multiqc -f .

# Collate mapping stats
out="mapped_percent.txt"
echo "INFO:  Generating Wochenende mapping stats to $out"
echo "Wochenende mapping stats" > $out
for z in `ls *flagstat.txt`
	do
	echo -n $z "\t"  >> $out
	grep "mapped (" $z >> $out
done


echo "INFO:  Filtering and sorting  Wochenende output bam.txt files"
## Filter: $3>= means column3 (no. reads assigned to taxon) must have 20 or more reads
# Sorted descending in column 3 for within experiment clarity 
for i in `find . -name "*.bam.txt"`
        do
	awk -F "\t" '$3>=20' $i | sort -t$'\t' -k3 -nr > $i.filt.sort.csv
done

echo "INFO:  Script runbatch_metagen_awk_filter.sh completed"
