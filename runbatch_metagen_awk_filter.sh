#!/bin/bash
# Colin Davenport, April 2020
# Run multiqc report
# Collect mapping stats from flagstat
# Run filter: Keep all lines in the bam.txt where column 3 (reads aligned)
# is greater than X (here probably 20, or 200). Good for idxstats files i.e. bam.txt files from Wochenende

# Quick cleanup, only keep *.dup.bam.txt
rm *.01mm.bam.txt *.01mm.dup.calmd.bam.txt *.mq30.bam.txt


# Run multiqc
multiqc -f .

# Get mapping stats
out=mapped_percent.txt
echo "INFO:  Generating Wochenende mapping stats to $out"
echo "Wochenende mapping stats" > $out
for z in `ls *flagstat.txt`
	do
	echo -n $z "\t"  >> $out
	grep "mapped (" $z >> $out
done


echo "INFO:  Filtering and sorting  Wochenende output bam.txt files"
## Filter:
# Sorted descending in column 3 for within experiment clarity
for i in `find . -name "*.bam.txt"`
        do
	awk -F "\t" '$3>=20' $i | sort -t$'\t' -k3 -nr > $i.filt.sort.csv
done

## Filter2:
# Unsorted to maintain row comparability over experiments.
# Paste output together with paste -d"" unsort1.csv unsort2.csv > summary.csv
#for i in `find . -name "*.bam.txt"`
#        do
#        awk -F "\t" '$3>=20' $i > $i.filt.unsort.csv
#done


echo "INFO:  Script completed"
