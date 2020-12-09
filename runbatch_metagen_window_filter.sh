#!/bin/bash
# Colin Davenport, April-Nov 2020

# Operates on window output files from runbatch_sambamba_depth.sh
# Preprocesses files for the Wochenende plot tool
# sorted descending in col5 for within experiment clarity and to find regions with excessive counts 
# eg false positives, illumina adapters in assemblies etc
# Filters human and mouse chromosomes out of coverage files (starting with 1_1_1 or chr)

for i in `ls *window.txt`
        do

	# Change following line if needed to human or mouse
	species=human

	if [ $species == "human" ]
		then

		echo "Filtering windows, excluding" $species

		# filter out human reads
		grep -v "^1_1_1"  "$i" > "$i".filt.csv

		# filter out human reads and sort by total reads per window
		grep -v "^1_1_1"  "$i" | sort -t$'\t' -k5 -nr > "$i".filt.sort.csv


	fi
	if [ $species == "mouse" ]
		then

		echo "Filtering windows, excluding" $species

		# filter out mouse reads
        	grep -v "^chr"  "$i" > "$i".filt.csv

		# filter out mouse reads and sort by total reads per window
        	grep -v "^chr"  "$i" | sort -t$'\t' -k5 -nr > "$i".filt.sort.csv
	fi

done

