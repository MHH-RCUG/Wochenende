#!/bin/bash

# Prepare Wochenende bam files for the bed to csv script. Makes  bed out of Wochenende bam files. Filter out mouse human and mito.
# Author: Sophia Poertner, 2021

# covert bam to bed
for bam in *.bam
do
	./bedtools bamtobed -i "$bam" > "${bam%.bam}".unfiltered.bed
done

for bed in *.bed
do
	#exclude mouse, human, mito
  grep -v "^chr" "$bed" | grep -v "^1_1_1" > "${bed%.unfiltered.bed}".bed
	rm "$bed"  # remove temp file
	echo "starting bed to csv for file $bed"
	python3 bed_to_pos_csv.py -i "${bed%.unfiltered.bed}".bed -p .
	echo "completed file $bed"
  # cleanup
  if [[ ! -d "${bed%.unfiltered.bed}_subsamples" ]]
    then
    mkdir "${bed%.unfiltered.bed}_subsamples"
  fi

  csv_count=$(ls -1 "${bed%.unfiltered.bed}"*.csv 2>/dev/null | wc -l)
  if [[ $csv_count != 0 ]]
      then
      mv "${bed%.unfiltered.bed}"*.csv "${bed%.unfiltered.bed}_subsamples"
  fi
done
