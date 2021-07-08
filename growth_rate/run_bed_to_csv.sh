#!/bin/bash

# Prepare Wochenende bam files for the bed to csv script. Makes bed out of Wochenende bam files. Filter out mouse human and mito.
# Author: Sophia Poertner, 2021

echo "Version 0.11"

#Changelog
#0.11 -link in bam and bai files, unlink later
#0.1 - remove bedtools binary and use conda bedtools

#link bam and bai files found by ls command to current directory
for i in `ls ../*calmd.bam.bai`
        do
        ln -s $i .
done
for i in `ls ../*calmd.bam`
        do
        ln -s $i .
done

ls *

# convert *calmd.bam to bed
bam_count=$(ls -1 *.bam 2>/dev/null | wc -l)
if [[ $bam_count != 0 ]]
  then
  for bam in *.bam
  do
    #bedtools bamtobed -i "$bam" > "${bam%.bam}".unfiltered.bed
    bedtools bamtobed -i "$bam" > "${bam%.bam}.unfiltered.bed"
  done

  for bed in *.unfiltered.bed
  
  do
    #exclude mouse, human, mito chromosomes
    grep -v "^chr" "$bed" | grep -v "^1_1_1" > "${bed%.unfiltered.bed}".bed
    echo "INFO: Starting bed to csv for file $bed"
    python3 bed_to_pos_csv.py -i "${bed%.unfiltered.bed}".bed -p .
    echo "INFO: Completed file $bed"

    # cleanup
    #rm "$bed"  # remove temp file
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
fi


#unlink bam and bai files
for i in `ls *calmd.bam.bai`
        do
        unlink $i
done
for i in `ls *calmd.bam`
        do
        unlink $i
done
