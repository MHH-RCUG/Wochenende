#!/bin/bash

# Prepare Wochenende bam files for the bed to csv script (input for reproduction determiner)
# Links input files from one higher directory. Converts bam to bed.
# Filter out mouse human and mito chromosomes.
# Author: Sophia Poertner, 2021
# Usage: bash run_bed_to_csv.sh input.bam

echo "Version 0.13"

# Changelog
# 0.13 - add usage, correct runbatch_bed_to_csv.sh SLURM submission
# 0.12 - get input file from sbatch script for speedup
# 0.11 - link in bam, bam.txt and bai files, unlink later
# 0.10 - remove bedtools binary and use conda bedtools

# check if input exists
count_bam=`ls -1 ../*calmd.bam 2>/dev/null | wc -l`
count_bai=`ls -1 ../*calmd.bam.bai 2>/dev/null | wc -l`
count_bam_txt=`ls -1 ../*calmd.bam.txt 2>/dev/null | wc -l`

if [ $count_bam != 0 ]  && [ $count_bai != 0 ]  && [ $count_bam_txt != 0 ]
  then
  # link bam, bai and bam.txt files found by ls command to current directory
  ln -s $1 .
  ln -s ${1%bam}bam.txt .
  ln -s ${1%bam}bam.bai .
  ls *
  bam=${1/..\//}

  bedtools bamtobed -i "$bam" > "${bam%.bam}.base.bed"

  bed="${bam%.bam}.base.bed"
  # filter - exclude mouse, human, mito chromosomes
  grep -v "^chr" "$bed" | grep -v "^1_1_1" > "${bed%.bed}.filt.bed"
  echo "INFO: Starting bed to csv for file $bed"
  python3 bed_to_pos_csv.py -i "${bed%.bed}.filt.bed" -p .
  echo "INFO: Completed file $bed"

  # cleanup
  #rm "$bed"  # remove temp file
  if [[ ! -d "${bed%.bed}_subsamples" ]]
    then
    mkdir "${bed%.bed}_subsamples"
  fi

  csv_count=$(ls -1 "${bed%.bed}"*.csv 2>/dev/null | wc -l)
  if [[ $csv_count != 0 ]]
      then
      mv "${bed%.bed}"*.csv "${bed%.bed}_subsamples"
  fi

  echo "cleanup: unlink bam, bai and bam.txt files"
  # unlink bam, txt and bai files
  unlink $bam
  unlink ${bam%bam}bam.txt
  unlink ${bam%bam}bam.bai
else
  echo "no bam.txt and bai found for input bam. Can't convert to pos.csv"
fi
