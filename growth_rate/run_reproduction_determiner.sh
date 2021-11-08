#!/bin/bash
# Run the growth rate reproduction determiner module
# Sophia Poertner 2021

csv_count=$(ls -d *_subsamples 2>/dev/null | wc -l)
if [[ $csv_count != 0 ]]
  then
  for csv in $(ls -d *_subsamples);do
    echo "working on sample ${csv%_subsample}"
    python3 main.py -p "$csv"
  done
fi