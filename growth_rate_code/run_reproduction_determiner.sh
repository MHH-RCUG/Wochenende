#!/bin/bash

for csv in $(ls -d *_subsamples);do
  echo "working on sample "${csv%_subsample}""
  python3 main.py -p "$csv"
done