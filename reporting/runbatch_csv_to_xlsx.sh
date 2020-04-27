#!/bin/bash
set -eo pipefail
shopt -s nullglob

# Fabian Charly Friedrich
# Start from a folder containing csv files from the nextflow blast pipeline
# bash runbatch_csv_to_xlsx.sh

# activate conda env
. /mnt/ngsnfs/tools/miniconda3/etc/profile.d/conda.sh
conda activate wochenende >> /dev/null

# wochenende output files
for file in ./*rep.*.csv
do
        echo "$file"
	      srun -p short -c 1 python3 csv_to_xlsx_converter.py "$file" &
	wait
done

# kraken2 and krakenuniq output files 
for file in ./*.txt
do
        echo "$file"
              srun -p short -c 1 python3 csv_to_xlsx_converter.py "$file" &
	wait
done

# Nextflow blast output files
for file in ./*annot.csv
do
	conda activate nextflow >> /dev/null
        echo "$file"
              srun -p short -c 1 python3 csv_to_xlsx_converter.py "$file" &
        wait
done

