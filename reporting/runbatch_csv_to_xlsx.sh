#!/bin/bash
set -eo pipefail
shopt -s nullglob

# Fabian Charly Friedrich and Colin Davenport
# Start from a folder containing csv files 
# bash runbatch_csv_to_xlsx.sh

# Setup SLURM using data parsed from config.yaml
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)
# Setup job scheduler
# use SLURM job scheduler (yes, no)
if [[ "${USE_CUSTOM_SCHED}" == "yes" ]]; then
    #echo USE_CUSTOM_SCHED set"
    scheduler=$CUSTOM_SCHED_CUSTOM_PARAMS_SINGLECORE
fi
if [[ "${USE_SLURM}" == "yes" ]]; then
    #echo USE_SLURM set"
    scheduler=$SLURM_CUSTOM_PARAMS_SINGLECORE
fi


# wochenende output files
for file in ./*rep.*.csv
do
        echo "$file"
	      $scheduler python3 csv_to_xlsx_converter.py "$file" &
	wait
done

# kraken2 and krakenuniq output files 
for file in ./*.txt
do
        echo "$file"
              $scheduler python3 csv_to_xlsx_converter.py "$file" &
	wait
done

# Nextflow blast output files
for file in ./*annot.csv
do
	conda activate nextflow >> /dev/null
        echo "$file"
              $scheduler python3 csv_to_xlsx_converter.py "$file" &
        wait
done

