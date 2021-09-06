#!/bin/bash

# Setup conda and directories using data parsed from config.yaml
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)
haybaler_dir=$HAYBALER_DIR
wochenende_dir=$WOCHENENDE_DIR
# Set and activate existing conda env
. $CONDA_SH_PATH
conda activate $WOCHENENDE_CONDA_ENV_NAME


csv_count=$(ls -d *_subsamples 2>/dev/null | wc -l)
if [[ $csv_count != 0 ]]
  then
  for csv in $(ls -d *_subsamples);do
    echo "working on sample ${csv%_subsample}"
    python3 main.py -p "$csv"
  done
fi