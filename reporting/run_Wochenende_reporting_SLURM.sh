#!/bin/bash
## Supply the Wochenende bam.txt input as arg1, bash run_Wochenende_reporting_SLURM.sh in.bam.txt


echo "Input bam: " $1
bamtxt=$1

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

# Set and activate existing conda env
. $CONDA_SH_PATH
conda activate $WOCHENENDE_CONDA_ENV_NAME


# read reference path from ref.tmp file in the reporting (current) directory 
ref=$(<ref.tmp)
echo "$ref"

if [[ $ref == "" ]]; then
	echo "INFO: Error in reporting - reference file is empty or not present. Rerun Wochenende to correct this."
	exit 1
fi

# Run script


#Former script runbatch_Wochenende_reporting was removed and functionality included here
echo "INFO: Starting batch reporting"

for bamtxt in `ls *.bam.txt`
        do
        # run 
        $scheduler python3 basic_reporting.py --input_file $bamtxt --reference $ref --sequencer illumina --output_name $bamtxt >/dev/null 2>&1 &
        #bash run_Wochenende_reporting_SLURM.sh $bamtxt >/dev/null 2>&1 &
done

echo "INFO: Waiting for SLURM scripts to complete"
wait

echo "INFO: Completed batch reporting"
