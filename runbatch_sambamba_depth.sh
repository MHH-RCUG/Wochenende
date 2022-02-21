#!/bin/bash
# Erik Wessels, Colin Davenport Jan 2020 - Nov 2021
# Check window coverage on Wochenende sorted dup.bam output
# Use output for Python script to check coverage distribution

# Setup SLURM using data parsed from config.yaml
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)
# Setup job scheduler
# use SLURM job scheduler (yes, no)
if [[ "${USE_CUSTOM_SCHED}" == "yes" ]]; then
    #echo USE_CUSTOM_SCHED set"
    scheduler=$CUSTOM_SCHED_CUSTOM_PARAMS
fi
if [[ "${USE_SLURM}" == "yes" ]]; then
    #echo USE_SLURM set"
    scheduler=$SLURM_CUSTOM_PARAMS
fi


# Actually run for each BAM file
for i in `ls *calmd.bam`; do
	input=$i
	sec_input=${input%%.bam}
	#sec_in_bam=${input%%.bam}

	window=100000
	overlap=50000
	covMax=999999999
	
	# Get coverage depth in windows
	# x threads, Windows 100000, overlap 50000, -c minimum coverage. 
	# SLURM
	$scheduler $path_sambamba depth window -t $THREADS --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${sec_input}.bam > ${sec_input}_cov_window.txt &

done
wait
