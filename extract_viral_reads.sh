#!/bin/bash
# Colin Davenport Jan 2020 - March 2021
# Gets only header and aligned reads from a bam file, writes to extract/infile.filt.bam
# Extracted viral reads are typically used for high resolution plotting
# Use script get_wochenende.sh to set up the extract directory
# Usage: bash extract_viral_reads.sh in a directory containing eg *calmd.bam files

# get current dir containing Wochenende BAM and bam.txt output
bamDir=$(pwd)

# Refresh Wochenende plot files on calling this script
cp plots/*plot* extract/
cp runbatch_metagen_window_filter.sh extract/

taxaToKeep="extract/viruses_2021_02.bed"

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

for i in *calmd.bam; do

	input=$i
	sec_input=${input%%.bam}
	
	# small windows for high resolution profiling of small references
	window=1000
	overlap=500
	covMax=999999999

	echo "INFO: Create new BAM file with reads from references specified in BED file $taxaToKeep only"
	
	#Aligned reads in a region specified by a BED file
    $scheduler $path_samtools view -@ $THREADS -b -h -L $taxaToKeep -o extract/$sec_input.filt.bam $input 
    $scheduler $path_samtools index extract/$sec_input.filt.bam


	echo "INFO: Extracting data from extracted specified BAMs"
	window_input_bam=extract/$sec_input.filt.bam	
	window_input_bam_prefix=${window_input_bam%%.bam}
	# Get coverage depth in tiny windows (eg 1kbp, set above) for small ref seqs, eg viruses
	# SLURM
	$scheduler sambamba depth window -t $THREADS --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${window_input_bam_prefix}.bam > ${window_input_bam_prefix}_cov_window.txt &



done
wait


echo "INFO: Prepare filtered data for plotting"
cd extract/
bash runbatch_metagen_window_filter.sh


echo "INFO: Making plots"
runbatch_wochenende_plot.sh >/dev/null 2>&1
cd $bamDir