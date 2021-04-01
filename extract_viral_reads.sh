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


for i in *calmd.bam; do

	input=$i
	sec_input=${input%%.bam}
	
	# small windows for high resolution profiling of small references
	window=1000
	overlap=500
	covMax=999999999
	threads=8
	queue=short

	echo "INFO: Create new BAM file with reads from references specified in BED file $taxaToKeep only"
	
	# SLURM
    #Aligned reads in a region specified by a BED file
    srun -c $threads samtools view -@ $threads -b -h -L $taxaToKeep -o extract/$sec_input.filt.bam $input 
    srun -c 1 samtools index extract/$sec_input.filt.bam

	# Direct submission, not SLURM
	#samtools view -@ $threads -b -h -L $taxaToKeep -o extract/$sec_input.filt.bam $input 
    #samtools index extract/$sec_input.filt.bam



	echo "INFO: Extracting data from extracted specified BAMs"
	window_input_bam=extract/$sec_input.filt.bam	
	window_input_bam_prefix=${window_input_bam%%.bam}
	# Get coverage depth in tiny windows (eg 1kbp, set above) for small ref seqs, eg viruses
	# SLURM
	srun -c $threads -p $queue sambamba depth window -t $threads --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${window_input_bam_prefix}.bam > ${window_input_bam_prefix}_cov_window.txt &

	# Direct submission, not SLURM
	#sambamba depth window -t $threads --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${window_input_bam_prefix}.bam > ${window_input_bam_prefix}_cov_window.txt &



done
wait


echo "INFO: Prepare filtered data for plotting"
cd extract/
bash runbatch_metagen_window_filter.sh


echo "INFO: Making plots"
runbatch_wochenende_plot.sh >/dev/null 2>&1
cd $bamDir