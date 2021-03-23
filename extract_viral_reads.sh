#!/bin/bash
# Colin Davenport Jan 2020 - March 2021
# Gets only header and aligned reads from a bam file, writes to extract/infile.filt.bam
# Extracted viral reads are typically used for high resolution plotting
# Use script get_wochenende.sh to set up the extract directory
# Usage: bash extract_viral_reads.sh in a directory containing eg *calmd.bam files

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

	# Create new BAM file with reads from references specified in a BED file only
	
	# SLURM
    #Aligned reads in a region specified by a BED file
    srun -c $threads samtools view -@ $threads -b -h -L $taxaToKeep -o extract/$sec_input.filt.bam $input 
    srun -c 1 samtools index extract/$sec_input.filt.bam

	# Direct submission, not SLURM

	window_input_bam=extract/$sec_input.filt.bam	
	window_input_bam_prefix=${window_input_bam%%.bam}


	# Get coverage depth in tiny windows for small ref seqs, eg viruses
	# 8 threads, Windows 1000, overlap 500, -c minimum coverage. 
	# SLURM
	srun -c $threads -p $queue sambamba depth window -t $threads --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${window_input_bam_prefix}.bam > ${window_input_bam_prefix}_cov_window.txt &

	# Direct submission, not SLURM
	#sambamba depth window -t $threads --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${window_input_bam_prefix}.bam > ${window_input_bam_prefix}_cov_window.txt &


done
wait