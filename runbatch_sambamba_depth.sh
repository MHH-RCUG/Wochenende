#!/bin/bash
# Erik Wessels, Colin Davenport Jan 2020 - March 2021
# Check window coverage on Wochenende sorted dup.bam output
# Use output for Python script to check coverage distribution

for i in *calmd.bam; do
	input=$i
	sec_input=${input%%.bam}
	#sec_in_bam=${input%%.bam}

	window=100000
	overlap=50000
	covMax=999999999
	threads=1
	queue=short

	# Get coverage depth in windows
	# 8 threads, Windows 100000, overlap 50000, -c minimum coverage. 
	# SLURM
	srun -c $threads -p $queue sambamba depth window -t $threads --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${sec_input}.bam > ${sec_input}_cov_window.txt &

	# Direct submission, not SLURM
	#sambamba depth window -t $threads --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${sec_input}.bam > ${sec_input}_cov_window.txt &

done
wait
