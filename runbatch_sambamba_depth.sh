#!/bin/bash
# Erik Wessels, Colin Davenport Jan 2020
# Check window coverage on Wochenende sorted dup.bam output
# Use output for Python script to check coverage distribution

for i in *dup.bam; do
	input=$i
	sec_input=${input%%.bam}
	#sec_in_bam=${input%%.bam}

	window=100000
	overlap=50000
	covMax=999999999
	threads=8
	queue=short

	# Sort - coo means coordinate sorting
        #sambamba sort $input -o ${sec_input}_coo.bam

	# Get coverage depth in windows
	# 8 threads, Windows 100000, overlap 50000, -c minimum coverage. 
	# SLURM
	srun -c $threads -p $queue sambamba depth window -t $threads --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${sec_input}.bam > ${sec_input}_cov_window.txt &

	# Direct submission, not SLURM
	#sambamba depth window -t $threads --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${sec_input}.bam > ${sec_input}_cov_window.txt &

	# if using sorting step _coo from above
	#srun -c 8 sambamba depth window -t 8 -w 10000 --overlap=5000 -c 0.00001 ${sec_input}_coo.bam > ${sec_input}_cov_window.txt

done
wait


#echo "Starting sleep phase 1500s. If jobs done on SLURM then the wait command worked, and you can Ctrl-C out of this"
#sleep 1500
