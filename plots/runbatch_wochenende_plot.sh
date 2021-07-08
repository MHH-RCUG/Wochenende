#!/bin/bash
# Colin Davenport, Jan 2020 - Jul 2021
# Run Wochenende plotting as batch, add SLURM

# Set variables
cpus=1
#with SLURM srun, default
slurm_srun="srun -c $cpus"

# Save output log in directory containing bams and preprocess script
output_log="plot_"$(date +%s)".log"
echo "INFO: output_log: " $output_log

# mv existing images to a backup directory
images_backup="images_old_"$(date +%s)
mv images $images_backup

for i in $(ls *cov_window.txt.filt.csv)
        do
		
		# filter out human reads and sort by total reads per window
		#python wochenende_plot.py "$i" >> $output_log
		# slurm version
		$slurm_srun python wochenende_plot.py "$i" >> $output_log &

done
wait

