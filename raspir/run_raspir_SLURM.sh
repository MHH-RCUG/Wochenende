#!/bin/bash

# run_raspir_SLURM.sh
# last updated: 28 July 2021

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 30000

# set run x cpus
#SBATCH --cpus-per-task 8

# set name of job
#SBATCH --job-name=raspir_run

# Add miniconda3 to PATH. TODO - autodetection
. /mnt/ngsnfs/tools/miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate raspir_env >> /dev/null

for i in `ls *raspir.csv`

        do
                echo $i
                python raspir.py $i ${i%.csv}

done
