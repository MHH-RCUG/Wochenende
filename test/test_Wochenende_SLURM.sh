#!/bin/bash
## A SLURM sbatch script which is part of the Wochenende Pipeline https://github.com/MHH-RCUG/Wochenende for testing

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 40000

# set run x cpus
#SBATCH --cpus-per-task 4

# set max wallclock time
#SBATCH --time=47:00:00

# set name of job
#SBATCH --job-name=Wochenende_test

# Add miniconda3 to PATH
. /mnt/ngsnfs/tools/miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate wochenende >> /dev/null

pytest --basetemp=/ngsssd1/scheitht/tmp/ -n 4 test/
