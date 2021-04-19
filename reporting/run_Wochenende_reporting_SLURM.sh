#!/bin/bash
## Supply the Wochenende bam.txt input as arg1, bash run_Wochenende_reporting_SLURM.sh in.bam.txt

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 40000

# set run on bigmem node only
#SBATCH --cpus-per-task 4

# share node
#SBATCH --share

# set max wallclock time
# SBATCH --time=4:00:00

# set name of job
#SBATCH --job-name=Wo_report

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<mailaddress>

echo "Input bam: " $1
bamtxt=$1

# Source miniconda3
. /mnt/ngsnfs/tools/miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate wochenende


# read reference path from ref.tmp file in the reporting (current) directory 
ref=$(<ref.tmp)
echo "$ref"



# Run script
srun python3 basic_reporting.py --input_file $bamtxt --reference $ref --sequencer illumina --output_name $bamtxt

#conda activate wochenende
#python3 basic_reporting.py --input_file tmp_R1.ndp.lc.trm.s.mq30.01mm.dup.bam.txt --refseq_file /lager2/rcug/seqres/metagenref/2016_06_PPKC_metagenome_test_1p_spec_change_cln.fa --sequencer illumina --sample_name test

