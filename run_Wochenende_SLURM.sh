#!/bin/bash
## Supply the FASTA input as arg1, bash run_Wochenende_slurm.sh in.fastq

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 40000

# set run on bigmem node only
#SBATCH --cpus-per-task 12

# share node
#SBATCH --share

# set max wallclock time
# SBATCH --time=47:00:00

# set name of job
#SBATCH --job-name=Wochenende

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<mail_address_for_slurm>

echo "Input file: " $1
fastq=$1

# Add miniconda3 to PATH
PATH=/mnt/ngsnfs/tools/miniconda3/bin:$PATH

# Activate env on cluster node
source activate wochenende


### IMPORTANT
# Remember to check specified a) refseq b) threads c) adapters
##########

# Run script - Paired end reads R2 will be calculated by replacing R1 with R2
#python3 run_Wochenende.py --metagenome 2016_06_1p_spec_corrected --threads 12  --readType SE --debug $fastq
python3 run_Wochenende.py --metagenome 2016_06_1p_spec_corrected --threads 12 --aligner bwamem --no_abra --mq30 --readType SE --debug $fastq
#python3 run_Wochenende.py --metagenome PA14 --threads 36  --debug $fastq
#python3 run_Wochenende.py --metagenome hg19 --threads 55 --longread --debug $fastq
#python3 run_Wochenende.py --metagenome mm10 --threads 36 --readType PE --debug $fastq
#python3 run_Wochenende.py --metagenome GRCh38-mito --threads 16 --readType PE $fastq


