#!/bin/bash
## Supply the FASTA input as arg1, bash run_Wochenende_slurm.sh in.fastq

# set partition
#SBATCH -p old

# set run on x MB node only
#SBATCH --mem 20000

# set run on bigmem node only
#SBATCH --cpus-per-task 1

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

# Add miniconda3 to PATH
PATH=/mnt/ngsnfs/tools/miniconda3/bin:$PATH

# Activate env on cluster node
conda activate wochenende

ref2016=/lager2/rcug/seqres/metagenref/bwa/2016_06_PPKC_metagenome_test_1p_spec_change_cln.fa

# Run script
python3 basic_reporting.py --input_file $bamtxt --refseq_file $ref2016 --sequencer illumina --sample_name $bamtxt

#conda activate wochenende
#python3 basic_reporting.py --input_file tmp_R1.ndp.lc.trm.s.mq30.01mm.dup.bam.txt --refseq_file /lager2/rcug/seqres/metagenref/2016_06_PPKC_metagenome_test_1p_spec_change_cln.fa --sequencer illumina --sample_name test

