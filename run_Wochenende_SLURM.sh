#!/bin/bash
## Supply the FASTA input as arg1, bash run_Wochenende_slurm.sh in.fastq

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 40000

# set run x cpus
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
. /mnt/ngsnfs/tools/miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate wochenende


### IMPORTANT
# Remember to check specified a) refseq b) threads c) adapters
##########

# All options
# --threads 12
#--longread
#--aligner bwamem
#--aligner minimap2
#--no_abra
#--mq30
#--readType SE
#--readType PE
#--debug
#--force_restart
#--remove_mismatching
#--no_duplicate_removal
#--no_prinseq
#--no_fastqc
#--testWochenende
#--fastp - fastp is recommended as a trimmer for SOLiD data which can auto find adapters

cpus=12

# Run script - Paired end reads R2 will be calculated by replacing R1 with R2

# Test wochenende pipeline. Run tests with sbatch run_Wochenende_SLURM.sh testdb/reads_R1.fastq
#python3 run_Wochenende.py --metagenome testdb --threads $cpus --testWochenende --aligner bwamem --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq

# Test wochenemde with fastp - the auto tests here will fail here since fastp not trm is in the filename
#python3 run_Wochenende.py --metagenome testdb --fastp --threads $cpus --testWochenende --aligner bwamem --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq

## 2019 10 October metagenomes
#python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
# without prinseq
#python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --no_fastqc --no_prinseq --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq


## 2019 01 January metagenomes
#python3 run_Wochenende.py --metagenome 2019_01_meta --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_ASF_OMM --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_ASF --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_OMM --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq

## 2020 03 NCI viruses only - with mq30
#python3 run_Wochenende.py --metagenome nci_viruses --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
## 2020 03 NCI viruses only - no mq30
python3 run_Wochenende.py --metagenome nci_viruses --threads $cpus --aligner bwamem --no_abra --remove_mismatching --readType SE --debug --force_restart $fastq

#2016 2016_06_PPKC_metagenome_test_1p_spec_change
#python3 run_Wochenende.py --metagenome 2016_06_1p_spec_corrected --threads $cpus  --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2016_06_1p_spec_corrected --threads $cpus --aligner bwamem --no_abra --mq30 --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2016_06_1p_spec_corrected --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq

# Genomes
#python3 run_Wochenende.py --metagenome PA14 --threads $cpus  --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome hg19 --threads $cpus --longread --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome hg19 --threads $cpus --readType PE --no_duplicate_removal --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome GRCh37 --threads $cpus --readType PE --aligner bwamem --no_duplicate_removal --no_abra --force_restart $fastq
#python3 run_Wochenende.py --metagenome mm10 --threads $cpus --readType PE --aligner bwamem --force_restart $fastq
#python3 run_Wochenende.py --metagenome mm10 --threads $cpus --no_duplicate_removal --no_abra --readType SE --longread --aligner minimap2 --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome rat_1AR1_ont --threads $cpus --no_duplicate_removal --readType SE --aligner bwa --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome ss11 --threads $cpus --readType PE --no_duplicate_removal --no_abra --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome GRCh38-mito --threads $cpus --readType PE --force_restart $fastq
#python3 run_Wochenende.py --metagenome GRCh38-noalt --threads $cpus --readType SE --aligner minimap2 --longread --no_duplicate_removal --no_abra --force_restart $fastq
#python3 run_Wochenende.py --metagenome GRCh38-noalt --threads $cpus --readType SE --aligner bwamem --no_duplicate_removal --no_abra --force_restart $fastq
#python3 run_Wochenende.py --metagenome GRCh38-noalt --threads $cpus --readType PE --aligner bwamem --no_duplicate_removal --no_abra --remove_mismatching --mq30 --force_restart $fastq
