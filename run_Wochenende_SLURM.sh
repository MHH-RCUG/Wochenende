#!/bin/bash
## A SLURM sbatch script which is part of the Wochenende Pipeline https://github.com/MHH-RCUG/Wochenende
## Supply the FASTQ read input as arg1, bash run_Wochenende_slurm.sh in_R1.fastq
## Reference fasta sequences and adapters should be configured in run_Wochenende.py

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 40000

# set run x cpus
#SBATCH --cpus-per-task 12

# set max wallclock time
#SBATCH --time=47:00:00

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
conda activate wochenende >> /dev/null


### IMPORTANT
# Remember to check specified a) refseq b) threads c) adapters
##########

# All options
# --threads 12
#--longread (implies aligner is not bwa-mem, no prinseq, no duplicate removal)
#--aligner bwamem
#--aligner minimap2
#--aligner ngmlr
#--nextera  - remove Nextera adapters with Trimmomatic, not default Ultra II / Truseq adapters
#--no_abra  - no read realignment
#--mq30     - remove reads with a mapping quality of less than 30
#--readType SE
#--readType PE
#--debug
#--force_restart
#--remove_mismatching
#--no_duplicate_removal
#--no_prinseq   - do not filter out initial reads using prinseq
#--no_fastqc
#--testWochenende - runs the test scripts with test reads vs a testDB and checks if all seems well.
#--fastp - fastp is recommended as an alternative trimmer to Trimmomatic if you are having adapter problems

cpus=12

# Run script - Paired end reads R2 will be calculated by replacing R1 with R2
# Uncomment/adapt the only line you want to run





# 2020_03 reference
#python3 run_Wochenende.py --metagenome 2020_03_meta_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2020_03_meta_human --threads $cpus --aligner bwamem --nextera --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2020_03_meta_human --threads $cpus --aligner minimap2 --longread --no_abra --mq30 --readType SE --debug --force_restart $fastq
python3 run_Wochenende.py --metagenome 2020_03_meta_human --threads $cpus --aligner ngmlr --longread --no_abra --mq30 --readType SE --debug --force_restart $fastq



## 2019 10 October metagenomes
#python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
# without prinseq
#python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --no_fastqc --no_prinseq --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
# with fastp
#python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --fastp --no_prinseq --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
# longread with ngmlr aligner
#python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --longread --no_prinseq --aligner ngmlr --mq30 --readType SE --debug --force_restart $fastq
# longread with minimap2 aligner
#python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --longread --no_prinseq --aligner minimap2 --mq30 --readType SE --debug --force_restart $fastq

## 2019 10 October metagenomes with Univec contamination
#python3 run_Wochenende.py --metagenome 2019_10_meta_human_univec --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq


## 2019 01 January metagenomes
#python3 run_Wochenende.py --metagenome 2019_01_meta --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_ASF_OMM --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_ASF --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_OMM --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq

## 2020 03 NCI viruses only - with mq30. Use with unmapped reads only after removing human + bact
#python3 run_Wochenende.py --metagenome nci_viruses --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq
## 2020 03 NCI viruses only - no mq30.  Use with unmapped reads only after removing human + bact
#python3 run_Wochenende.py --metagenome nci_viruses --threads $cpus --aligner bwamem --no_abra --remove_mismatching --readType SE --debug --force_restart $fastq
# EZV viruses
#python3 run_Wochenende.py --metagenome ezv_viruses --threads $cpus --aligner bwamem --no_abra --readType SE --no_prinseq --debug --force_restart $fastq


# Univec added for exclusion of contamination
#python3 run_Wochenende.py --metagenome 2019_10_meta_human_univec --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq



#2016 2016_06_PPKC_metagenome_test_1p_spec_change
#python3 run_Wochenende.py --metagenome 2016_06_1p_spec_corrected --threads $cpus  --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2016_06_1p_spec_corrected --threads $cpus --aligner bwamem --no_abra --mq30 --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2016_06_1p_spec_corrected --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq

# Genomes
#python3 run_Wochenende.py --metagenome PA14 --threads $cpus --readType SE --debug --force_restart $fastq
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

# Test wochenende pipeline. Run tests with sbatch run_Wochenende_SLURM.sh testdb/reads_R1.fastq
#python3 run_Wochenende.py --metagenome testdb --threads $cpus --testWochenende --aligner bwamem --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq

# Test wochenemde with fastp - the auto tests here will fail here since fastp not trm is in the filename
#python3 run_Wochenende.py --metagenome testdb --fastp --threads $cpus --testWochenende --aligner bwamem --mq30 --remove_mismatching --readType SE --debug --force_restart $fastq

# Test new ngmlr aligner
#python3 run_Wochenende.py --metagenome testdb --longread --threads $cpus --aligner ngmlr --readType SE --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome testdb --longread --threads $cpus --aligner minimap2 --readType SE --debug --force_restart $fastq

