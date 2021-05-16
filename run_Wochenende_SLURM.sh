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

# Add miniconda3 to PATH
. /mnt/ngsnfs/tools/miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate wochenende >> /dev/null

echo "Input file: " $1
fastq=$1


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
#--mq20     - remove reads with a mapping quality of less than 20. Less stringent than MQ30, required for raspir https://github.com/mmpust/raspir
#--mq30     - remove reads with a mapping quality of less than 30
#--runAlignerBoost - run the AlignerBoost software to adjust mapping quality.
#--readType SE - single ended reads
#--readType PE - paired end reads
#--debug
#--force_restart
#--remove_mismatching 3 (remove those reads with 3 or more mismatches) # questionable for very long reads
#--remove_mismatching 250 (remove those reads with 250 or more mismatches) # for long reads set to ~10% of median read length.    
#--no_duplicate_removal  - do not remove duplicate reads
#--no_prinseq   - do not filter out low complexity initial reads using prinseq (default in this file after 2020_11)
#--no_fastqc
#--testWochenende - runs the test scripts with test reads vs a testDB and checks if all seems well.
#--fastp - fastp is recommended as an alternative trimmer to Trimmomatic if you are having adapter problems
#--trim_galore - trim_galore is the best adapter trimmer for nextera reads

cpus=12

# Run script - Paired end reads R2 will be calculated by replacing R1 with R2
# Uncomment/adapt the only line you want to run




# 2020_09 massive reference - requires ca 120GB RAM per process, set #SBATCH --mem 450000 and #SBATCH --cpus-per-task 24
#python3 run_Wochenende.py --metagenome 2020_09_massiveref_human --threads $cpus --aligner bwamem --no_abra        --remove_mismatching 3 --readType SE --no_prinseq --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2020_09_massiveref_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --no_prinseq --debug --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2020_09_massiveref_human --threads $cpus --aligner bwamem --no_abra --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2020_09_massiveref_human --threads $cpus --aligner minimap2 --no_abra --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq

# 2021_02 reference update of 2020_03, with fungi, can better detect C. diff, and B. subtilis. No UNVERIF sp like one Achromobacter
# blacklister masked
python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.p --metagenome 2021_02_meta_fungi_human_masked --threads $cpus --aligner bwamem  --nextera --trim_galore --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_masked --threads $cpus --aligner minimap2 --longread --no_abra --mq30 --remove_mismatching 250 --readType SE --debug --no_prinseq --force_restart $fastq
# not blacklister masked
#python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_unmasked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_unmasked --threads $cpus --aligner bwamem  --nextera --trim_galore --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_unmasked --threads $cpus --aligner minimap2 --longread --no_abra --mq30 --remove_mismatching 250 --readType SE --debug --no_prinseq --force_restart $fastq

#Alignerboost test
#python3 run_Wochenende.py --runAlignerBoost --no_fastqc --metagenome 2021_02_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq


# 2020_05 reference (blacklister-masked version of 2020_03 below)
#python3 run_Wochenende.py --metagenome 2020_05_meta_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2020_05_meta_human --threads $cpus --aligner bwamem --nextera --trim_galore --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2020_05_meta_human --threads $cpus --aligner minimap2 --longread --no_abra --mq30 --remove_mismatching 250 --readType SE --debug --no_prinseq --force_restart $fastq



# 2020_03 reference (unmasked, Achromobacter problem and other masking issues, use 2020_05 above or 2021_02 instead)
#python3 run_Wochenende.py --metagenome 2020_03_meta_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2020_03_meta_human --threads $cpus --aligner bwamem --nextera --trim_galore --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2020_03_meta_human --threads $cpus --aligner minimap2 --longread --no_abra --mq30 --remove_mismatching 250 --readType SE --debug --no_prinseq --force_restart $fastq



## 2019 10 October metagenomes
#python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
# with fastp
#python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --fastp --no_prinseq --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --force_restart $fastq
# longread with minimap2 aligner
#python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --longread --no_prinseq --aligner minimap2 --mq30  --remove_mismatching 250 --readType SE --debug --force_restart $fastq

## 2019 10 October metagenomes with Univec contamination
#python3 run_Wochenende.py --metagenome 2019_10_meta_human_univec --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --no_prinseq --debug --force_restart $fastq


## 2019 01 January metagenomes
#python3 run_Wochenende.py --metagenome 2019_01_meta --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_ASF_OMM --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_ASF --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_OMM --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq

## 2020 03 NCI viruses only - with mq30. Use with unmapped reads only after removing human + bact
#python3 run_Wochenende.py --metagenome nci_viruses --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
## 2020 03 NCI viruses only - no mq30.  Use with unmapped reads only after removing human + bact
#python3 run_Wochenende.py --metagenome nci_viruses --threads $cpus --aligner bwamem --no_abra --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
# EZV viruses
#python3 run_Wochenende.py --metagenome ezv_viruses --threads $cpus --aligner bwamem --no_abra --readType SE --no_prinseq --debug --no_prinseq --force_restart $fastq


# Univec added for exclusion of contamination
#python3 run_Wochenende.py --metagenome 2019_10_meta_human_univec --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq


# Genomes
#python3 run_Wochenende.py --metagenome strept_halo --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome k_variicola --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome k_oxytoca --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome clost_bot --threads $cpus --readType PE  --aligner bwamem --mq30 --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome clost_bot_e --threads $cpus --readType PE  --aligner bwamem --mq30 --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome clost_perf --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome clost_diff --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome citro_freundii --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome ecoli --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome PA14 --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome hg19 --threads $cpus --longread --remove_mismatching 250 --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome hg19 --threads $cpus --readType PE --no_duplicate_removal --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome GRCh37 --threads $cpus --readType PE --aligner bwamem --no_duplicate_removal --no_abra --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome mm10 --threads $cpus --readType PE --aligner bwamem --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome mm10 --threads $cpus --no_duplicate_removal --no_abra --readType SE --longread --aligner minimap2 --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome rat_1AR1_ont --threads $cpus --no_duplicate_removal --readType SE --aligner bwa --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome ss11 --threads $cpus --readType PE --no_duplicate_removal --no_abra --debug --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome GRCh38-mito --threads $cpus --readType PE --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome GRCh38-noalt --threads $cpus --readType SE --aligner minimap2 --longread --no_duplicate_removal --no_abra --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome GRCh38-noalt --threads $cpus --readType SE --aligner bwamem --no_duplicate_removal --no_abra --force_restart --no_prinseq $fastq
#python3 run_Wochenende.py --metagenome GRCh38-noalt --threads $cpus --readType PE --aligner bwamem --no_duplicate_removal --no_abra --remove_mismatching 3 --mq30 --force_restart --no_prinseq $fastq

# Test wochenende pipeline. Run tests with sbatch run_Wochenende_SLURM.sh testdb/reads_R1.fastq
#python3 run_Wochenende.py --metagenome testdb --threads $cpus --testWochenende --aligner bwamem --mq30 --remove_mismatching 3 --readType SE --debug --force_restart $fastq

# Test wochenemde with fastp - the auto tests here will fail here since fastp not trm is in the filename
#python3 run_Wochenende.py --metagenome testdb --fastp --threads $cpus --testWochenende --aligner bwamem --mq30 --remove_mismatching 3 --readType SE --debug --force_restart $fastq


