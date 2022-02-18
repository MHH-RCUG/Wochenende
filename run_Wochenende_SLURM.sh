#!/bin/bash
## A SLURM sbatch script which is part of the Wochenende Pipeline https://github.com/MHH-RCUG/Wochenende
## Supply the FASTQ read input as arg1, bash run_Wochenende_slurm.sh in_R1.fastq
## Reference fasta sequences and adapters should be configured in run_Wochenende.py

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 50000

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
#--aligner minimap2short   (for Illumina reads)
#--aligner minimap2long    (for nanopore reads)
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


# Set variables
cpus=12
#with SLURM srun, default. Change in config.yaml

# Setup SLURM using data parsed from config.yaml
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)
# Setup job scheduler
# use SLURM job scheduler (yes, no)
if [[ "${USE_CUSTOM_SCHED}" == "yes" ]]; then
    #echo USE_CUSTOM_SCHED set"
    scheduler=$CUSTOM_SCHED_CUSTOM_PARAMS
fi
if [[ "${USE_SLURM}" == "yes" ]]; then
    #echo USE_SLURM set"
    scheduler=$SLURM_CUSTOM_PARAMS
fi


# Run script - Paired end reads R2 will be calculated by replacing R1 with R2
# Uncomment/adapt the only line you want to run

#2021_12, minor update for 2021_10 ref. 
$scheduler python3 run_Wochenende.py --metagenome 2021_12_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_12_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 5 --readType PE --debug --no_prinseq --force_restart $fastq

#2021_10. Minor update of 2021_09. Streptococcus agalactiae double removed.
#$scheduler python3 run_Wochenende.py --metagenome 2021_10_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_10_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq


#2021_09. Update of 2021_02 ref. Archaea, Gemella, better fungi, includes common cold viruses. 
#$scheduler python3 run_Wochenende.py --metagenome 2021_09_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_09_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq

# 2021_02 reference update of 2020_03, with fungi, can better detect C. diff, and B. subtilis. No UNVERIF sp like one Achromobacter
# blacklister masked
#$scheduler python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.p --metagenome 2021_02_meta_fungi_human_masked --threads $cpus --aligner bwamem  --nextera --trim_galore --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_masked --threads $cpus --aligner minimap2short --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_masked --threads $cpus --aligner minimap2short --no_abra --mq30 --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_masked --threads $cpus --aligner minimap2long --longread --no_abra --mq30 --remove_mismatching 250 --readType SE --debug --no_prinseq --force_restart $fastq
# not blacklister masked
#$scheduler python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_unmasked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_unmasked --threads $cpus --aligner bwamem  --nextera --trim_galore --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_02_meta_fungi_human_unmasked --threads $cpus --aligner minimap2long --longread --no_abra --mq30 --remove_mismatching 250 --readType SE --debug --no_prinseq --force_restart $fastq


# 2021_08 test: 2021_08_meta_fungi_human_masked. Warning - poor representation of P. aeruginosa and other common pathogens due to genomic crowding
#$scheduler python3 run_Wochenende.py --metagenome 2021_08_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq

# 2020_09 massive reference - requires ca 120GB RAM per process, set #SBATCH --mem 450000 and #SBATCH --cpus-per-task 24
#$scheduler python3 run_Wochenende.py --metagenome 2020_09_massiveref_human --threads $cpus --aligner bwamem --no_abra        --remove_mismatching 3 --readType SE --no_prinseq --debug --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2020_09_massiveref_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --no_prinseq --debug --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2020_09_massiveref_human --threads $cpus --aligner bwamem --no_abra --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2020_09_massiveref_human --threads $cpus --aligner minimap2short --no_abra --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq


#Alignerboost test
#$scheduler python3 run_Wochenende.py --runAlignerBoost --no_fastqc --metagenome 2021_02_meta_fungi_human_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq


# 2020_05 reference (blacklister-masked version of 2020_03 below)
#$scheduler python3 run_Wochenende.py --metagenome 2020_05_meta_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2020_05_meta_human --threads $cpus --aligner bwamem --nextera --trim_galore --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2020_05_meta_human --threads $cpus --aligner minimap2long --longread --no_abra --mq30 --remove_mismatching 250 --readType SE --debug --no_prinseq --force_restart $fastq




# 2020_03 reference (unmasked, Achromobacter problem and other masking issues, use 2020_05 above or 2021_02 instead)
#$scheduler python3 run_Wochenende.py --metagenome 2020_03_meta_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2020_03_meta_human --threads $cpus --aligner bwamem --nextera --trim_galore --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2020_03_meta_human --threads $cpus --aligner minimap2long --longread --no_abra --mq30 --remove_mismatching 250 --readType SE --debug --no_prinseq --force_restart $fastq




## 2019 10 October metagenomes
#$scheduler python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
# with fastp
#$scheduler python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --fastp --no_prinseq --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --force_restart $fastq
# longread with minimap2 aligner
#$scheduler python3 run_Wochenende.py --metagenome 2019_10_meta_human --threads $cpus --longread --no_prinseq --aligner minimap2long --mq30  --remove_mismatching 250 --readType SE --debug --force_restart $fastq


## 2019 10 October metagenomes with Univec contamination
#$scheduler python3 run_Wochenende.py --metagenome 2019_10_meta_human_univec --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --no_prinseq --debug --force_restart $fastq

# 2021_11 mouse
#$scheduler python3 run_Wochenende.py --metagenome 2021_11_meta_fungi_mouse_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_11_meta_fungi_mouse_masked --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq


# 2021_07 mouse - deprecated
#$scheduler python3 run_Wochenende.py --metagenome 2021_07_meta_fungi_mouse_masked --threads $cpus --aligner bwamem --no_abra --mq20 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_07_meta_fungi_mouse_masked --threads $cpus --aligner bwamem --no_abra --mq20 --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq


## 2019 01 January metagenomes
#$scheduler python3 run_Wochenende.py --metagenome 2019_01_meta --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2019_01_meta_mouse --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2019_01_meta_mouse --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType PE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_ASF_OMM --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_ASF --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2019_01_meta_mouse_OMM --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq

## 2020 03 NCI viruses only - with mq30. Use with unmapped reads only after removing human + bact
#$scheduler python3 run_Wochenende.py --metagenome nci_viruses --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
## 2020 03 NCI viruses only - no mq30.  Use with unmapped reads only after removing human + bact
#$scheduler python3 run_Wochenende.py --metagenome nci_viruses --threads $cpus --aligner bwamem --no_abra --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq
# EZV viruses
#$scheduler python3 run_Wochenende.py --metagenome ezv_viruses --threads $cpus --aligner bwamem --no_abra --readType SE --no_prinseq --debug --no_prinseq --force_restart $fastq


# Univec added for exclusion of contamination
#$scheduler python3 run_Wochenende.py --metagenome 2019_10_meta_human_univec --threads $cpus --aligner bwamem --no_abra --mq30 --remove_mismatching 3 --readType SE --debug --no_prinseq --force_restart $fastq

# Seqins from Anaquini metasequin_sequences_3.0.fa. If present, final eg calmd.bams probably have MB size files, but check alignments.
#$scheduler python3 run_Wochenende.py --metagenome seqins_v3 --threads $cpus --readType SE --debug --force_restart --no_abra --no_fastqc --mq30 --remove_mismatching 3 --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome 2021_12_meta_fungi_human_masked_seqins_v3 --threads $cpus --readType SE --debug --force_restart --no_abra --no_fastqc --mq30 --remove_mismatching 3 --no_prinseq $fastq

#plasmids
#$scheduler python3 run_Wochenende.py --metagenome ssplasmid1 --threads $cpus --aligner minimap2long --longread --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome ssplasmid2 --threads $cpus --aligner minimap2long --longread --debug --force_restart --no_prinseq $fastq

# Genomes - remember higher --remove_mismatching and also maybe --nextera
#$scheduler python3 run_Wochenende.py --metagenome strept_halo --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome k_variicola --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome k_oxytoca --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome clost_bot --threads $cpus --readType PE  --aligner bwamem --mq30 --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome clost_bot_e --threads $cpus --readType PE  --aligner bwamem --mq30 --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome clost_perf --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome clost_diff --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome citro_freundii --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome ecoli --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome PA14 --threads $cpus --readType SE --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome hg19 --threads $cpus --longread --remove_mismatching 250 --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome hg19 --threads $cpus --readType PE --no_duplicate_removal --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome GRCh37 --threads $cpus --readType PE --aligner bwamem --no_duplicate_removal --no_abra --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome mm10 --threads $cpus --readType PE --aligner bwamem --no_abra --remove_mismatching 4 --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome mm10 --threads $cpus --no_duplicate_removal --no_abra --readType SE --longread --aligner minimap2long --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome rat_1AR1_ont --threads $cpus --no_duplicate_removal --readType SE --aligner bwa --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome ss11 --threads $cpus --readType PE --no_duplicate_removal --no_abra --debug --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome GRCh38-mito --threads $cpus --readType PE --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome GRCh38-noalt --threads $cpus --readType SE --aligner minimap2short --no_duplicate_removal --no_abra --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome GRCh38-noalt --threads $cpus --readType SE --aligner minimap2long --longread --no_duplicate_removal --no_abra --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome GRCh38-noalt --threads $cpus --readType SE --aligner bwamem --no_duplicate_removal --no_abra --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome GRCh38-noalt --threads $cpus --readType PE --aligner bwamem --no_duplicate_removal --no_abra --remove_mismatching 3 --mq30 --force_restart --no_prinseq $fastq
#$scheduler python3 run_Wochenende.py --metagenome T2T_v1_1 --threads $cpus --readType SE --aligner minimap2long --longread --no_duplicate_removal --no_abra --force_restart --no_prinseq $fastq

# Test wochenende pipeline. Run tests with sbatch run_Wochenende_SLURM.sh testdb/reads_R1.fastq
#$scheduler python3 run_Wochenende.py --metagenome testdb --threads $cpus --testWochenende --aligner bwamem --mq30 --remove_mismatching 3 --readType SE --debug --force_restart $fastq

# Test wochenemde with fastp - the auto tests here will fail here since fastp not trm is in the filename
#$scheduler python3 run_Wochenende.py --metagenome testdb --fastp --threads $cpus --testWochenende --aligner bwamem --mq30 --remove_mismatching 3 --readType SE --debug --force_restart $fastq


