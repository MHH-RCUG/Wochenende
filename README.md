Wochenende - Whole Genome/Metagenome Sequencing Alignment Pipeline
Wochenende was created by Dr. Colin Davenport and Tobias Scheithauer
version: 1.0

usage: run_Wochenende.py [-h] [--readType {PE,SE}]
                         [--metagenome {2016_06_1p_genus,2016_06_1p_spec,hg19,GRCh38-45GB,GRCh38-noalt,GRCh38-mito,mm10,rn6,zf10,PA14}]
                         [--aligner ALIGNER] [--threads THREADS] [--fastp]
                         [--debug] [--longread] [--no_duplicate_removal]
                         [--force_restart]
                         fastq

positional arguments:
  fastq                 _R1.fastq Input read1 fastq file

optional arguments:
  -h, --help            show this help message and exit
  --readType {PE,SE}    Single end or paired end data, default= SE
  --metagenome {2016_06_1p_genus,2016_06_1p_spec,hg19,GRCh38-45GB,GRCh38-noalt,GRCh38-mito,mm10,rn6,zf10,PA14}
                        Metagenome to use
  --aligner ALIGNER     Aligner to use, default= bwamem. Usage of minimap2
                        optimized for ONT data only.
  --threads THREADS     Number of cores, default = 16
  --fastp               Use fastp instead of fastqc and trimmomatic
  --debug               Report all files
  --longread            Only do steps relevant for long PacBio/ONT reads eg.
                        no trimming, alignment & bam conversion
  --no_duplicate_removal
                        Skips steps for duplicate removal. Recommended for
                        amplicon sequencing.
  --force_restart       Force restart, without regard to existing progress

We recommend using bioconda for the installation of the tools. Remember to run
'source activate <environment name>' before you start if you are using
bioconda. Details about the installation are available on
https://github.com/MHH-RCUG/Wochenende#Installation