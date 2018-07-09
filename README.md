# Wochenende - A whole Genome/Metagenome Sequencing Alignment Pipeline

## Usage 

### SLURM usage

1. Copy all the run_Wochenende* files to your directory with your FASTQ files
`cp /path/to/wochenende/run_Wochenende* .`
2. Adjust settings in the script
`nano run_Wochenende_slurm.sh`
3. Run the pipeline
`sbatch run_Wochenende_SLURM.sh x.fastq`

### General usage

```
Wochenende - Whole Genome/Metagenome Sequencing Alignment Pipeline
Wochenende was created by Dr. Colin Davenport and Tobias Scheithauer
version: 1.0

usage: run_Wochenende.py [-h] [--readType {PE,SE}]
                         [--metagenome {2016_06_1p_genus,2016_06_1p_spec,hg19,GRCh38-45GB,GRCh38-noalt,GRCh38-mito,mm10,rn6,zf10,PA14}]
                         [--aligner ALIGNER] [--threads THREADS] [--fastp]
                         [--debug] [--longread] [--no_duplicate_removal]
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
```

## Installation

We recommend using [Bioconda](https://bioconda.github.io/) for installation of the tools used by our pipeline.

1. Clone or download the repository to your local machine.
`git clone https://github.com/MHH-RCUG/wochenende.git`
OR
`wget https://github.com/MHH-RCUG/wochenende/archive/master.zip`
2. Create a conda environment for the pipeline.
`conda create -n wochenende -c conda-forge -c bioconda bwa trimmomatic prinseq samtools=1.8 ncurses r-base64 sambamba fastuniq fastqc ea-utils bbmap`
3. Install all the other tools. (tools marked with an asterik are not necessarily needed)
   - afterqc *
   - fastp *
   - perldup
   - ABRA2
4. Edit the configuration section of run_Wochenende.py to set the paths to the tools and reference sequences.
5. Remember to activate the conda environment before trying to use the pipeline.

## List of Tools

- [ABRA2](https://github.com/mozack/abra2)
- [AfterQC](https://github.com/OpenGene/AfterQC)
- [BWA](https://github.com/lh3/bwa)
- [fastp](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [ea-utils](https://github.com/ExpressionAnalysis/ea-utils)
- [FastUniq](https://sourceforge.net/projects/fastuniq/)
- [Minimap2](https://github.com/lh3/minimap2)
- [perldup](https://github.com/richardmleggett/scripts/blob/master/remove_pcr_duplicates.pl) Already copied to dependencies folder. Developed by [Richard Leggett](https://github.com/richardmleggett).
- [PRINSEQ](http://prinseq.sourceforge.net/)
- [sambamba](https://github.com/biod/sambamba)
- [samtools](https://github.com/samtools/samtools)
- [trimmomatic](https://github.com/timflutre/trimmomatic)

# Old Version (internal only)

PPKC: A whole genome/metagenome alignment sequencing analysis pipeline
  * Dr. Colin Davenport
  * Tobias Scheithauer

Hannover Medical School, 2017-2018

See documentation on Wiki:
http://172.17.189.147/doku.php?id=install_ppkc_pipeline
