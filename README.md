# Wochenende - A whole Genome/Metagenome Sequencing Alignment Pipeline

Wochenende runs alignment of either single ended or paired end short reads against a reference sequence. It is simple (python script), portable and has many optional steps. 

Features include (see programs listed below at the bottom of this page)
- QC (Fastqc)
- pre alignment duplicate removal (perldup)
- pre alignment poor sequence removal (Prinseq)
- trimming (Trimmomatic, fastp)
- alignment (bwa mem, minimap2)
- SAM-> BAM conversion (samtools, sambamba)
- Post-alignment duplicate removal (Picard)
- Realignment (Abra2)
- MD tag marking (Samtools)




## Usage 

### SLURM usage

1. Copy all the run_Wochenende* files to your directory with your FASTQ files
`cp /path/to/wochenende/run_Wochenende* .`
2. Adjust settings in the script
`nano run_Wochenende_SLURM.sh`
3. Run the pipeline using SLURM
`sbatch run_Wochenende_SLURM.sh x.fastq`
4. Optional reporting step to normalize the extracted read counts  (see reporting below)
`sbatch run_Wochenende_reporting_SLURM.sh`

### General usage

```
Wochenende - Whole Genome/Metagenome Sequencing Alignment Pipeline
Wochenende was created by Dr. Colin Davenport and Tobias Scheithauer
version: 1.1


usage: run_Wochenende.py [-h] [--aligner {bwamem,minimap2}]
                         [--readType {PE,SE}]
                         [--metagenome {zf10,GRCh38-mito,mm10,rn6,2016_06_1p_spec,GRCh38-noalt,PA14,2016_06_1p_spec_corrected,2016_06_1p_genus,GRCh38-45GB,hg19}]
                         [--threads THREADS] [--fastp] [--debug] [--longread]
                         [--no_duplicate_removal] [--no_abra] [--mq30]
                         [--remove_mismatching] [--force_restart]
                         fastq

positional arguments:
  fastq                 _R1.fastq Input read1 fastq file

optional arguments:
  -h, --help            show this help message and exit
  --aligner {bwamem,minimap2}
                        Aligner to use, either bwamem or minimap2. Usage of
                        minimap2 optimized for ONT data only.
  --readType {PE,SE}    Single end or paired end data
  --metagenome {zf10,GRCh38-mito,mm10,rn6,2016_06_1p_spec,GRCh38-noalt,PA14,2016_06_1p_spec_corrected,2016_06_1p_genus,GRCh38-45GB,hg19}
                        Meta/genome reference to use
  --threads THREADS     Number of cores, default = 16
  --fastp               Use fastp instead of fastqc and trimmomatic
  --debug               Report all files
  --longread            Only do steps relevant for long PacBio/ONT reads eg.
                        no trimming, alignment & bam conversion
  --no_duplicate_removal
                        Skips steps for duplicate removal. Recommended for
                        amplicon sequencing.
  --no_abra             Skips steps for Abra realignment. Recommended for
                        metagenome and amplicon analysis.
  --mq30                Remove reads with mapping quality less than 30.
                        Recommended for metagenome and amplicon analysis.
  --remove_mismatching  Remove reads with 2 or more mismatches (via the NM bam
                        tag)
  --force_restart       Force restart, without regard to existing progress

We recommend using bioconda for the installation of the tools. Remember to run
'source activate <environment name>' before you start if you are using
bioconda. Details about the installation are available on
https://github.com/MHH-RCUG/Wochenende#installation
```

## Installation

We recommend using [Bioconda](https://bioconda.github.io/) for installation of the tools used by our pipeline.

1. Clone or download the repository to your local machine.
`git clone https://github.com/MHH-RCUG/wochenende.git`
OR
`wget https://github.com/MHH-RCUG/wochenende/archive/master.zip`
2. Create a conda environment for the pipeline.
`conda create -n wochenende -c conda-forge -c bioconda bwa trimmomatic prinseq samtools=1.8 ncurses r-base64 sambamba fastuniq fastqc ea-utils bbmap fastp minimap2 bamtools`
3. Install all the other tools.
   - ABRA2
4. Edit the configuration section of `run_Wochenende.py` to set the paths to the tools and reference sequences.
5. Activate the conda environment before running the pipeline.
`source activate wochenende`

## List of Tools used or optional in the pipeline

- [ABRA2](https://github.com/mozack/abra2)
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
- [bamtools](https://github.com/pezmaster31/bamtools)


### Running the metagenomic reporting scripts

This tool reports length, GC content of the sequence, read counts attributed to the species and various normalized read count parameters. 
Normalizations are for:

a) reads normalized to the idealized length of a bacterial chromosome (normalization to 1 million base pairs)

b) total reads in the sequencing library (normalization to 1 million reads)

c) the above two normalizations combined (RPMM)


See the subfolder reporting in the repository.

```
conda activate 
conda activate wochenende
python3 basic_reporting.py --input_file tmp_R1.ndp.lc.trm.s.mq30.01mm.dup.bam.txt --refseq_file /lager2/rcug/seqres/metagenref/2016_06_PPKC_metagenome_test_1p_spec_change_cln.fa --sequencer illumina --sample_name test
```

### Wochenende output

Wochenende produces many output files, many of which are superseded by later output files and can be removed.

- MB_aero_S2_R1.fastq # Input file Read1. Note the form R1.fastq is required, R1_001.fastq will not work well.
- MB_aero_S2_R1.fastqprogress.tmp # Temporary file with pipeline stage progress
- MB_aero_S2_R1.trm.bam
- MB_aero_S2_R1.trm.fastq
- MB_aero_S2_R1.trm.s.bam
- MB_aero_S2_R1.trm.s.bam.bai
- MB_aero_S2_R1.trm.s.mq30.01mm.bam
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam.bai
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam.txt
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam.txt.filt.sort.csv
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam.txt.reporting.sorted.csv
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam.txt.reporting.unsorted.csv
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.calmd.bam
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.calmd.bam.bai
- MB_aero_S2_R1.trm.s.mq30.bam
- MB_aero_S2_R2.fastq
- MB_aero_S2_R2.trm.fastq



