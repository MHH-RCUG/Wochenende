# Wochenende - A whole Genome/Metagenome Sequencing Alignment Pipeline

Wochenende runs alignment of either single ended or paired end short reads against a reference sequence. It is simple (python script), portable and has many optional steps. 

Features include (see programs listed below at the bottom of this page)
- QC (Fastqc)
- pre alignment duplicate removal (perldup)
- pre alignment poor sequence removal (Prinseq - used for single ended reads only)
- trimming (Trimmomatic, fastp)
- alignment (bwa mem, minimap2)
- SAM-> BAM conversion (samtools, sambamba)
- Report % aligned reads (samtools)
- Output unmapped reads as fastq (samtools)  (from v1.4)
- Post-alignment duplicate removal (Picard)
- Realignment (Abra2)
- MD tag marking (Samtools)
- Normalization (to Reads per Million Reads etc, see Reporting below for details)
- Visualization (chromosome coverage, intended for bacteria in metagenomics projects) (from v1.4)




## Usage 

You can just run the pipeline as a normal Python3 script. However, we also offer a template for the job scheduler SLURM below. This template can also be used with Bash to run the commands at the bottom of the SLURM pipeline while ignoring any SLURM specific parameters.

### SLURM usage

1. Copy all the run_Wochenende* files to your directory with your FASTQ files
`cp /path/to/wochenende/run_Wochenende* .`
2. Adjust settings in the script
`nano run_Wochenende_SLURM.sh`
3. Run the pipeline using SLURM (the "_R1" is important)
`sbatch run_Wochenende_SLURM.sh sample_R1.fastq`
4. Optional reporting step to normalize the extracted read counts  (see reporting below)
`sbatch run_Wochenende_reporting_SLURM.sh`

### General usage

```
Wochenende - Whole Genome/Metagenome Sequencing Alignment Pipeline
Wochenende was created by Dr. Colin Davenport and Tobias Scheithauer
version: 1.4 - Feb 2020

usage: run_Wochenende.py [-h] [--aligner {bwamem,minimap2}]
                         [--readType {PE,SE}]
                         [--metagenome {2019_01_meta_mouse,2019_01_meta_mouse_ASF_OMM,2019_01_meta_mouse_OMM,2019_10_meta_human,ss11,testdb,2019_01_meta_mouse_ASF,GRCh38-45GB,GRCh37,mm10,GRCh38-noalt,hg19,PA14,2016_06_1p_spec_corrected,rat_1AR1_ont,2016_06_1p_genus,2019_01_meta,zf10,2016_06_1p_spec,GRCh38-mito,rn6}]
                         [--threads THREADS] [--fastp] [--debug] [--longread]
                         [--no_duplicate_removal] [--no_prinseq] [--no_fastqc]
                         [--no_abra] [--mq30] [--remove_mismatching]
                         [--force_restart] [--testWochenende]
                         fastq

positional arguments:
  fastq                 _R1.fastq Input read1 fastq file

optional arguments:
  -h, --help            show this help message and exit
  --aligner {bwamem,minimap2}
                        Aligner to use, either bwamem or minimap2. Usage of
                        minimap2 optimized for ONT data only.
  --readType {PE,SE}    Single end or paired end data
  --metagenome {2019_01_meta_mouse,2019_01_meta_mouse_ASF_OMM,2019_01_meta_mouse_OMM,2019_10_meta_human,ss11,testdb,2019_01_meta_mouse_ASF,GRCh38-45GB,GRCh37,mm10,GRCh38-noalt,hg19,PA14,2016_06_1p_spec_corrected,rat_1AR1_ont,2016_06_1p_genus,2019_01_meta,zf10,2016_06_1p_spec,GRCh38-mito,rn6}
                        Meta/genome reference to use
  --threads THREADS     Number of cores, default = 16
  --fastp               Use fastp instead of fastqc and trimmomatic
  --debug               Report all files
  --longread            Only do steps relevant for long PacBio/ONT reads eg.
                        no trimming, alignment & bam conversion
  --no_duplicate_removal
                        Skips steps for duplicate removal. Recommended for
                        amplicon sequencing.
  --no_prinseq          Skips prinseq step for low_complexity sequence
                        removal.
  --no_fastqc           Skips FastQC quality control step.
  --no_abra             Skips steps for Abra realignment. Recommended for
                        metagenome and amplicon analysis.
  --mq30                Remove reads with mapping quality less than 30.
                        Recommended for metagenome and amplicon analysis.
  --remove_mismatching  Remove reads with 2 or more mismatches (via the NM bam
                        tag)
  --force_restart       Force restart, without regard to existing progress
  --testWochenende      Run pipeline tests vs testdb, needs the subdirectory
                        testdb, default false

We recommend using bioconda for the installation of the tools. Remember to run
'conda activate <environment name>' before you start if you are using
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
`conda create -n wochenende -c conda-forge -c bioconda bwa trimmomatic prinseq samtools=1.8 ncurses r-base64 sambamba=0.6.6 fastuniq fastqc ea-utils bbmap fastp minimap2 bamtools`
3. Install all the other tools.
   - ABRA2
4. Edit the configuration section of `run_Wochenende.py` to set the paths to the tools and reference sequences.
5. (Wochenende_plot only). Install the python dependencies for visualization by pip. The works on Ubuntu 1604: `pip3 install --user numpy==1.17.4 pandas==0.24.2 matplotlib==3.0.3`
6. Activate the conda environment before running the pipeline.
`conda activate wochenende`
7. Optional: run the tests, see below.

## List of Tools used or optional in the pipeline

- [ABRA2](https://github.com/mozack/abra2)
- [BWA](https://github.com/lh3/bwa)
- [fastp](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [ea-utils](https://github.com/ExpressionAnalysis/ea-utils)
- [FastUniq](https://sourceforge.net/projects/fastuniq/)
- [Minimap2](https://github.com/lh3/minimap2)
- [perldup](https://github.com/richardmleggett/scripts/blob/master/remove_pcr_duplicates.pl) Already copied to dependencies folder with permission. Developed by [Richard Leggett](https://github.com/richardmleggett).
- [PRINSEQ](http://prinseq.sourceforge.net/)
- [sambamba](https://github.com/biod/sambamba)
- [samtools](https://github.com/samtools/samtools)
- [trimmomatic](https://github.com/timflutre/trimmomatic)
- [bamtools](https://github.com/pezmaster31/bamtools)







### Wochenende output

Wochenende produces many output files, many of which are superseded by later output files and can be removed.

```
Initial quality checks and read filtering.
- MB_AERO_044_S70_R1.ndp.fastq                  # Fastq after removal of duplicates by Perldup
- MB_AERO_044_S70_R1.ndp.lc.fastq               # Fastq after removal of low-complexity sequences by Prinseq
- MB_AERO_044_S70_R1.ndp.lc_seqs.fq.fastq       # The removed low-complexity output sequences from Prinseq
- MB_AERO_044_S70_R1.ndp.lc.trm.s.bam.unmapped.fastq    # Unmapped reads in FASTQ format. Can be further analysed, eg with alternative programs such as nextflow-blast, kraken, centrifuge etc

BAMs, Mapping Quality (MQ), Duplicate filtering (dup) and mismatch (mm) filtering results
- MB_aero_S2_R1.fastq               # Input file Read1. Note the form R1.fastq is required, R1_001.fastq will not work well.
- MB_aero_S2_R1.fastqprogress.tmp   # Temporary file with pipeline stage progress
- MB_aero_S2_R1.trm.bam             # Initial, unsorted BAM. Can usually be deleted !
- MB_aero_S2_R1.trm.fastq           # Trimmed FASTQ.
- MB_aero_S2_R1.trm.s.bam           # Sorted BAM output file
- MB_aero_S2_R1.ndp.lc.trm.s.bam.unmapped.fastq.gz
- MB_aero_S2_R1.trm.s.mq30.bam                    # BAM where only well mapped reads with Mapping Quality 30 are retained.
- MB_aero_S2_R1.trm.s.mq30.01mm.bam               # Reads with more than 0 or 1 mismatches (ie 2+) have been excluded
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam           # Duplicates excluded
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam.txt       # Important: input for simple runbatch_metagen_awk_filter.sh and complex Wochenende reporting
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam.txt.filt.sort.csv           # Filtered and sorted BAM.txt read output
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam.txt.reporting.sorted.csv    # Output from Wochenende reporting step
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam.txt.reporting.unsorted.csv  # Output from Wochenende reporting step
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.calmd.bam     # MD tags have been calculated. Suitable for viewing SNVs in JBrowse etc
- MB_aero_S2_R2.fastq                 # Read 2 file
- MB_aero_S2_R2.trm.fastq             # Trimmed Read 2 file




BAM Indices
- MB_aero_S2_R1.trm.s.bam.bai       # Index
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.bam.bai       # Index
- MB_aero_S2_R1.trm.s.mq30.01mm.dup.calmd.bam.bai # Index
```


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

### Running Wochenende_plot

Either run the three scripts all together with wochenende_posthoc_filter.sh
```
bash wochenende_posthoc_filter.sh
```

Or run each stage manually:

First generate the data files for wochenende_plot.py

```
# in a directory full of *dup.bam files
bash runbatch_sambamba_depth.sh

Then 
bash runbatch_metagen_window_filter.sh

```

```
MB_AERO_044_S70_R1.ndp.lc.trm.s.mq30.01mm.dup_cov_window.txt              # Coverage per window in each BAM
MB_AERO_044_S70_R1.ndp.lc.trm.s.mq30.01mm.dup_cov_window.txt.filt.csv     # Filtered (regions have at least 1+ reads) coverage per window in each BAM
MB_AERO_044_S70_R1.ndp.lc.trm.s.mq30.01mm.dup_cov_window.txt.filt.sort.csv  # Filtered and sorted (descending) coverage per window
```

### Wochenende_plot output




### Running tests

From the Wochenende directory
```
Using SLURM:
sbatch run_Wochenende_SLURM.sh testdb/reads_R1.fastq
Or without a scheduler:
python3 run_Wochenende.py --metagenome testdb --threads 4 --testWochenende --aligner bwamem --mq30 --remove_mismatching --readType SE --debug --force_restart testdb/reads_R1.fastq
```
You should be able to see in the SLURM outfile or standard out if the tests passed or not. Failed tests may be due to program versions or pipeline configuration issues.





### Contributors

Thanks to:

@B1T0 Main programmer, testing, evaluation, documentation

@Nijerik Wochenende reporting

@konnosif Plots visualisation

@colindaven Concept, programming, updates, integration, maintenance, evaluation, documentation

@Colorstorm Programming, testing, maintenance
