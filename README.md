# Wochenende - A whole Genome/Metagenome Sequencing Alignment Pipeline

Wochenende runs alignment of short reads (eg Illumina) or long reads (eg Oxford Nanopore) against a reference sequence. It is simple (python script), portable and has many optional steps. 

Features include (see programs listed below at the bottom of this page)
- QC (Fastqc)
- pre alignment duplicate removal (perldup)
- pre alignment poor sequence removal (Prinseq - used for single ended reads only)
- trimming (Trimmomatic or fastp or trim galore or ea-utils)
- alignment (bwa mem, minimap2 or ngmlr)
- SAM-> BAM conversion (samtools and sambamba)
- AlignerBoost Mapping Quality recalculation 
- Report % aligned reads (samtools)
- Output unmapped reads as fastq (samtools)  (from v1.4)
- Post-alignment duplicate removal (Samtools from v1.7.8, Sambamba)
- Removal reads with x mismatches (bamtools), adjustable from v1.7.3
- Realignment (Abra2)
- MD tag marking (Samtools)
- Normalization (to Reads per Human cell, Reads Per Million reads per Million bases etc, see Reporting below for details)
- Visualization (chromosome coverage, intended for bacteria in metagenomics projects) (from v1.4)

Project Haybaler https://github.com/MHH-RCUG/haybaler allows postprocessing of Wochenende results:
- collation/integration of multiple reports (reporting csv or bam.txt files) using Python Pandas
- prepare results for heatmaps
- create heatmaps using multiple different R libraries

### Did you know ? 
Wochenende means weekend in German. The original developer, Tobias, called the pipeline Wochenende, because you can start it running and go off to enjoy your weekend early.

![Alt text](dependencies/wochenende_2.png?raw=true "Wochenende schematic")


## Platforms

Wochenende has currently (July 2020) only been tested on Ubuntu Linux 20.04 and 16.04 64bit. We advise against any attempts on MacOS or Windows. An appropriate conda environment, BASH and Python3.6+ is critical to get things working.


## Usage 

You can just run the pipeline as a normal Python3 script. However, we also offer a template for the job scheduler SLURM below. This template can also be used with Bash to run the commands at the bottom of the SLURM pipeline while ignoring any SLURM specific parameters.

### SLURM usage

1. Copy all the run_Wochenende* files and prerequisite subfolders to your directory with your FASTQ files
`cp /path/to/wochenende/get_wochenende.sh .` `bash get_wochenende.sh`
2. Adjust settings in the script
`nano run_Wochenende_SLURM.sh`
3. Run the pipeline using SLURM (the "_R1" is important)
`sbatch run_Wochenende_SLURM.sh sample_R1.fastq`
4. Optional reporting step to normalize the extracted read counts  (see reporting below)
`sbatch run_Wochenende_reporting_SLURM.sh`

### Tutorial

Once you've got the tools installed and tested, you can look at or run the commands in the tutorial in the subdirectory `tutorial`. https://github.com/MHH-RCUG/Wochenende/blob/master/tutorial/tutorial.txt

### General usage

```
Warning, this usage is just an example and might be slightly out of date. 

Run this with: 
python3 run_Wochenende.py

Wochenende - Whole Genome/Metagenome Sequencing Alignment Pipeline
Wochenende was created by Dr. Colin Davenport, Tobias Scheithauer and Fabian Friedrich with help from many further contributors https://github.com/MHH-RCUG/Wochenende/graphs/contributors
version: 1.9.1 - Mar 2021

usage: run_Wochenende.py [-h] [--aligner {bwamem,minimap2,ngmlr}]
                         [--readType {PE,SE}]
                         [--metagenome {2021_02_meta_fungi_human_masked,2021_02_meta_fungi_human_unmasked,2020_09_massiveref_human,2020_05_meta_human,2020_03_meta_human,2019_01_meta,2019_10_meta_human,2019_10_meta_human_univec,2019_01_meta_mouse,2019_01_meta_mouse_ASF_OMM,2019_01_meta_mouse_ASF,2019_01_meta_mouse_OMM,hg19,GRCh37,GRCh38-45GB,GRCh38-noalt,GRCh38-mito,mm10,rn6,rat_1AR1_ont,zf10,ss11,PA14,ecoli,nci_viruses,ezv_viruses,testdb,strept_halo,k_variicola,k_oxytoca,clost_bot,clost_bot_e,clost_diff,clost_perf,citro_freundii}]
                         [--threads THREADS] [--fastp] [--nextera]
                         [--trim_galore] [--debug] [--longread]
                         [--no_duplicate_removal] [--no_prinseq] [--no_fastqc]
                         [--no_abra] [--mq20] [--mq30]
                         [--remove_mismatching REMOVE_MISMATCHING]
                         [--force_restart] [--testWochenende]
                         fastq

positional arguments:
  fastq                 _R1.fastq Input read1 fastq file

optional arguments:
  -h, --help            show this help message and exit
  --aligner {bwamem,minimap2,ngmlr}
                        Aligner to use, either bwamem, ngmlr or minimap2.
                        Usage of minimap2 and ngmlr currently optimized for
                        nanopore data only.
  --readType {PE,SE}    Single end or paired end data
  --metagenome {2021_02_meta_fungi_human_masked,2021_02_meta_fungi_human_unmasked,2020_09_massiveref_human,2020_05_meta_human,2020_03_meta_human,2019_01_meta,2019_10_meta_human,2019_10_meta_human_univec,2019_01_meta_mouse,2019_01_meta_mouse_ASF_OMM,2019_01_meta_mouse_ASF,2019_01_meta_mouse_OMM,hg19,GRCh37,GRCh38-45GB,GRCh38-noalt,GRCh38-mito,mm10,rn6,rat_1AR1_ont,zf10,ss11,PA14,ecoli,nci_viruses,ezv_viruses,testdb,strept_halo,k_variicola,k_oxytoca,clost_bot,clost_bot_e,clost_diff,clost_perf,citro_freundii}
                        Meta/genome reference to use
  --threads THREADS     Number of threads to use
  --fastp               Use fastp trimmer instead of fastqc and trimmomatic
  --nextera             Attempt to remove Illumina Nextera adapters and
                        transposase sequence (default is Illumina Ultra II
                        adapters, but Illumina Nextera more common in future)
  --trim_galore         Use trim_galore read trimmer. Effective for Nextera
                        adapters and transposase sequence
  --debug               Report all files
  --longread            Only do steps relevant for long PacBio/ONT reads eg.
                        no dup removal, no trimming, just alignment and bam
                        conversion
  --no_duplicate_removal
                        Skips steps for duplicate removal. Recommended for
                        amplicon sequencing.
  --no_prinseq          Skips prinseq step (low_complexity sequence removal)
  --no_fastqc           Skips FastQC quality control step.
  --no_abra             Skips steps for Abra realignment. Recommended for
                        metagenome and amplicon analysis.
  --mq20                Remove reads with mapping quality less than 20.
                        Recommended for metagenome and amplicon analysis. Less
                        stringent than MQ30.
  --mq30                Remove reads with mapping quality less than 30.
                        Recommended for metagenome and amplicon analysis.
  --remove_mismatching REMOVE_MISMATCHING
                        Remove reads with less than x mismatches (via the NM
                        bam tag). Default 3. Argument required.
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
2. Create a conda environment for the pipeline. You should have first installed miniconda 64-bit Linux.
```
cd wochenende
conda env create -f env.wochenende.yml
# or if you have problems remove the previous and try this:
conda env create -f env.wochenende.minimal.yml
```
3. Install all the other tools.
   - [ABRA2](https://github.com/mozack/abra2)
4. Important! Edit the configuration section of `run_Wochenende.py` to set the paths to the tools, tmp directory and reference sequences.
5. (Wochenende_plot only). Install the python dependencies for visualization by pip. The works on Ubuntu 1604: `pip3 install --user numpy==1.12.1 numpy-base==1.15.0 pandas==0.23.4 matplotlib==2.2.2`
6. Activate the conda environment before running the pipeline.
`conda activate wochenende`
7. Optional: run the tests, see below.

### Update conda enviroment
If there is already a conda enviroment named wochenende:
```
conda env update -f env.wochenende.yml
```

## List of Tools used or optional in the pipeline

- [Alignerboost](https://github.com/Grice-Lab/AlignerBoost). GPL3, in dependencies folder.
- [ABRA2](https://github.com/mozack/abra2)
- [bamtools](https://github.com/pezmaster31/bamtools)
- [BWA](https://github.com/lh3/bwa)
- [ea-utils](https://github.com/ExpressionAnalysis/ea-utils)
- [fastp](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [FastUniq](https://sourceforge.net/projects/fastuniq/)
- [Minimap2](https://github.com/lh3/minimap2)
- [NGMLR](https://github.com/philres/ngmlr)
- [perldup](https://github.com/richardmleggett/scripts/blob/master/remove_pcr_duplicates.pl) Already copied to dependencies folder with permission. Developed by [Richard Leggett](https://github.com/richardmleggett).
- [PRINSEQ](http://prinseq.sourceforge.net/)
- [sambamba](https://github.com/biod/sambamba)
- [samtools](https://github.com/samtools/samtools)
- [trim_galore](https://github.com/FelixKrueger/TrimGalore)
- [trimmomatic](https://github.com/timflutre/trimmomatic)

Postprocessing
- [Haybaler](https://github.com/MHH-RCUG/haybaler)








### Wochenende output

Wochenende produces many output files, many of which are superseded by later output files and can be removed.

```
Initial quality checks and read filtering.
- MB_AERO_044_S70_R1.ndp.fastq                  # Fastq after removal of duplicates by Perldup
- MB_AERO_044_S70_R1.ndp.lc.fastq               # Fastq after removal of low-complexity sequences by Prinseq
- MB_AERO_044_S70_R1.ndp.lc_seqs.fq.fastq       # The removed low-complexity output sequences from Prinseq
- MB_AERO_044_S70_R1.ndp.lc.trm.s.bam.unmapped.fastq    # Unmapped reads in FASTQ format. Can be further analysed, eg with alternative programs such as nextflow-blast, kraken, centrifuge etc

BAMs, Mapping Quality (MQ), Duplicate filtering (dup) and mismatch (mm) filtering results
- MB_aero_S2_R1.fastq               # Input file Read1. Note the form R1.fastq is required, R1_001.fastq will not work.
- MB_aero_S2_R1.fastqprogress.tmp   # Temporary file with pipeline stage progress
- MB_aero_S2_R1.trm.bam             # Initial, unsorted BAM. Can usually be deleted !
- MB_aero_S2_R1.trm.fastq           # Trimmed FASTQ.
- MB_aero_S2_R1.trm.s.bam           # Sorted BAM output file
- MB_aero_S2_R1.ndp.lc.trm.s.bam.unmapped.fastq   # unmapped FASTQ reads
- MB_aero_S2_R1.trm.s.mq30.bam                    # BAM where only well mapped reads with Mapping Quality 30 are retained
- MB_aero_S2_R1.trm.s.mq30.mm.bam               # Reads with 3 or more mismatches (default, can be changed) have been excluded
- MB_aero_S2_R1.trm.s.mq30.mm.dup.bam           # Duplicate reads were excluded by Picard
- MB_aero_S2_R1.trm.s.mq30.mm.dup.calmd.bam     # MD tags have been calculated to enable SNV visualization in JBrowse etc
- MB_aero_S2_R2.fastq                 # Read 2 input file
- MB_aero_S2_R2.trm.fastq             # Trimmed Read 2 file

# Wochenende reporting input and output
- MB_aero_S2_R1.trm.s.mq30.mm.dup.bam.txt       # Important: input for simple runbatch_metagen_awk_filter.sh and complex Wochenende reporting
- MB_aero_S2_R1.trm.s.mq30.mm.dup.bam.txt.filt.sort.csv           # Filtered and sorted BAM.txt read output
- MB_aero_S2_R1.trm.s.mq30.mm.dup.bam.txt.reporting.sorted.csv    # Output from Wochenende reporting step
- MB_aero_S2_R1.trm.s.mq30.mm.dup.bam.txt.reporting.unsorted.csv  # Output from Wochenende reporting step


# Wochenende_plot.py input (.filt.csv) and output (png images)
- MB_AERO_044_S70_R1.ndp.lc.trm.s.mq30.mm.dup_cov_window.txt              # Coverage per window in each BAM
- MB_AERO_044_S70_R1.ndp.lc.trm.s.mq30.mm.dup_cov_window.txt.filt.csv     # Filtered (regions have at least 1+ reads) coverage per window in each BAM
- MB_AERO_044_S70_R1.ndp.lc.trm.s.mq30.mm.dup_cov_window.txt.filt.sort.csv  # Filtered and sorted (descending) coverage per window

# Wochenende_plot.py output (png images)
- wochenende_png_files/
- wochenende_png_files/sample1.dup_cov_window.txt.filt.csv/
- wochenende_png_files/sample1.dup_cov_window.txt.filt.csv/perhaps_present/  
- wochenende_png_files/sample1.dup_cov_window.txt.filt.csv/probably_present/


```

## Run the postprocessing automatically (recommended!)

After a successful Wochenende run, make sure you check that all bams have been created and are sized as expected eg `ls -lh *.bam`

Now start the postprocessing script `bash wochenende_postprocess.sh`
- run sambamba depth to get read coverage of all configured BAM files in the current directory
- run the Wochenende plot to create coverage diagrams
- run Wochenende reporting to count and normalize all read data
- run the Haybaler report integration tool (provided it is installed and configured)
- clean up files

This script requires haybaler to be installed, and will otherwise fail.


### Running the metagenomic reporting scripts manually 

This tool (which requires python v3.6+) reports length, GC content of the sequence, read counts attributed to the species and various normalized read count parameters. 
Normalizations are for:

a) reads normalized to the idealized length of a bacterial chromosome (normalization to 1 million base pairs)

b) total reads in the sequencing library (normalization to 1 million reads)

c) the above two normalizations combined (RPMM)

d) reads per human cell (only works for metagenomes from human hosts)

See the subfolder `reporting` in the repository.

```
conda activate 
conda activate wochenende
python3 basic_reporting.py --input_file tmp_R1.ndp.lc.trm.s.mq30.01mm.dup.bam.txt --reference /lager2/rcug/seqres/metagenref/2016_06_PPKC_metagenome_test_1p_spec_change_cln.fa --sequencer illumina --output_name test
```

#### Usage
```
Usage: basic_reporting.py [OPTIONS]

  This script can be used to report the results of the Wochenende pipeline.
  The .bam.txt file as input is recommended. The .bam file will take longer
  and  generate more information.

  The column reads_per_human_cell is only for metagenomes from human hosts.

  Reports for solid sequencing data are not supported, a special
  normalisation model has to be implemented first.

Options:
  -i, --input_file TEXT   File in .bam.txt or .bam format from the Wochenende
                          pipeline output
  -r, --reference TEXT    File in .fasta format has to be the reference used
                          by the Wochenende pipeline
  -s, --sequencer TEXT    Sequencer technology used only solid and illumina
                          are available, only illumina is supported, default:
                          illumina
  -o, --output_name TEXT  Name for the output file(sample name), default
                          report
  --help                  Show this message and exit.
```
### Running Wochenende_plot

#### Preparing the data from BAM files

First generate the data files for wochenende_plot.py

```
# in a directory full of *dup.bam files
bash runbatch_sambamba_depth.sh

Then 
bash runbatch_metagen_window_filter.sh
```
The result should be a series of output files containing the keyword *window*

#### Run the actual plotting 

Finally, run the actual wochenende_plot.py script or the helper bash script.

```
bash runbatch_wochenende_plot.sh
```

wochenende_plot.py usage:
```
python3 wochenende_plot.py
usage: wochenende_plot.py [-h] [--minMeanCov MINMEANCOV]
                          [--createAllPngs CREATEALLPNGS] [--sclim SCLIM]
                          [--minWindows MINWINDOWS]
                          filename1
wochenende_plot.py: error: the following arguments are required: filename1
```


### Wochenende_plot output

Wochenende_plot creates one subdirectory per input file. These contain png images of taxa which are probabably (high score, largely based on consistent evenness of coverage and high mean coverage) or perhaps present (need manual review). Confident attributions to taxa depends strongly on the number of reads assigned to bacterial taxa (low in airway metagenomes, higher in for example stool samples).  


### Known bugs

RPMM bug: fixed in v1.7.8. In October 2020 a bug in the Wochenende_reporting script was found which calculated the RPMM column incorrectly. Please recalculate your reporting statistics if you use this feature. Thanks to @sannareddyk and colleagues at Leibniz University Hannover.

### Running software tests

These tests test the software and installation on your server with known read set and database. If all is working correctly, then an output message should be created that these tests passed.

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

@B1T0 Original programmer, testing, evaluation, documentation

@colindaven Concept, programming, updates, integration, maintenance, evaluation, documentation

@Colorstorm Programming, testing, maintenance

@konnosif Plots visualisation

@Nijerik Wochenende reporting

@sannareddyk Bug testing, updates, evaluation

@poer-sophia Code review, testing, maintenance, programming (haybaler and more)
