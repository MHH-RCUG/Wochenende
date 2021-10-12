![Alt text](raspir_image.jpg?raw=true)
# Background
In shotgun metagenomic sequencing experiments, the total DNA is extracted from complex samples. The DNA of rare species is scarce, and poor coverage of these genomes is often observed. Scientists commonly define thresholds to exclude 1-10% of the least abundant species in a given sample from further analyses. On the one hand, this filtering step allows for robust investigations of core communities. On the other hand, valuable information on the community structure is lost. The rare biosphere harbours more species than the core microbial community and hence provides the environment of interest with high functional flexibility.

However, if only a few short DNA reads are detected that are unique to species A, there are at least three explanations: <br>
a) Sample contamination; <br>
b) Rare species A was present in the environment of interest, so it is a true positive species. In this case, the reads are expected to spread across the entire reference genome in a fairly uniform manner due to random DNA sequencing; or <br>
c) Rare species A was absent (false positive) but rare species B was present, which acquired genes of species A during past events. In this case, the reads are expected to cluster at specific sides of the reference genome of species A. <br>

The raspir tool calculates a position-domain signal based on the distances of reads aligning to a circular reference genome and converts the information into a frequency signal through Discrete Fourier Transforms (DFT). In addition, a reference frequency signal is constructed with the same number of reads, but with an ideal uniform distribution of reads across the reference genome. Both frequency signals are compared using Pearson's correlation measures. 

# Implementation
Raspir is implemented in Python 3.7. Using this tool, a real-world dataset (5 GB) containing information on hundreds of species can be processed on a single node server. The input data must be structured in the following manner: genome length of the corresponding reference genome, organism, read position, read depth. 

See the following section for further information on how to convert .FASTQ files into the .CSV input files to sucessfully execute raspir.

# Get started
## Set up the environment
### Install conda packages 

```bash
conda create --name raspir_env
conda activate raspir_env

# Install trimmomatic [1] 
conda install -c bioconda trimmomatic
# Install samtools [2] 
conda install -c bioconda samtools

# Install your alignment tool of choice
# Burrows-Wheeler aligner [3]
conda install -c bioconda bwa
# Bowtie2 [4]
conda install -c bioconda/label/cf201901 bowtie2

# Install python packages for raspir
conda install pandas
conda install -c conda-forge statsmodels
conda install -c conda-forge matplotlib
```


### Create working directory

```bash
mkdir raspir/
cd raspir/

YOURPATH=$PWD
echo "$YOURPATH"

mkdir reference_database/
mkdir run_raspir/
```

### Set-up the reference database <br>
The reference database has to be downloaded and indexed only once. 
The database contains only complete bacterial genomes with one strain per species.
Note: You may also use a customised reference database. It is however strongly recommended to avoid draft or low-quality reference genomes and use complete sequences of circular microorganisms only.

```bash
# Load database into your working directory
cd reference_database/
# Get database from
https://sync.academiccloud.de/index.php/s/DcKLpKtoKMHKvcn/download
# Unzip the reference database 
gunzip complete_bacterialRefSeqs_201910_3.fasta

# Generate an index of the reference fasta depending on the alignment tool of your choice
samtools faidx complete_bacterialRefSeqs_201910_3.fasta
bwa index complete_bacterialRefSeqs_201910_3.fasta
bowtie2-build complete_bacterialRefSeqs_201910_3.fasta complete_bacterialRefSeqs_201910_3
cd ..
```

### Sorting, indexing & final clean-up

This process can be performed on SAMs or BAMs with the included script run_SLURM_file_prep.sh

```bash


# Activate env on cluster node
conda activate raspir_env >> /dev/null

cpus=8
human_chr_pattern="1_1_1"

# Only run SAM section if SAM files exist in the current directory
count=`ls -1 *.sam 2>/dev/null | wc -l`
if [ $count != 0 ]
    then
    # run SAM conversion

	# For pipelines which start from SAM. Many pipelines will start from BAM.
	for items in *.sam
		do
			echo $items
			fname=$(echo ${items} | sed 's/.sam//')
			echo $fname


			# Remove reads with low mapping quality
			samtools view -hM -q 20 $items > ${items%.sam}.mq20.sam

			# Convert file from SAM to BAM format
			samtools view -h -b -S ${items%.sam}.mq20.sam  > ${fname}.bam

			# Discard unmapped sequences 
			samtools view -b -F 4 $items > ${fname}_1.bam

			# Sort bam file
			samtools sort @ $cpus -${fname}_1.bam -o ${fname}.sorted.bam
		done
fi

# Start from sorted mq20 BAM with no unmapped sequences  
for items in *.bam
	do
		echo $items
		fname=$(echo ${items} | sed 's/.bam//')
		echo $fname

   		# Obtain coverage information
   		samtools depth ${fname}.bam | grep -v $human_chr_pattern  > ${fname}.raspir1.csv

   		# Add genome size, pipe in a BAM header only
   		samtools view -H ${items} | sed 's/LN://g' > ${fname}.genomeSize_1.csv
   		sed -i 's/SN://g' ${fname}.genomeSize_1.csv
   		cut -f2- ${fname}.genomeSize_1.csv > ${fname}.genomeSize.csv

   		# Add column with genome size
		awk -v FS="\t" -v OFS="\t" 'FNR==NR{a[$1]=$2;next;} {if(a[$1]) {print a[$1], $0} else {print "NA",$0}}' \
		${fname}.genomeSize.csv ${fname}.raspir1.csv > ${fname}.raspir.csv

		# Convert into a comma-separated file
		sed -i 's/\t/,/g' ${fname}.raspir.csv
		# Add header
		sed -i '1iGenomeLength,Organism,Position,Depth\' ${fname}.raspir.csv

		# Remove intermediate files
		rm ${fname}.raspir1.csv ${fname}.genomeSize_1.csv ${fname}.genomeSize.csv
	done

```

# Run raspir
Download the python script into the run_raspir/ folder.

```python
python raspir.py input.csv output_prefix 
```

# Output
A table is generated (.CSV format). The assignment output has 6 columns.

| Species | r_value  | p_value  | stError | euclidean | distribution |
| :---:   | :-: | :-: | :-: | :-: | :-: | 
| Pseudomonas aeruginosa | 0.99 | 0.0 | 0.00019 | 0.01 | uniform |
| Streptococcus salivarius | 0.97 | 0.0 | 0.00016 | 0.002 | uniform |
| Rothia mucilaginosa | 0.99 | 0.0 | 0.000002 | 0.0001 | uniform | 


# Contributors
@mmpust author

@colindaven updates

@nick-youngblut updates 

# References
[1] Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170. <br>
[2] Li H., Handsaker B., Wysoker A. et al. (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. <br>
[3] Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. <br>
[4] Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359. <br>

# Cite the tool
Pust, MM., Tümmler, B. Identification of core and rare species in metagenome samples based on shotgun metagenomic sequencing, Fourier transforms and spectral comparisons. ISME COMMUN. 1, 2 (2021). https://doi.org/10.1038/s43705-021-00010-6
