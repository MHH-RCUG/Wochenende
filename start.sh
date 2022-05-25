#!/bin/bash
# bash start.sh input.fastq

#git pull && nextflow run longread_alignment.nf
ref="/mnt/omics/data/projects/26890/uploaded/214815-ZR_S706_v1_chromosomes.fasta"
fastq=$1

git pull -q 

# check args
if [ -z $fastq ]
        then
        echo ""
        echo "## Usage: Input fastq required. bash start.sh in.fastq ## "
        exit
fi


# run specifying fasta reference as arg
nextflow run longread_alignment.nf -with-timeline -with-report -with-dag flowchart.png --fasta $ref --fastq $fastq --test yes -resume
