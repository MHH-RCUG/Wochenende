#!/bin/bash
# Colin Davenport, Lisa Hollstein 2022_01
set -u pipefail

echo "Removing non essential files such as unsorted BAMs and moving logs to subfolders"

base_dir=$(pwd)
echo $base_dir

mkdir -p stats
mkdir -p fastqc
mkdir -p logs

rm *.tmp
rm *.trm.fastq
rm *.trm.bam
#rm *.trm.s.bam
#rm *.trm.s.bam.bai
rm *.ns.bam
rm *.fix.bam
mv *.bam.stats stats
mv *.bam.txt stats
mv *_out fastqc
mv *.out logs

echo "INFO: removing png files from growth_rate"
# a find exec rm is dangerous though, so go more targeted
cd growth_rate
#find . -n
cd $base_dir


echo "INFO: gzipping fastq with pigz"
#assumes SLURM present, else remove all before pigz
srun -c 56 pigz -p 56 *.fastq &


