#!/bin/bash
# Colin Davenport, Lisa Hollstein 2022_01

echo "Removing non essential files such as unsorted BAMs and moving logs to subfolders"


mkdir stats
mkdir fastqc
mkdir logs

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
#assumes SLURM present, else remove all before pigz
srun -c 56 pigz -p 56 *.fastq &


