mkdir stats
mkdir fastqc


rm *.tmp
rm *.trm.fastq
rm *.trm.bam
rm *.trm.s.bam
rm *.trm.s.bam.bai
mv *.bam.stats stats
mv *.bam.txt stats
mv *_out fastqc
srun -c 56 pigz -p 56 *.fastq &


