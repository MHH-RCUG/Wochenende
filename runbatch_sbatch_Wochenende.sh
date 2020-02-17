#!/bin/bash

## Colin, Dec. 2017
## Submit all FASTQ in directory to sbatch run_Wochenende

# Remember to check specified a) refseq b) threads c) adapters

for i in `ls *R1.fastq`

        do
                echo $i
                sbatch -c 16 run_Wochenende_SLURM.sh $i

done


