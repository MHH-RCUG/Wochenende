#!/bin/bash

## Colin, Dec. 2017
## Submit all FASTQ in directory to sbatch run_Wochenende

# Remember to check specified a) refseq b) threads c) adapters

for i in `ls *_R1.fastq`

        do
                echo $i
                sbatch run_Wochenende_slurm.sh $i

done


