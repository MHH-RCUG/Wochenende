#!/bin/bash
## Detects all bam.txt files in the current folder and runs the reporting SLURM script


for bamtxt in `ls *.bam.txt`
        do
        # run local
        bash run_Wochenende_reporting_SLURM.sh $bamtxt &
done
