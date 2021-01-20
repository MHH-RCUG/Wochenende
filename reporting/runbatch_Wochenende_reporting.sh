#!/bin/bash
## Detects all bam.txt files in the current folder and runs the reporting SLURM script

echo "INFO: Starting batch reporting"

for bamtxt in `ls *.bam.txt`
        do
        # run local
        bash run_Wochenende_reporting_SLURM.sh $bamtxt >/dev/null 2>&1 &
done

echo "INFO: Waiting for SLURM scripts to complete"
wait

echo "INFO: Completed batch reporting"

