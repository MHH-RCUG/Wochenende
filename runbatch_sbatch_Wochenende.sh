#!/bin/bash

## Colin Davenport, Dec. 2017 - Nov 2021
## Submit all R1.fastq in directory to sbatch run_Wochenende

echo "INFO: Starting Wochenende"

echo "INFO: First checking which analysis mode is selected - only ONE should be selected."
egrep '^\$scheduler python3' run_Wochenende_SLURM.sh
analysisCount="$(egrep '^\$scheduler python3 run_Wochenende.py' run_Wochenende_SLURM.sh | wc -l )"
echo "INFO: Argument count: " ${analysisCount}
if [ $analysisCount != "1" ]; then
        echo "#############################################"
        echo "ERROR: Too many run_Wochenende.py jobs uncommented in run_Wochenende_SLURM.sh"
        echo "ERROR: Stop the script now and check run_Wochenende_SLURM.sh"
        echo "#############################################"
        exit 1
fi

# Remember to check specified a) refseq b) threads c) adapters

echo "INFO: Pipeline check OK, starting job submission "


for i in `ls *R1.fastq`

        do
                echo $i
                sbatch run_Wochenende_SLURM.sh $i

done


