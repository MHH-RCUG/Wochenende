#!/bin/bash

## Colin, Dec. 2017 - Dec 2020
## Submit all FASTQ in directory to sbatch run_Wochenende

echo "INFO: Starting Wochenende"

echo "INFO: First checking which analysis mode is selected - only ONE should be selected."
egrep ^python3 run_Wochenende_SLURM.sh
#$cmd
echo "INFO: Argument count: "

analysisCount="$(egrep ^'python3 run_Wochenende.py' run_Wochenende_SLURM.sh | wc -l )"
echo "${analysisCount}"
if [ $analysisCount != "1" ]; then
        echo "#############################################"
        echo "ERROR: Too many run_Wochenende.py jobs uncommented in run_Wochenende_SLURM.sh"
        echo "ERROR: Stop the script now and check run_Wochenende_SLURM.sh"
        echo "#############################################"
        exit
fi

# Remember to check specified a) refseq b) threads c) adapters

echo "INFO: Pipeline check OK, starting job submission "


for i in `ls *R1.fastq`

        do
                echo $i
                sbatch -c 16 run_Wochenende_SLURM.sh $i

done


