#!/bin/bash

## Set variables
cpus=1
## with SLURM srun, default
slurm_srun="srun -c $cpus"

# check if input exists
count_bam=`ls -1 ../*calmd.bam 2>/dev/null | wc -l`

echo "INFO: Starting growth rate analysis module"
if [ "$count_bam" != 0 ]
  then
  for i in ../*calmd.bam
    do
    $slurm_srun bash run_bed_to_csv.sh "$i"
  done
else
  echo "ERROR: Bam files ../*calmd.bam not found. Cannot run growth rate module and convert to pos.csv"
fi
