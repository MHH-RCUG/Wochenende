#!/bin/bash

## Set variables
#cpus=1
## with SLURM srun, default
slurm_srun="srun -c $cpus"

# check if input exists
count_bam=`ls -1 ../*calmd.bam 2>/dev/null | wc -l`

if [ $count_bam != 0 ]
  then
  for i in `ls ../*calmd.bam`
    do
    $slurm_srun -c 2 run_bed_to_csv.sh $i
  done
else
  echo "No bam file found. Can't convert to pos.csv"
fi
