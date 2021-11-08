#!/bin/bash
# Growth rate module. 

# Setup SLURM using data parsed from config.yaml
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)
# Setup job scheduler
# use SLURM job scheduler (yes, no)
if [[ "${USE_CUSTOM_SCHED}" == "yes" ]]; then
    #echo USE_CUSTOM_SCHED set"
    scheduler=$CUSTOM_SCHED_CUSTOM_PARAMS_SINGLECORE
fi
if [[ "${USE_SLURM}" == "yes" ]]; then
    #echo USE_SLURM set"
    scheduler=$SLURM_CUSTOM_PARAMS_SINGLECORE
fi

# check if input bams exist
count_bam=`ls -1 ../*calmd.bam 2>/dev/null | wc -l`

echo "INFO: Starting growth rate analysis module"
if [ "$count_bam" != 0 ]
  then
  for i in ../*calmd.bam
    do
    $scheduler bash run_bed_to_csv.sh "$i"
  done
else
  echo "ERROR: Bam files ../*calmd.bam not found. Cannot run growth rate module and convert to pos.csv"
fi

echo "Growth rate: Remove temp data after samples completed"
# check if files exist
count_bed=`ls -1 *.bed 2>/dev/null | wc -l`
count_subsamp=`ls -1 ../*_subsamples 2>/dev/null | wc -l`

if [ $count_bed != 0 ]
  then
  rm *.bed
fi
if [ $count_subsamp != 0 ]
  then
  rm -rf *_subsamples/
fi
