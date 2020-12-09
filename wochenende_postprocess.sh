#!/bin/bash
# Colin Davenport, December 2020

# Operates on Wochenende output files (eg BAM, bam.txt files)
# Postprocess all files 
# Run following tools
# - Wochenende reporting
# - Haybaler
# - cleanup directories 
# - Wochenende plot 


# Setup directories
haybaler_dir=/mnt/ngsnfs/tools/dev/haybaler/
wochenende_dir=/mnt/ngsnfs/tools/dev/Wochenende/

# Use existing conda env
conda activate wochenende

bamDir=$(PWD)

echo "INFO: Starting Wochenende_postprocess"
echo "INFO: Current directory" $bamDir

echo "INFO: Started Sambamba depth"
bash runbatch_metagen_awk_filter.sh
wait
bash runbatch_sambamba_depth.sh
wait
echo "INFO: Sleeping 120s"
sleep 120
echo "INFO: Completed Sambamba depth"


# Run reporting and haybaler
echo "INFO: Started Wochenende reporting"
cd reporting
cp ../*.bam.txt .
srun bash runbatch_Wochenende_reporting.sh &
wait
echo "INFO: Sleeping 1000s"
sleep 1000
echo "INFO: Completed Wochenende reporting"

echo "INFO: Start Haybaler"
cp $haybaler_dir/*.sh .
cp $haybaler_dir/*.py .
mv output output_$RANDOM 
bash run_haybaler.sh
wait
echo "INFO: Sleeping 120s"
sleep 120
bash runbatch_csv_to_xlsx.sh
wait
echo "INFO: Sleeping 120s"
sleep 120
echo "INFO: Completed Haybaler"





# Plots
echo "INFO: Started Wochenende plot"
cd $bamDir
cd plots
cp *window.txt .
bash runbatch_wochenende_plot.sh
#wait
echo "INFO: Sleeping 120s"
sleep 120
cd $bamDir
echo "INFO: Completed Wochenende plot"



echo "INFO: Start cleanup reporting"
cd reporting
rm -rf txt csv xlsx
mkdir txt csv xlsx
mv *.txt txt
mv *.csv csv
mv *.xlsx xlsx
cd $bamDir
echo "INFO: Completed cleanup reporting"

echo "INFO: Completed Wochenende_postprocess"

