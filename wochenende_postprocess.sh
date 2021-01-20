#!/bin/bash
# Colin Davenport

version="0.1 Jan 2021"
echo "Version: " $version
echo "Remember to run this using the wochenende conda environment"
echo "WORK IN PROGRESS ! . May not completely work for all steps, still useful."


# Operates on Wochenende output files (eg BAM and bam.txt files)
# Postprocess all files 
# Makes sure directories plots and reporting exist
# Run following tools
# - sambamba depth
# - Wochenende reporting
# - Haybaler
# - cleanup directories 
# - Wochenende plot 


# Setup conda and directories
haybaler_dir=/mnt/ngsnfs/tools/dev/haybaler/
wochenende_dir=/mnt/ngsnfs/tools/dev/Wochenende/
# Use existing conda env
. /mnt/ngsnfs/tools/miniconda3/etc/profile.d/conda.sh
conda activate wochenende

# Setup sleep duration
sleeptimer=12
#sleeptimer=120


### Check if required directories exist ###
if [ ! -d "reporting" ] 
then
    echo "INFO: Creating directory reporting" 
    mkdir reporting
fi
if [ ! -d "plots" ] 
then
    echo "INFO: Creating directory plots" 
    mkdir plots
fi

# get current dir containing Wochenende BAM and bam.txt output
bamDir=$(pwd)

echo "INFO: Starting Wochenende_postprocess"
echo "INFO: Current directory" $bamDir
sleep 3

echo "INFO: Started Sambamba depth"
bash runbatch_metagen_awk_filter.sh
wait
bash runbatch_sambamba_depth.sh >/dev/null 2>&1
wait
echo "INFO: Sleeping for " $sleeptimer
sleep $sleeptimer
runbatch_metagen_window_filter.sh >/dev/null 2>&1
wait
echo "INFO: Completed Sambamba depth and filtering"


# Run reporting and haybaler
echo "INFO: Started Wochenende reporting"
cd reporting
cp ../*.bam.txt .
bash runbatch_Wochenende_reporting.sh >/dev/null 2>&1
wait
echo "INFO: Sleeping for " $sleeptimer
sleep $sleeptimer
#echo "INFO: Sleeping for " $sleeptimer " * 10"
#sleep $((sleeptimer*10))


echo "INFO: Completed Wochenende reporting"

echo "INFO: Start Haybaler"
conda activate ont
cp $haybaler_dir/*.sh .
cp $haybaler_dir/*.py .
mv output output_$RANDOM 
bash run_haybaler.sh
wait
echo "INFO: Sleeping for " $sleeptimer
sleep $sleeptimer

echo "INFO: Start csv to xlsx conversion"
bash runbatch_csv_to_xlsx.sh >/dev/null 2>&1
wait
echo "INFO: Sleeping for " $sleeptimer
sleep $sleeptimer
echo "INFO: Completed Haybaler"





# Plots
echo "INFO: Started Wochenende plot"
cd $bamDir
cd plots
cp ../*window.txt ../*window.txt.filt.csv .
bash runbatch_wochenende_plot.sh >/dev/null 2>&1
#wait
echo "INFO: Sleeping for " $sleeptimer
sleep $sleeptimer
cd $bamDir
echo "INFO: Completed Wochenende plot"



echo "INFO: Start cleanup reporting"
cd $bamDir
cd reporting
mkdir reporting_$RANDOM
mv txt csv xlsx reporting_$RANDOM 
mkdir txt csv xlsx
mv *.txt txt
mv *.csv csv
mv *.xlsx xlsx
cd $bamDir
echo "INFO: Completed cleanup reporting"

echo "INFO: ########### Completed Wochenende_postprocess #############"

