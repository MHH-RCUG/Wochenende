#!/bin/bash
# Automated postprocessing of results from the Wochenende pipeline
# Authors: Colin Davenport, Sophia Poertner

version="0.22, May 2021"

#Changelog
#0.22 - add heat trees
#0.21 - attempt recovery for second runs to copy data from csv or txt subdirs into haybaler dir
#0.20 - add haybaler heat tree support
#0.19 - update haybaler copying and add double square brackets for bash ifs
#0.1xx - TODO Use environment variables for haybaler and wochenende installations
#0.18 - make wochenende_plot optional with --no-plots
#0.17 - check directories and improve haybaler integration
#0.16 - use Haybaler update runbatch_heatmaps.sh
#0.15 - check for files to cleanup before moving
#0.14 - add haybaler env and use this
#0.13 - copy filled directory plots and reporting into current directory, if missing
#0.12 - make prerequisite docs clearer
#0.11 - add variable for random 
#0.10 - initial commits


echo "INFO: Postprocess Wochenende BAM and bam.txt files for plotting, reporting and haybaler integration" 
echo "INFO: Version: " $version
echo "INFO: Usage: bash wochenende_postprocess.sh args"
echo "INFO: Remember to run this using the haybaler conda environment if available - we attempt to load this in the script"
echo "INFO: Usage: Make sure the directories plots/ and reporting/ exist and are filled"
echo "INFO:  eg. run: bash get_wochenende.sh to get the relevant files"
echo "INFO:  ####### "
echo "INFO:  Runs following stages"
echo "INFO:  - sambamba depth"
echo "INFO:  - Wochenende plot (disable with --no-plots argument)"
echo "INFO:  - Wochenende reporting"
echo "INFO:  - Haybaler and heatmaps in R (Haybaler, and R required)"
echo "INFO:  - Haybaler taxonomy and heat-trees in R (Haybaler, pytaxonkit, and R required)"
echo "INFO:  - cleanup directories "

if [[ $1 == "--no-plots" ]]
then
    echo "INFO: Found --no-plots argument: Plot mode disabled"
fi

# Setup conda and directories
haybaler_dir=/mnt/ngsnfs/tools/dev/haybaler/
wochenende_dir=/mnt/ngsnfs/tools/dev/Wochenende/
# Use existing conda env
. /mnt/ngsnfs/tools/miniconda3/etc/profile.d/conda.sh
conda activate wochenende

# Setup sleep duration. Might be useful to set higher for some big projects, where the wait command may fail for some SLURM jobs.
sleeptimer=12
#sleeptimer=120

# Cleanup previous results to a directory with a random name which includes a number, calculated here.
rand_number=$RANDOM

### Check if required directories/files exist, copy if missing ###
if [[ ! -d "reporting" ]] 
then
    echo "INFO: Copying directory reporting, as it was missing!" 
    cp -R $wochenende_dir/reporting .
fi
if [[ ! -d "plots" ]] 
then
    echo "INFO: Copying directory plots, as it was missing!" 
    cp -R $wochenende_dir/plots .
fi
if [[ ! -f "reporting/ref.tmp" ]] 
then
    echo "INFO: Missing file reporting/ref.tmp , attempting to copy ./ref.tmp to reporting/ref.tmp" 
    cp ref.tmp reporting/ref.tmp
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
bash runbatch_metagen_window_filter.sh >/dev/null 2>&1
wait
echo "INFO: Completed Sambamba depth and filtering"


# Plots
if [[ $1 == "--no-plots" ]] 
then
    echo "INFO: Found --no-plots argument: Plot mode disabled"
else

    echo "INFO: Started Wochenende plot"
    cd $bamDir
    cd plots
    cp ../*_window.txt . 
    cp ../*_window.txt.filt.csv .

    bash runbatch_wochenende_plot.sh >/dev/null 2>&1
    #wait
    echo "INFO: Sleeping for " $sleeptimer
    sleep $sleeptimer
    cd $bamDir
    echo "INFO: Completed Wochenende plot"
fi

# Run reporting 
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


# Run haybaler
echo "INFO: Start Haybaler"
if [[ ! -d "haybaler" ]]
    then
    mkdir haybaler
fi
count_mq30=`ls -1 *mq30.bam*us*.csv 2>/dev/null | wc -l`
count_dup=`ls -1 *dup.bam*us*.csv 2>/dev/null | wc -l`
count=`ls -1 *.bam*us*.csv 2>/dev/null | wc -l`
if [[ $count_mq30 != 0 ]]
    then
    cp *mq30.bam*us*.csv haybaler
elif [[ $count_dup != 0 ]]
    then
    cp *dup.bam*us*.csv haybaler
elif [[ $count != 0 ]]
    then
    cp *.bam*us*.csv haybaler
else
    echo "WARNING: No bam*us*.csv found to process for haybaler"
    echo "INFO: Attempting to find and copy bam*us*.csv to haybaler"
    cp txt/*.bam*us*.csv haybaler
    cp csv/*.bam*us*.csv haybaler
fi
cd haybaler
conda activate haybaler
cp $haybaler_dir/*.sh .
cp $haybaler_dir/*.py .
cp $haybaler_dir/*.R .
bash run_haybaler.sh $haybaler_dir >/dev/null 2>&1
wait
cp $haybaler_dir/runbatch_heatmaps.sh haybaler_output/ && cp $haybaler_dir/*.R haybaler_output/
cp $haybaler_dir/*tax* haybaler_output/
cp $haybaler_dir/*tree* haybaler_output/

echo "INFO: Attempting to filter results and create heatmaps. Requires R installation." 
cd haybaler_output
bash runbatch_heatmaps.sh  
echo "INFO: Attempting to add taxonomy. Requires R installation and pytaxonkit." 
bash run_haybaler_tax.sh
echo "INFO: Attempting create heat-trees. Requires R installation and packages: packages = c("metacoder", "taxa", "dplyr", "tibble", "ggplot2")." 
bash run_heattrees.sh
cd ..
cd ..
echo "INFO: Sleeping for " $sleeptimer
sleep $sleeptimer

echo "INFO: Start csv to xlsx conversion"
bash runbatch_csv_to_xlsx.sh >/dev/null 2>&1
wait
echo "INFO: Sleeping for " $sleeptimer
sleep $sleeptimer
echo "INFO: Completed Haybaler"





echo "INFO: Start cleanup reporting"
cd $bamDir
cd reporting
# create backup, move folders from previous reporting run to a directory
mkdir reporting_$rand_number
mv txt/ csv/ xlsx/ reporting_$rand_number 
# make and fill current folders from this run
mkdir txt csv xlsx

# cleanup .txt, .csv and .xlsx files if they exists in directory
count=`ls -1 *.txt 2>/dev/null | wc -l`
if [[ $count != 0 ]]
    then 
    mv *.txt txt
fi 
count=`ls -1 *.csv 2>/dev/null | wc -l`
if [[ $count != 0 ]]
    then 
    mv *.csv csv
fi 
count=`ls -1 *.xlsx 2>/dev/null | wc -l`
if [[ $count != 0 ]]
    then 
    mv *.xlsx xlsx
fi 

cd $bamDir
echo "INFO: Completed cleanup reporting"

echo "INFO: ########### Completed Wochenende_postprocess #############"

