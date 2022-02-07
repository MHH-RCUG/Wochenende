#!/bin/bash
# Automated postprocessing of results from the Wochenende pipeline, with wochenende reporting and haybaler.
# Authors: Colin Davenport, Sophia Poertner

version="0.36, Jan 2022"

#Changelog
#0.36 - add random number to allow concurrent independent fileprep jobs
#0.35 - update raspir and attempt raspir samtools depth parallel speedup
#0.34 - resolve bug with ordering of conda envs with -a option (thanks @irosenboom, @vangreuj )
#0.33 - add -a all option, test refactoring with global scheduler setup 
#0.32 - add command line args
#0.31 - remove --no-plots
#0.30 - remove growth rate temp files
#0.29 - solve growth rate problems
#0.28 - check important env vars are set.
#0.27 - add growth rate estimation by Tom Wehrbein @Leibniz University Hannover
#0.26 - don't save voluminous sambamba depth output to log
#0.25 - add viral read extraction
#0.24 - use bash config.yaml parsing
#0.23 - handle mq20 output files
#0.22 - add heat trees
#0.21 - attempt recovery for second runs to copy data from csv or txt subdirs into haybaler dir
#0.20 - add haybaler heat tree support
#0.19 - update haybaler copying and add double square brackets for bash ifs
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
echo "INFO: Usage: bash wochenende_postprocess.sh -r -h -s -p -g"
echo "INFO: Remember to run this using the haybaler conda environment if available - we attempt to load this in the script"
echo "INFO:  ####### "
echo "INFO:  Usage: Make sure the directories plots/ and reporting/ exist and are filled"
echo "INFO:  eg. run: bash get_wochenende.sh to get the relevant files"
echo "INFO:  ####### "



# Setup conda, directories and SLURM using data parsed from config.yaml
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)
haybaler_dir=$HAYBALER_DIR
wochenende_dir=$WOCHENENDE_DIR
# Set and activate existing conda env
. $CONDA_SH_PATH
conda activate $WOCHENENDE_CONDA_ENV_NAME
# Setup job scheduler
# use SLURM job scheduler (yes, no)
if [[ "${USE_SLURM}" == "yes" && "${USE_CUSTOM_SCHED}" == "yes" ]]; then
    echo "Config warning, both USE_SLURM and USE_CUSTOM_SCHED set. Defaulting to SLURM"
fi
# Setup job scheduler
# use SLURM job scheduler (yes, no)
if [[ "${USE_CUSTOM_SCHED}" == "yes" ]]; then
    #echo USE_CUSTOM_SCHED set"
    scheduler=$CUSTOM_SCHED_CUSTOM_PARAMS
fi
if [[ "${USE_SLURM}" == "yes" ]]; then
    #echo USE_SLURM set"
    scheduler=$SLURM_CUSTOM_PARAMS
fi

# check if env variables could be defined.
if [[ -z "${WOCHENENDE_DIR}" || -z "${HAYBALER_DIR}" ]]; then
    echo "ERROR: WOCHENENDE_DIR or HAYBALER_DIR was not found. Use setup.sh in the Wochenende project to set the directory properly. Exiting! "
    exit 1
fi

# Setup variable sleep duration. Might be useful to set higher for some big projects where the filesystem IO becomes too much.
# The wait command may fail for some SLURM jobs in these cases simply because files have not yet been written in time.
sleeptimer=12
#sleeptimer=120

# get current dir containing Wochenende BAM and bam.txt output
bamDir=$(pwd)

# Cleanup previous results to a directory with a random name which includes a number, calculated here.
rand_number=$RANDOM
# Save output log in directory containing bams and preprocess script
output_log=$bamDir"/postprocess_"$(date +%s)".log"
echo "INFO: output_log: " $output_log


### Check if required directories/files exist, copy if missing ###
if [[ ! -d "reporting" ]] 
then
    echo "INFO: Copying directory reporting, as it was missing! Use get_wochenende.sh to set up Wochenende properly." 
    cp -R $wochenende_dir/reporting .
fi
if [[ ! -d "plots" ]] 
then
    echo "INFO: Copying directory plots, as it was missing!" 
    cp -R $wochenende_dir/plots .
fi
if [[ ! -d "extract" ]] 
then
    echo "INFO: Copying directory extract, as it was missing!" 
    cp -R $wochenende_dir/extract .
fi
if [[ ! -d "growth_rate" ]] 
then
    echo "INFO: Copying directory growth_rate, as it was missing!" 
    cp -R $wochenende_dir/growth_rate .
fi
if [[ ! -f "reporting/ref.tmp" ]] 
then
    echo "INFO: Missing file reporting/ref.tmp , attempting to copy ./ref.tmp to reporting/ref.tmp" 
    cp ref.tmp reporting/ref.tmp
fi


## Set command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -a) runAll="1";  ;;
        -r) runReporting="1";  ;;
        -s) runRaspir="1";  ;;
        -p) runPlotting="1";  ;;
        -g) runGrowth="1";  ;;
        -h) runHaybaler="1";  ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ $runAll == "1" ]] 
then
    runReporting="1"
    runRaspir="1"
    runPlotting="1"
    runGrowth="1"
    runHaybaler="1"
fi
echo "###################################################"
echo "INFO: Selected stages via command line args  (1 means activated)"
echo "INFO: Arguments Run all stages (-a)   : $runAll"
echo "INFO: Arguments Run reporting (-r)    : $runReporting"
echo "INFO: Arguments Run haybaler (-h)     : $runHaybaler"
echo "INFO: Arguments Run Raspir (-s)       : $runRaspir"
echo "INFO: Arguments Run plotting (-p)     : $runPlotting"
echo "INFO: Arguments Run growth rate (-g)  : $runGrowth"

echo "INFO: Starting Wochenende_postprocess" 
echo "INFO: Current directory" $bamDir >>$output_log 2>&1
echo "INFO: Current directory" $bamDir
echo "INFO: Current directory" $bamDir >>$output_log 2>&1
sleep 3


# Run simple filtering and sorting script and multiqc
bash runbatch_metagen_awk_filter.sh
wait

# Check args
if [[ $runHaybaler == "1" ]]; then
    #echo "INFO: Haybaler requires Wochenende reporting"
    runReporting="1"
fi

# Run reporting 
if [[ $runReporting == "1" ]]; then
    echo "INFO: Started Wochenende reporting"
    echo "INFO: Started Wochenende reporting" >>$output_log 2>&1
    cd $bamDir
    cd reporting
    cp ../*.bam.txt .
    bash run_Wochenende_reporting_SLURM.sh >>$output_log 2>&1
    wait
    echo "INFO: Sleeping for "$sleeptimer "to allow writes to complete."
    sleep $sleeptimer
    echo "INFO: Completed Wochenende reporting"
fi

# Run haybaler
if [[ $runHaybaler == "1" ]]; then
    echo "INFO: Start Haybaler"
    echo "INFO: Start Haybaler" >>$output_log 2>&1
    if [[ ! -d "haybaler" ]]
        then
        mkdir haybaler
    fi
    # count files and only copy files that exist to avoid missing files errors
    count_mq20=`ls -1 *mq20.bam*us*.csv 2>/dev/null | wc -l`
    count_mq30=`ls -1 *mq30.bam*us*.csv 2>/dev/null | wc -l`
    count_dup=`ls -1 *dup.bam*us*.csv 2>/dev/null | wc -l`
    count=`ls -1 *.bam*us*.csv 2>/dev/null | wc -l`
    if [[ $count_mq30 != 0 ]]
        then
        cp *mq30.bam*us*.csv haybaler
    elif [[ $count_mq20 != 0 ]]
        then
        cp *mq20.bam*us*.csv haybaler
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
    conda activate $HAYBALER_CONDA_ENV_NAME
    cp $haybaler_dir/*.sh .
    cp $haybaler_dir/*.py .
    cp $haybaler_dir/*.R .
    bash run_haybaler.sh $haybaler_dir >>$output_log 2>&1
    wait
    cp $haybaler_dir/runbatch_heatmaps.sh haybaler_output/ && cp $haybaler_dir/*.R haybaler_output/
    cp $haybaler_dir/*tax* haybaler_output/
    cp $haybaler_dir/*tree* haybaler_output/

    echo "INFO: Attempting to filter results and create heatmaps. Requires R installation." 
    cd haybaler_output
    bash runbatch_heatmaps.sh >>$output_log 2>&1
    echo "INFO: Attempting to add taxonomy. Requires pytaxonkit." 
    bash run_haybaler_tax.sh >>$output_log 2>&1
    echo "INFO: Attempting to create heat-trees. Requires R installation and packages: packages = c("metacoder", "taxa", "dplyr", "tibble", "ggplot2")." 
    bash run_heattrees.sh >>$output_log 2>&1
    cd ..
    cd ..
    echo "INFO: Sleeping for "$sleeptimer "to allow writes to complete."
    sleep $sleeptimer

    echo "INFO: Start csv to xlsx conversion"
    bash runbatch_csv_to_xlsx.sh >>$output_log 2>&1
    wait
    echo "INFO: Sleeping for "$sleeptimer "to allow writes to complete."
    sleep $sleeptimer
    echo "INFO: Completed Haybaler"
fi

# leave haybaler conda env, reactivate Wochenende conda env
conda deactivate
conda activate $WOCHENENDE_CONDA_ENV_NAME


if [[ $runReporting == "1" ]]; then

    echo "INFO: Start cleanup reporting"
    cd $bamDir
    cd reporting
    # create backup, move folders from previous reporting run to a directory (if the txt directory exists already)
    mkdir reporting_$rand_number
    if [[ -d "txt" ]]
        then
        mv txt/ csv/ xlsx/ reporting_$rand_number
    fi

    # make and fill current folders from this run
    mkdir txt csv xlsx

    # cleanup .txt, .csv and .xlsx files if they exist in directory
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
fi

if [[ $runPlotting == "1" ]] 
    then
    echo "INFO: Started Sambamba depth"
    bash runbatch_sambamba_depth.sh >>/dev/null 2>&1
    wait
    echo "INFO: Sleeping for "$sleeptimer "to allow writes to complete."
    sleep $sleeptimer
    bash runbatch_metagen_window_filter.sh >>$output_log 2>&1
    wait
    echo "INFO: Completed Sambamba depth and filtering"
fi


# Growth rate
if [[ $runGrowth == "1" ]]; then
    echo "INFO: Started bacterial growth rate analysis"
    echo "INFO: Started bacterial growth rate analysis" >>$output_log 2>&1
    cd $bamDir
    cd growth_rate/
    echo "INFO: Cleanup original results" >>$output_log 2>&1
    rm -rf fit_results >>/dev/null 2>&1
    rm -rf *_subsamples >>/dev/null 2>&1
    bash runbatch_bed_to_csv.sh  >>$output_log 2>&1 
    wait
    echo "INFO: Sleeping for "$sleeptimer "to allow writes to complete."
    sleep $sleeptimer

    bash run_reproduction_determiner.sh  >>$output_log 2>&1
    wait
    echo "INFO: Sleeping for "$sleeptimer "to allow writes to complete."
    sleep $sleeptimer
    echo "INFO: Files produced by growth rate"   >>$output_log 2>&1
    ls "fit_results/output/*" >>$output_log 2>&1
    cd $bamDir
    echo "INFO: Completed bacterial growth rate analysis, see growth_rate/fit_results/output for results"
fi

# Plots
if [[ $runPlotting == "1" ]]; then
    echo "INFO: Started Wochenende plot"
    echo "INFO: Started Wochenende plot" >>$output_log 2>&1
    cd $bamDir
    cd plots
    cp ../*_window.txt . 
    cp ../*_window.txt.filt.csv .

    bash runbatch_wochenende_plot.sh >>$output_log 2>&1
    #wait
    echo "INFO: Sleeping for "$sleeptimer "to allow writes to complete."
    sleep $sleeptimer
    cd $bamDir
    echo "INFO: Completed Wochenende plot"

    echo "INFO:  - Extracting selected human viral pathogen reads"
    echo "INFO:  - Extracting selected human viral pathogen reads" >>$output_log 2>&1
    cd $bamDir
    bash extract_viral_reads.sh >>$output_log 2>&1
    wait
    echo "INFO: Sleeping for "$sleeptimer "to allow writes to complete."
    sleep $sleeptimer
fi

# raspir
if [[ $runRaspir == "1" ]]; then
    echo "INFO: Run raspir by M. Pust"
    echo "INFO: Run raspir by M. Pust"  >>$output_log 2>&1
    conda activate $HAYBALER_CONDA_ENV_NAME
    cd $bamDir/raspir
    #cleanup
    rm -f *raspir*.csv >/dev/null 2>/dev/null
    echo "INFO: link BAM files in"  >>$output_log 2>&1
    bash batch_create_links.sh  >>$output_log 2>&1
    echo "INFO: Start preparing the files for raspir. Now with SLURM loop"  >>$output_log 2>&1
    # Use random number to create unique fileprep job names and avoid clashes when multiple fileprep jobs are running
    rand=$RANDOM
    for input_bam in `ls *.bam`
        do      
        if [[ "${USE_SLURM}" == "yes" ]]; 
        then
            scheduler=$SLURM_CUSTOM_PARAMS
            # SLURM job scheduler- srun will not work here, need sbatch
            sbatch -J fileprep$rand run_SLURM_file_prep.sh $input_bam >>$output_log 2>&1
        else
            # local job submission
            bash run_SLURM_file_prep.sh $input_bam >>$output_log 2>&1
        fi
    done
    echo "INFO: waiting for raspir file prep jobs to complete"
    srun --dependency=singleton --job-name=fileprep$rand sleep $sleeptimer
    wait
    echo "INFO: Run raspir"  >>$output_log 2>&1
    sbatch run_raspir_SLURM.sh  >>$output_log 2>&1
    echo "INFO: Remove soft linked BAM files"  >>$output_log 2>&1
    bash batch_remove_links.sh  >>$output_log 2>&1
    cd $bamDir
    echo "INFO: Raspir module completed"
fi


echo "INFO: Remember to check the output log for errors at: " $output_log
echo "INFO: ########### Completed Wochenende_postprocess #############"

